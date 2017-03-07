#pragma once

#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <array>
#include <set>
#include <utility>
#include <mutex>
#include <chrono> // time required for temporary filename construction

#include "sequence.h"

//#include "SAM.h" //markus : not needed? 
//#include "sqlite.h" //markus not needed?
#include "pack.h"
#include "ksw.h" // Later we should create our own aligner
#include "configure.h"
#include "container.h"
#include <boost/python.hpp>


extern Configuration xGlobalConfiguration;

/* Some flag positions in the context of MEM_OPTIONS.
 */
#define MEM_F_NO_MULTI  0x10
#define MEM_F_SELF_OVLP 0x40
#define MEM_F_ALL       0x80 
#define MEM_F_SOFTCLIP  0x200

/* Options for aligners.
 * FIX ME: Find a better location for this class than FM_index.h
 */
class MEM_OPTIONS
{
	/* Matrix for 5 residues (A, C, G, T, N).
	 */
	void initScoreMatrix( int iMatchScore, int iMismatchPenalty, int8_t mat[25] )
	{
		int i, j, k;
		for ( i = k = 0; i < 4; ++i )
		{
			for ( j = 0; j < 4; ++j )
				mat[k++] = i == j ? iMatchScore : -iMismatchPenalty;
			mat[k++] = -1; // ambiguous base
		}
		for ( j = 0; j < 5; ++j ) mat[k++] = -1;
	} // method

public:
	/* Alignment parameter
	 */
	int iMatchScore = 2;			// Score for match (normally 1)
	int iMismatchPenalty = 4;		// Mismatch penalty (represented as positive value) (normally 4)
	int o_del = 6;					// Penalty for starting a deletion (normally 6)
	int o_ins = 6;					// Penalty for starting an insertion (normally 6)
	int e_del = 1;					// Penalty for extending a deletion (normally 1)
	int e_ins = 1;					// Penalty for extending an insertion (normally 1)

	// int pen_unpaired;			// phred-scaled penalty for unpaired reads (only used for paired alignments!)
	int pen_clip5 = 5;				// clipping penalty. This score is not deducted from the DP score.
	int	pen_clip3 = 5;					
	int w = 100;					// band width (maximals size of insertions and deletions)
	int zdrop = 100;				// Z-dropoff used in ksw_extend (SW extension)

	int iOutputScoreThreshold = 1;	// output score threshold; only affecting output (value 30 taken from original code)
	int flag = 0;					// see MEM_F_* macros
	unsigned int min_seed_len = 16;	// minimum seed length (the higher the seed length the faster the aligner, but for the price of worse quality)
	uint32_t min_chain_weight = 0;		// chains with a weight equal or below the minimum are dropped (weight is the number of bases of the query represented in a chain)
	int max_chain_extend = 1<<30;	// maximum number of chains that survive the chain filtering and that become part of the extension process
	float split_factor = 1.5;		// split into a seed if MEM is longer than min_seed_len * split_factor
	int split_width = 10;			// split into a seed if its occurrence is smaller than this value
	int max_occ = 500;				// skip a seed if its occurrence is larger than this value
	int max_chain_gap = 10000;      // do not chain a seed if it is max_chain_gap-bp away from the closest seed
	// int n_threads;				// number of threads
	// int chunk_size;				// process chunk_size-bp sequences in a batch
	float mask_level = 0.5;			// regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
	float drop_ratio = 0.5;			// drop a chain if its seed coverage is below drop_ratio times the seed coverage of a better chain overlapping with the small chain
	// float XA_drop_ratio;			// when counting hits for the XA tag, ignore alignments with score < XA_drop_ratio * max_score; only effective for the XA tag
	float mask_level_redun = 0.95f;	// used for determining if a hit (alignment) is considered redundant.

	int mapQ_coef_len = 50;			// Used in the context of the mapping quality computation
	float mapQ_coef_fac;			// Used in the context of the mapping quality computation (initialized in the constructor)

	// int max_ins;					// when estimating insert size distribution, skip pairs with insert longer than this value
	// int max_matesw;				// perform maximally max_matesw rounds of mate-SW for each end
	// int max_hits;				// if there are max_hits or fewer, output them all
	int8_t mat[25];					// scoring matrix; mat[0] == 0 if unset

#if INTERNAL_PARAMETER_PARSING == 0
	// bool copy_comment;		// optimization switch that avoid the copy of annotations.
	// std::string rg_line;

	MEM_OPTIONS()
		// : rg_line()	// initialize rg_line as empty string
	{
		mapQ_coef_fac = (float)log(mapQ_coef_len);
		initScoreMatrix( iMatchScore, iMismatchPenalty, mat );
	} // constructor
#endif
}; // struct

typedef uint64_t bwt64bitCounter;
typedef int64_t t_bwtIndex; // IMPORTANT: We can have -1 in the context of occurrence counting.

#ifndef WORK_WITH_SINGLE_REFERENCE
	#define WORK_WITH_SINGLE_REFERENCE ( 0 )
#endif // !WORK_WITH_SINGLE_REFERENCE

/* The original code worked with calloc. We use STL vectors now.
 */
#define confWORK_WITH_DEPRECTED_SA ( 0 )
#define confWORK_WITH_DEPRECTED_BWT ( 0 )

/* requirement: (OCC_INTERVAL%16 == 0); please DO NOT change this line because some part of the code assume OCC_INTERVAL=0x80
 * TO DO: Put these as constants into class FM_index.
 */
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

/* The below configuration switch is only used if we have confWORK_WITH_DEPRECTED_SA == 1
 */
#define confBWT_USE_NEW_OPERATOR ( 1 )

/* Overloading of the operator "<<" for arrays, so that we get automatic formated output for arrays.
 * From: http://stackoverflow.com/questions/19152178/printing-an-stdarray
 * TO DO: Move me to the section for meta programming.
 */
template <class T, std::size_t N>
std::ostream& operator<< ( std::ostream& rxOstream, const std::array<T, N> &raArray )
{
	std::copy( raArray.cbegin(), raArray.cend(), std::ostream_iterator<T>( rxOstream, " " ) );
	return rxOstream;
} // template function

/* The BWT interval class represents intervals with respect the suffix array!
 * See http://en.wikipedia.org/wiki/FM-index (intervals within the F-L matrix )
 */
class SA_IndexInterval {
public :
	/* x[0] = start position of interval in BWT
	 * x[1] = start position of interval in "complement BWT"
	 * x[2] = size of interval in BWT.
	 */
	t_bwtIndex x[3]; // uint64
	
	/* Query sequence related information:
	 * The unsigned 64 bit value consists of two 32 bit values (high: start point of match in query, low: end point of match in query)
	 */
	uint64_t info;

	/* The length of all matches in the current interval with respect to the query sequence.
	 */
	inline uint32_t sizeOfMatch() const
	{
		return (uint32_t)info - (info >> 32);
	} // method
	
	/* Start index for all matches in the interval with respect to the query sequence.
	 */
	inline uint32_t startIndexOfMatchInQuery() const
	{
		return (info >> 32);
	} // method

//markus
	inline void setMatchInQuery(const uint32_t uiStart, const uint32_t uiEnd)
	{
		info = ((uint64_t)uiStart)<<32 | uiEnd;
		assert(startIndexOfMatchInQuery() == uiStart);
		assert(endIndexOfMatchInQuery() == uiEnd);
	}//function
//end markus

	/* End index for all matches in the interval with respect to the query sequence.
	 * Extract the lower 32 bit of info.
	 */
	inline uint32_t endIndexOfMatchInQuery() const
	{
		return (uint32_t)info;
	} // method
	
	/* Sort order in original code: #define intv_lt(a, b) ((a).info < (b).info)
	 */
	inline bool operator < (const SA_IndexInterval& operand) const
	{
		return (info < operand.info);
	} // method
}; // class ( BWT_Interval )

/* We require the appropriate function in is.c
 * FIX ME: BWT construction should be move to C++ as well.
 */
typedef unsigned char ubyte_t;
extern "C" int is_bwt(ubyte_t *T, int n); // construction of small size BWT
extern std::tuple< uint64_t, uint64_t, std::vector<uint64_t>, std::vector<unsigned int > > bwtLarge( const char *pcFileNamePack ); // construction of large size BWT



/* Creation of FM-indices.
 * FM-indices are central with the alignment process.
 * The original data-structure was part of the BWA-code.
 */
class FM_Index
{
public :
	/* C(), cumulative counts (part of the FM index)
	 * Initialization happens in the constructor.
	 */
	std::array<uint64_t, 6> L2; // uint64_t L2[6]; // 0, A, C, G, T, N

	/* S^{-1}(0), or the primary index of BWT.
	 * BWT-index, position of $
	 */
	t_bwtIndex primary; 


protected :
	typedef int64_t bwtint_t; // changed from uint64_t to int64_t
	
	/* Count-table for 4 values of 2 bit.
	 * Table is initialized by _vInitCountTable
	 */
	uint32_t _aCountTable[256];

#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
	/* Burrow-Wheeler transform stored in 4 byte blocks.
	 * (Pointer to first element of array on heap.)
	 */
	uint32_t *bwt; // BWT

	/* Size of the bwt array expressed in 32 bit words.
	 * (Assertion: bwt_size = (seq_len + 15) >> 4)
	 */
	uint64_t bwt_size;
#else
	std::vector<uint32_t> bwt;
#endif

	/* Size of the reference sequence that was used for BWT construction.
	 * Length of the reference sequence (forward plus reverse strand)
	 */
	uint64_t uiRefSeqLength;

	/* suffix array
	 */
	int sa_intv;	// interval size (must be power of 2 and is normally 32, because of 64 bit / 2 == 32 bit)

#if ( confWORK_WITH_DEPRECTED_SA == 1 )
	bwtint_t n_sa;	// Number of intervals of the suffix array
	bwtint_t *sa;	// Suffix array itself
#else
	std::vector<bwtint_t> sa;
#endif

	/* Core BWT construction for short sequences.
	 * Initializes the BWT using the given nucleotide sequence.
	 * We have a serious problem with ambiguous bases over here, because the (compressed) BWT can't represent them.
	 * WARNING: Do not pass sequences comprising ambiguous symbols (e.g. symbol 'N').
	 */
	void bwt_pac2bwt_step1( const NucleotideSequence &fn_pac_arg, 
							bool bIncludeReverse // In the context of the BWT construction generate and include data for the reverse strand as well.
						  )
	{	
		NucleotideSequence fn_pac;
		fn_pac.vAppend( fn_pac_arg.fullSequenceAsSlice() );
		
		/* Size of the reference sequence
		 */
		uiRefSeqLength = fn_pac.length(); // bwt->seq_len = bwa_seq_len( fn_pac );
		
		/* The sequence length is doubled if we have the reverse sequence in th BWT.
		 */
		uiRefSeqLength = bIncludeReverse ? uiRefSeqLength * 2 : uiRefSeqLength;

		/* Buffer for BWT construction. The BWT will be finally inside the buffer.
		 */
		std::unique_ptr<uint8_t[]> buf( new uint8_t[uiRefSeqLength + 1] ); // original code: buf = (ubyte_t*)calloc( bwt->seq_len + 1, 1 );

		/* Initialize the buffer with the reference sequence and prepare the cumulative count.
		 */
		for ( uint64_t i = 0; i < ( bIncludeReverse ? uiRefSeqLength / 2 : uiRefSeqLength ); ++i )
		{
			buf[i] = fn_pac[i]; // buf2[i >> 2] >> ((3 - (i & 3)) << 1) & 3;
			++L2[1 + buf[i]];
		} // for

		if ( bIncludeReverse )
		{
			/* Old form of computing the reverse strand...
			 */
			for ( uint64_t i = 0; i < uiRefSeqLength / 2; ++i )
			{
				buf[uiRefSeqLength - i - 1] = 3 - fn_pac[i]; // buf2[i >> 2] >> ((3 - (i & 3)) << 1) & 3;
				++L2[1 + buf[uiRefSeqLength - i - 1]];
			} // for

			/* Give the reference sequence bigger size
			 */
			fn_pac.vAppend( fn_pac.fullSequenceAsSlice().makeComplementSequenceUsingReverseOrder()->fullSequenceAsSlice() );
		} // if
		
		/* Complete cumulative count
		 */
		for ( int i = 2; i <= 4; ++i )
		{
			L2[i] += L2[i - 1];
		} // for
		
		/* Burrows-Wheeler Transform
		 * Here we call external code (in is.c) for BWT construction. 
		 * The BWT construction happens in-place in buf. (This is why we allocate one byte additionally.)
		 */
	#ifdef _DIVBWT
		primary = divbwt( buf, buf, 0, uiRefSeqLength );
	#else
		primary = is_bwt( buf.get(), (int)uiRefSeqLength );
	#endif

#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
		/* Size of the bwt array expressed in 32 bit words.
		 */
		bwt_size = ( uiRefSeqLength + 15 ) >> 4; // ( seq_len + 15 ) / 16
		
		/* Allocates the memory for the compressed BWT.
		 * We should better use the new operator here ...
		 */
	#if ( defBWT_USE_NEW_OPERATOR == 1 )
		bwt = new uint32_t[bwt_size](); // zero initialized -> http://stackoverflow.com/questions/7546620/operator-new-initializes-memory-to-zero
	#else
		bwt = (uint32_t*)calloc( bwt_size, 4 );
	#endif
#else
		/* Compute the expected size of the bwt expressed in 32 bit blocks and resize the bwt vector appropriately.
		 */
		auto uiExpectedBWTSize = ( uiRefSeqLength + 15 ) >> 4;
		bwt.resize( uiExpectedBWTSize );
#endif
		/* Copy bwt from buffer into vector.
		 * FIX ME: The construction programs should already work with the buffer.
		 */
		for ( uint64_t i = 0; i < uiRefSeqLength; ++i )
		{	/* Byte to packed transformation. 
			 */
			bwt[i >> 4] |= buf[i] << ((15 - (i & 15)) << 1);
		} // for
	} // method

	/* Retrieve character at position k from the $-removed packed BWT without occurrence counter.
	 * Used in the context of step 2 (bwt_bwtupdate_core_step2) for retrieving characters.
	 * REMARK: Only correct in the case of OCC_INTERVAL==0x80
	 */
	#define bwt_B00(k) (bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

	/* Injection of counting blocks into the BWT.
	 * Inserts the required space for 4 uint_64 counters each 8 uint_32 blocks and initializes the counters according to the counting.
	 * FIX ME: If we get exception over here, we have a memory leak -> use unique poiters are vectors
	 */
	void bwt_bwtupdate_core_step2()
	{
		bwtint_t i, k, c[4], n_occ;

		/* Number of occurrence counter intervals.
		 * Normally OCC_INTERVAL is 0x80 = 128, so we count the symbols for blocks of 128 symbols.
		 */
		n_occ = (uiRefSeqLength + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;

#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
		/* New size of the BWT including the counting blocks. (Expressed in words of 32 bits.)
		 */
		bwt_size += n_occ * sizeof(bwtint_t);

		/* We reserve fresh memory for the enlarged BWT; will be the new bwt, now with the counting blocks sizeof(uint32_t)
		 */
	#if ( defBWT_USE_NEW_OPERATOR == 1 )
		uint32_t *buf = new uint32_t[bwt_size](); // zero initialized -> http://stackoverflow.com/questions/7546620/operator-new-initializes-memory-to-zero
	#else
		uint32_t *buf = (uint32_t*)calloc( bwt_size, 4 );
	#endif
#else
		/* Create a temporary vector for storage of the BWT with injected Occ-counting blocks.
		 */
		auto uiOccInjectedBWTExpectedSize = bwt.size() + n_occ * sizeof(bwtint_t);
		std::vector<uint32_t> xVectorWithOccInjections( uiOccInjectedBWTExpectedSize );
#endif

		c[0] = c[1] = c[2] = c[3] = 0; // temporary uint_64 counter for A, C, G, T

		/* We iterate over the existing BWT and insert each 128th element.
		 * k counts in 16 nucleotides.
		 * i iterates over the nucleotides.
		 * FIX ME: In the context of STL vectors the memcpy is not that nice anymore ...
		 */
		for ( i = k = 0; i < static_cast<bwtint_t>(uiRefSeqLength); ++i )
		{
			if ( i % OCC_INTERVAL == 0 )
			{
#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
				memcpy( buf + k, c, sizeof(bwtint_t) * 4 ); // set the 4 counter initially to the current value of c
#else
				memcpy( &xVectorWithOccInjections[0] + k, c, sizeof(bwtint_t) * 4 ); // set the 4 counter initially to the current value of c
#endif

				/* Increment k according to the required size. 4 * sizeof(bwtint_t) many bytes required.
				 * k counts in 4 bytes steps. ( in fact: sizeof(bwtint_t)=4*(sizeof(bwtint_t)/4) )
				 */
				k += sizeof(bwtint_t);
			} // if

			/* i iterates over the nucleotides, so each 16 nucleotides we have a uint32_t for copying.
			 */
			if ( i % 16 == 0 )
			{
#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
				buf[k++] = bwt[i / 16]; // 16 == sizeof(uint32_t)/2
#else
				xVectorWithOccInjections[k++] = bwt[i / 16]; // 16 == sizeof(uint32_t)/2
#endif
			} // if

			/* We extract the nucleotide (2 bits) on position i out of corresponding 16 nucleotide block and increment the correct counter for the symbol.
			 * ++c[ bwt->bwt[i >> 4] >>((~i & 0xf)<<1)&3 ];
			 */
			//// std::cout << bwt_B00( i ); (here we can pick up the BWT elements in the context of debugging before inserting the Counting Blocks)
			++c[bwt_B00( i )];
		} // for
		////std::cout << "\n";

		/* Save the last counter-block (element).
		 * Check whether everything is fine.
		 */
#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
		memcpy( buf + k, c, sizeof(bwtint_t)* 4 );
		assert( k + sizeof(bwtint_t) == bwt_size ); //  "inconsistent bwt_size" 

		// update bwt(free the old one, assign the new one)
	#if ( defBWT_USE_NEW_OPERATOR == 1 )
		delete[] bwt;
	#else
		free( bwt );
	#endif
		bwt = buf;
#else
		memcpy( &xVectorWithOccInjections[0] + k, c, sizeof( bwtint_t ) * 4 );
		assert( k + sizeof( bwtint_t ) == xVectorWithOccInjections.size() ); //  "inconsistent bwt_size" 

		/* Move the bwt with injections to the bwt of the FM-index.
		 */
		bwt = std::move( xVectorWithOccInjections );
#endif
	} // method

	/* Retrieve character at position k from the $-removed packed BWT, which comprises occurrence counter.
	 * (Note that bwt_t::bwt is not exactly the BWT string and therefore this macro is called bwt_B0 instead of bwt_B.)
	 * REMARK: Only correct in the case of OCC_INTERVAL==0x80
	 */
	#define bwt_bwt(k) (bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)])
	#define bwt_B0(k) (bwt_bwt(k)>>((~(k)&0xf)<<1)&3)
	
	/* The function Occ(c, k) is the number of occurrences of character c in the prefix L[1..k]. 
	 * Ferragina and Manzini showed[1] that it is possible to compute Occ(c, k) in constant time.
	 * ( Paolo Ferragina and Giovanni Manzini (2000). "Opportunistic Data Structures with Applications". 
	 *   Proceedings of the 41st Annual Symposium on Foundations of Computer Science. p.390. )
	 * y contains 32 nucleotides encoded in 2 bits each. c is a single character A/C/G/T = 0/1/2/3
	 * Returned value is the number of occurrences of this character in y.
	 */
	inline int __occ_aux( uint64_t y, int c )
	{
		// reduce nucleotide counting to bits counting
		y =	  ((c & 2)? y : ~y) >> 1 
			& ((c & 1)? y : ~y) 
			& 0x5555555555555555ull;
		
		// count the number of 1s in y
		y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
		
		return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
	} // method
	
	/* Counts the number of occurrences of the character c in the BWT until position k.
	 * bwt_occ counts for a single character, while bwt_occ counts for all 4 characters (A/C/G/T) simultaneously.
	 * k should be within [0..seqlen], the value of k shall be not equal to primary (if you want to get the value at primary position deliver (bwtint_t)(-1)
	 */
	inline bwtint_t bwt_occ( t_bwtIndex k, ubyte_t c )
	{
		bwtint_t n;
		uint32_t *p, *end;
	
		if ( k == (t_bwtIndex) uiRefSeqLength )
		{
			return L2[c+1] - L2[c];
		} // if
		
		/* assertion: k != primary 
		 */
		if ( k == (bwtint_t)(-1) )
		{
			return 0;
		} // if
		
		k -= (k >= primary); // because $ is not in bwt
	
		/* retrieve Occ at k/OCC_INTERVAL
		 */
		n = ((bwtint_t*)(p = bwt_occ_intv( k )))[c];
		p += sizeof(bwtint_t); // jump to the start of the first BWT cell
	
		/* calculate Occ up to the last k/32
		 */
		end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
		for (; p < end; p += 2)
		{
			n += __occ_aux( (uint64_t)p[0] << 32 | p[1], c);
		} // for
	
		// calculate Occ
		n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
		if (c == 0)
		{
			n -= ~k&31; // corrected for the masked bits
		} // if
	
		return n;
	} // method

	/* Computes the inverse CSA. (Back to front projection in BWT)
	 * @param k Position in BWT
	 * @return Position of the character at k in F
	 */
	inline bwtint_t bwt_invPsi( bwtint_t k )
	{	/* Adapt, whether we are before or after the primary.
		 */
		bwtint_t x = k - (k > primary);
		
		/* Character in BWT at primary adjusted position
		 */
		uint8_t c = bwt_B0( x );

		/* bwt_occ will adjust k one again!
		 */
		x = L2[c] + bwt_occ(k, c);
		
		return (k == primary) ? 0 : x;
	} // method
	
	/* Small helper function used in the context of some assertion.
	 */
	bool ispowerof2( unsigned int x )
	{
		return x && !(x & (x - 1));
	} // method

	/* Creation of the suffix array (SA is organized on the foundation of 64 bit counters)
	 * @param intv Interval size for suffix array (must be a power of 2)
	 */
	void bwt_cal_sa_step3( unsigned int intv )
	{
		assert( ispowerof2( intv ) );	//, "SA sample interval is not a power of 2.");
#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
		assert( bwt );					//, "bwt_t::bwt is not initialized.");
#endif
	
		// if (bwt->sa) free(bwt->sa);
		sa_intv = intv;

#if ( confWORK_WITH_DEPRECTED_SA == 1 )
		/* Number of intervals of the suffix array
		 */
		n_sa = (uiRefSeqLength + intv) / intv;
		
		/* reserve memory for the suffix array
		 */
	#if ( defBWT_USE_NEW_OPERATOR == 1 )
			this->sa = new bwtint_t[n_sa](); // zero initialized -> http://stackoverflow.com/questions/7546620/operator-new-initializes-memory-to-zero
	#else
			this->sa = (bwtint_t*)calloc(n_sa, sizeof(bwtint_t));
	#endif
#else
		/* FIX ME: resize initializes all elements with zero. Better would be reserve and later vector insertion
		 */
		size_t uiRequiredSAsize = (uiRefSeqLength + intv) / intv;
		sa.resize( uiRequiredSAsize );
#endif
		/* calculate SA values
		 */
		bwtint_t isa = 0; // "jumping" index position
		bwtint_t sa = uiRefSeqLength; // S(isa) = sa (sa at jumping index position); counts down to 1
		
		/* Implements an iteration scheme that counts over all suffix array elements.
		 * Basic idea BWT based suffix array reconstruction.
		 * Important observation: The suffix array can be always reconstructed using the BWT, because indirectly the BWT contains all info.
		 */
		for ( bwtint_t i = 0; i < static_cast<bwtint_t>(uiRefSeqLength); ++i ) 
		{
			/* Do we have an interval position?
			 */
			if (isa % intv == 0)
			{
				this->sa[isa/intv] = sa;
			} // if
			--sa;
			
			/* Core iteration step
			 */
			isa = bwt_invPsi( isa );
		} // for
		
		/* Do not forget the last block, if necessary.
		 */
		if (isa % intv == 0) 
		{
			this->sa[isa/intv] = sa;
		} // if
		
		this->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
	} // method

	/* BWT construction for packs. Is is possible to choose among two different algorithms.
	 * uiAlgorithmSelection : 0 manually selected algorithm for small inputs
	 *						  1 manually selected algorithm for large inputs
	 *						  2 automatic selection on foundation of input size
	 */
	void build_FM_Index( const BWACompatiblePackedNucleotideSequencesCollection &rxSequenceCollection, // the pack for which we compute a BWT
						 unsigned int uiAlgorithmSelection = 2 // 2 -> automatic algorithm selection
					   ) 
	{

		if ( uiAlgorithmSelection > 1) 
		{
			uiAlgorithmSelection = rxSequenceCollection.uiUnpackedSizeForwardPlusReverse() < 50000000 ? 0 : 1; // automatic selection depending on size of pack
		} // if
		
		/* Step 1: Create the basic BWT.
		 */
		if ( uiAlgorithmSelection == 0 )
		{	/* For small packs we transform the pack into a single sequence with reverse strand and apply the algorithm for small inputs.
			 */
			NucleotideSequence xSequence;
			rxSequenceCollection.vColletionWithoutReverseStrandAsNucleotideSequence( xSequence ); // unpack the pack into a single nucleotide sequence
			bwt_pac2bwt_step1( xSequence, true ); // construct with build in function
		} // if
		else
		{	/* In this case we build the BWT using the BWA C-code for large inputs. 
			 * For delivering the pack to the BWA code we have to rely on an external file currently.
			 * FIX ME: The current solution is quite inefficient with respect to memory aspects.
			 */
			uiRefSeqLength = rxSequenceCollection.uiUnpackedSizeForwardPlusReverse(); // set seq_len, so that it represents the size of forward plus reverse strand
#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
			bwt_size = (uiRefSeqLength + 15) >> 4; // bwt is organized in blocks of 32 bit.
#endif
			
			/* Check existence / create directory for storage of temporary data.
			 */
			boost::filesystem::create_directories( xGlobalConfiguration.getTemporaryDataFolder() ); 

			/* Create temporary filename for pack export.
			 */
			const auto xTempFileName =	   xGlobalConfiguration.getTemporaryDataFolder() 
										/= (  std::to_string( std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() ) ) // current time as string
											 .append("-")
											 .append(std::to_string( rxSequenceCollection.uiUnpackedSizeForwardStrand )) // pack size as string
										   ); // temporary filename construction
			const auto sTempFileNameWithSuffix = std::string( xTempFileName.string() ).append(".pac"); // pac file construction adds a suffix
			
			/* Store the dynamically created pack on the filesystem.
			 * The suffix .pac will be automatically added.
			 */
			rxSequenceCollection.vCreateAndStorePackForBWTProcessing( xTempFileName );
			
			/* Start the BWT construction.
			 */
			auto xBWTAsTriple = bwtLarge( sTempFileNameWithSuffix.c_str() ); // construct the BWT using the old BWA code. 

#if ( confWORK_WITH_DEPRECTED_BWT == 1 )			
			/* Copy the externally created BWT into local class attributes.
			 */
			const auto &rxBWTAsVector = std::get<3>( xBWTAsTriple ); // get the BWT as vector out of the tuple
			const auto uiBwtSize = rxBWTAsVector.size(); // get the size of the BWT

	#if ( defBWT_USE_NEW_OPERATOR == 1 )	
			bwt = new std::remove_reference<decltype(*bwt)>::type[uiBwtSize]; // allocate memory for BWT storage
			memcpy( bwt, &rxBWTAsVector[0], uiBwtSize * sizeof( decltype(*bwt) ) ); // copy BWT itself - dirty, but working ....
			//// std::copy( rxBWTAsVector.begin(), rxBWTAsVector.end(), bwt ); // would be better than above statement, but throws a warning.
	#else
			#error "Copy code missing"
	#endif	
			this->bwt_size = std::get<0>( xBWTAsTriple );
#else
			/* Move the externally constructed BWT to this BWT-object.
			 */
			bwt = std::move( std::get<3>( xBWTAsTriple ) );
#endif
			this->primary = std::get<1>( xBWTAsTriple ); // copy the primary
			std::copy_n( std::get<2>( xBWTAsTriple ).begin(), 4, L2.begin() + 1); // copy the 4 cumulative counters
			
			/* Remove the temporary file, because it is not needed any longer.
			 */
			boost::filesystem::remove( sTempFileNameWithSuffix ); 
		} // else

		/* Step 2 and step 3 of FM-index creation
		 */
		vPostProcessBWTAndCreateSA();
	} // method

	/* Builds up a FM_Index for the given input sequence. 
	 * REMARK: The sequence should be free of ambiguous bases (character N).
	 */
	void build_FM_Index( const NucleotideSequence &fn_pac )
	{
		/* Construction of core BWT.
		 */
		bwt_pac2bwt_step1( fn_pac, true );

		/* Step 2 and step 3 of FM-index creation
		 */
		vPostProcessBWTAndCreateSA();
	} // method

	/* Step 2 and step 3 of FM-index creation.
	 */
	void vPostProcessBWTAndCreateSA( void )
	{
		/* Insertion of occurrence data into outcome of step 1.
		 */
		bwt_bwtupdate_core_step2();

		/* Computation of the suffix index cache. (Requires outcome of step 2)
		 * We store the sa-value in intervals of 32 values.
		 */
		bwt_cal_sa_step3( 32 );
	} // method

	/* Initializes the array _aCountTable.
	 */
	void _vInitCountTable()
	{
		for ( unsigned int i = 0; i != 256; ++i )
		{
			uint32_t x = 0; // 4 byte x[3],x[2],x[1],x[0]

			/* Count the occurrences of the 4 2-bit symbols in i
			 */
			for ( unsigned int j = 0; j != 4; ++j )
			{
				x |= (
						 ((i      & 3) == j)
					   + ((i >> 2 & 3) == j)
					   + ((i >> 4 & 3) == j)
					   + ( i >> 6      == j)
					 ) << (j << 3); // (j << 3) is 0, 8, 16, 24
			} // for
			
			/* x[c] = number of c in i, where i is decomposed in groups of 2 bit and c in {0, 1, 2, 3}
			 * Example: i=$FF = 11111111b, so x[0] = 0, x[1] = 0, x[2] = 0, x[3] = 4
			 */
			_aCountTable[i] = x;
		} // for
	} // method

	/* Counts A,C,G,T (00, 01, 10, 11) in ux16symbols .
	 * ux16symbols represents 16 nucleotides encoded using 2 bit each.
	 * The returned uint32_t represents 4 unsigned 8 bit values.
	 */
	inline uint32_t __occ_aux4( uint32_t ux16symbols )
	{
		return   _aCountTable[(ux16symbols)       & 0xff] 
			   + _aCountTable[(ux16symbols) >> 8  & 0xff]	// parallel add	of 4 bytes
			   + _aCountTable[(ux16symbols) >> 16 & 0xff]	// parallel add of 4 bytes
			   + _aCountTable[(ux16symbols) >> 24];			// parallel add of 4 bytes
	} // method

	/* The following two lines are ONLY correct when OCC_INTERVAL==0x80 (OCC_INTERVAL == 128, 2^7)
	 * Delivers pointer to the first 32 bit-block of interval.
	 * FIX ME: In the context of the BWT as vector this is not nice at all. (Work with references over here)
	 */
	inline uint32_t* bwt_occ_intv( t_bwtIndex k )
	{
#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
		return bwt + ( ( k >> 7 ) << 4 ); // (k / 128) * 16  (we take 16, because we need twice the size)
#else
		return &bwt[0] + ( ( k >> 7 ) << 4 ); // (k / 128) * 16  (we take 16, because we need twice the size)
#endif
	} // method

	/* k (input) is the position of some nucleotide within the BWT.
	 * Represents count function within the BWT.
	 * cnt[4] (output) delivers the number of A's, C's, G's and T's that are in front of nucleotide number k.
	 * Special case k = -1: initialize cnt[4] with 0.
	 * See http://en.wikipedia.org/wiki/FM-index for more info. (However, counting for k seems to start with 0 instead of 1!)
	 */
	inline void bwt_occ4( t_bwtIndex k,				// position within the BWT (input)
						  bwt64bitCounter cnt[4]	// output ( 4 counter )
						)
 	{
		uint64_t x; // originally 64 bit, but should work with a 32 bit value as well.
		uint32_t *p, tmp, *end; // pointer into bwt.
		
		if ( k == (t_bwtIndex)(-1) )
		{
			std::memset( cnt, 0, 4 * sizeof(bwt64bitCounter) ); // TO DO: Use C++ sizeof for arrays
			return;
		} // if
		
		/* Adjust k, because $ is not in bwt
		 */
		k -= (k >= primary); 
		
		/* Step 1.
		 * Retrieve Occ at k/OCC_INTERVAL. (These are accumulated data, precomputed in the context of the construction process)
		 * p points to a counter block start.
		 */
		p = bwt_occ_intv( k );
		std::memcpy( cnt, p, 4 * sizeof(bwt64bitCounter) ); // TO DO: Better C++ sizeof for arrays

		//// std::cout << "bwt_size : " << bwt_size << " k: " << k << " " << (k < bwt_size) << "\n";
		//// std::cout << "cnt after copy" << cnt[0] << " " << cnt[1] << " " << cnt[2] << " " << cnt[3] << "\n"; 
		
		/* Step 2. Jump (4(==sizeof(uint_32)) * 8 bytes(==sizeof(uint_64))) to the start of the first BWT cell (nucleotide) (see bwt_bwtupdate_core)
		 * This line is dirty, works only if the counter a 8 byte values.
		 */
		p += sizeof(bwt64bitCounter);

		/* Pointer to the end (with respect to k) within the group
		 */
		end = p + ((k >> 4) - ((k & ~OCC_INTV_MASK) >> 4)); // this is the end point of the following loop
		
		/* Add up nucleotides within the blocks up to k, using the parallel addition approach.
		 * At the end, we have in x the counting for A, C, G, T 
		 */
		for ( x = 0; p < end; ++p )
		{
			x += __occ_aux4( *p );
		} // for
		
		/* Inspect the 16 nucleotide block comprising nucleotide number k.
		 */
		tmp = *p & ~((1U << ((~k & 15) << 1)) - 1);
		x += __occ_aux4( tmp ) - (~k & 15);
		
		/* Add up the element counted for the final 128 nucleotides block
		 */
		cnt[0] += x       & 0xff; // A
		cnt[1] += x >> 8  & 0xff; // C
		cnt[2] += x >> 16 & 0xff; // G
		cnt[3] += x >> 24;        // T

		//// std::cout << "GOT for k:" << k << " " << cnt[0] << " " << cnt[1] << " " << cnt[2] << " " << cnt[3] << "\n"; 
	} // method

	/* Writes the BWT of the FM-Index to the stream given as argument.
	 * WARNING: The given stream should be opened with flag ios::binary.
	 */
	void vSaveBWT( std::ostream &rxOutputStream )
	{	/* Write the primary, 4 elements of L2 array as well as the BWT itself to the file system.
		 */
		rxOutputStream.write( reinterpret_cast<const char *>( &primary ), sizeof( primary ) ); // BWA code: err_fwrite( &bwt->primary, sizeof(bwtint_t), 1, fp );
		rxOutputStream.write( reinterpret_cast<const char *>( &L2[1] ), sizeof( L2[0] ) * 4 ); //  BWA code: err_fwrite( bwt->L2 + 1, sizeof(bwtint_t), 4, fp );
#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
		rxOutputStream.write( reinterpret_cast<const char *>( bwt ), sizeof( bwt[0] ) * bwt_size ); // BWA code: err_fwrite( bwt->bwt, 4, bwt->bwt_size, fp );
#else
		/* Write the vector to the file system.
		 */
		rxOutputStream.write( reinterpret_cast<const char *>( &bwt[0] ), sizeof( bwt[0] ) * bwt.size() );
#endif
		rxOutputStream.flush();
	} // method

	/* Writes the Suffix Array to the given stream.
	 */
	void vSaveSuffixArray( std::ostream &rxOutputStream )
	{
		rxOutputStream.write( reinterpret_cast<const char *>( &primary ), sizeof( primary ) ); // BWA code: err_fwrite( &bwt->primary, sizeof(bwtint_t), 1, fp );
		rxOutputStream.write( reinterpret_cast<const char *>( &L2[1] ), sizeof( L2[0] ) * 4 ); //  BWA code: err_fwrite( bwt->L2 + 1, sizeof(bwtint_t), 4, fp );
		rxOutputStream.write( reinterpret_cast<const char *>( &sa_intv ), sizeof( sa_intv ) ); // BWA code: err_fwrite( &bwt->sa_intv, sizeof(bwtint_t), 1, fp );
		rxOutputStream.write( reinterpret_cast<const char *>( &uiRefSeqLength ), sizeof( uiRefSeqLength ) ); // BWA code: err_fwrite( &bwt->seq_len, sizeof(bwtint_t), 1, fp );
#if ( confWORK_WITH_DEPRECTED_SA == 1 )
		rxOutputStream.write( reinterpret_cast<const char *>( sa + 1 ), n_sa - 1 ); // BWA code: err_fwrite( bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp );
#else
		/* Write the suffix array to an file. rxOutputStream must be in binary mode!
		 * Taken from: http://stackoverflow.com/questions/12372531/reading-and-writing-a-stdvector-into-a-file-correctly
		 */
		rxOutputStream.write( reinterpret_cast<const char *>( &sa[1] ), (sa.size() - 1) * sizeof(bwtint_t) );
		//// The following apporach does not work: Why? (we have no iterator for uint_64) std::copy( sa.begin() + 1, sa.end(), std::ostreambuf_iterator<char>( rxOutputStream ) );
#endif
		rxOutputStream.flush();
	} // method

	/* Reads the BWT from the given stream.
	 * IMPROVEMENT: The original BWA design for BWT storage is quite crappy. Store the size of the BWT as part of the data.
	 */
	void vRestoreBWT( std::ifstream &rxInputStream )
	{
		/* Get the file size.
		 * FIX ME: Move this code to a function of its own. 
		 */
		rxInputStream.seekg( 0, std::ios::end ); // go to the end of stream
		auto xFileEndPos =  rxInputStream.tellg(); // read position
		rxInputStream.seekg( 0, std::ios::beg ); // go to the begin of stream
		auto xFileStartPos =  rxInputStream.tellg(); // read position
		auto uiStreamLength = ( xFileEndPos - xFileStartPos ); // get the size of the stream

#if ( confWORK_WITH_DEPRECTED_SA == 1 )		
		bwt_size = ( uiStreamLength - ( sizeof( primary ) + ( sizeof( L2[0] ) * 4 ) ) ) / sizeof( bwt[0] ); // compute the bwt size using the info about the stream size.
#endif

		if ( rxInputStream.fail() ) // check whether we could successfully open the stream
		{
			throw std::runtime_error( "Opening of BWT of FM index failed: cannot seek." );
		} // if

		rxInputStream.seekg( 0, std::ios::beg ); // go back to the beginning of the stream.
		rxInputStream.read( reinterpret_cast<char *>( &primary ), sizeof( primary ) ); // get the primary back into the data structure.
		rxInputStream.read( reinterpret_cast<char *>( &L2[1] ), sizeof( L2[0] ) * 4 ); // get the L2 values back from the filesystem

#if ( confWORK_WITH_DEPRECTED_SA == 1 )
		/* Allocate memory and read the core BWT.
		 */
	#if ( defBWT_USE_NEW_OPERATOR == 1 )
		bwt = new std::remove_reference<decltype(*bwt)>::type[bwt_size]; //    calloc( bwt_size, sizeof( bwt[0] ) ) );
	#else
		bwt = reinterpret_cast<decltype(bwt)>( calloc( bwt_size, sizeof( bwt[0] ) ) ); // FIX ME: Move away from the C approach to C++ (we should use a unique pointer over here)
	#endif		
		rxInputStream.read( reinterpret_cast<char *>( bwt ), sizeof( bwt[0] ) * bwt_size ); // get the BWT from the filesystem
#else
		/* Compute the bwt size using the info about the stream size.
		 */
		size_t uiExpectedBWTSize = ( uiStreamLength - ( sizeof( primary ) + ( sizeof( L2[0] ) * 4 ) ) ) / sizeof( bwt[0] ); 
		
		/* This resize is quite expensive, because all elements are initialized with 0.
		 */
		bwt.resize( uiExpectedBWTSize );

		/* Restore vector by reading from file system.
		 */
		rxInputStream.read( reinterpret_cast<char *>( &bwt[0] ), sizeof( bwt[0] ) * uiExpectedBWTSize ); 

		if ( rxInputStream.fail() )
		{	/* We should have no bad stream over here, or the input file did not comprise the expected number of bytes.
			 */
			throw std::runtime_error( std::string( "Unexpected fail after reading BWT from stream. ") );
		} // if
#endif
		
		uiRefSeqLength = L2[4]; // restore seq_len;
		L2[0] = 0; // clean L2[0]
	} // method

	/* Restore suffix array by using the given input stream.
	 */
	void vRestoreSuffixArray( std::istream &rxInputStream )
	{
		decltype( primary ) uiPrimaryReferenceValue;
		decltype( uiRefSeqLength ) uiSeqLenReferenceValue;

		rxInputStream.read( reinterpret_cast<char *>( &uiPrimaryReferenceValue ), sizeof( primary ) ); // get the primary back into the data structure.
		if ( primary != uiPrimaryReferenceValue )
		{
			throw std::runtime_error( "BWT and suffix array have different primary." );
		} // if

		char aSkipped[ sizeof( L2[0] ) * 4 ]; // small buffer for dummy reading
		rxInputStream.read( reinterpret_cast<char *>( aSkipped ), sizeof( L2[0] ) * 4 ); // get the L2 values but ignore them.
		rxInputStream.read( reinterpret_cast<char *>( &sa_intv ), sizeof( sa_intv ) ); // read suffix array interval size.
		rxInputStream.read( reinterpret_cast<char *>( &uiSeqLenReferenceValue ), sizeof( uiSeqLenReferenceValue ) ); // read the sequence length for verification purposes
		if ( this->uiRefSeqLength != uiSeqLenReferenceValue )
		{
			throw std::runtime_error( "SA-BWT inconsistency: suffix array has non matching sequence length stored." );
		} // if
		
#if ( confWORK_WITH_DEPRECTED_SA == 1 )		
		n_sa = (uiRefSeqLength + sa_intv) / sa_intv; // compute the size of the suffix array on foundation of sequence length and interval distance
	#if ( defBWT_USE_NEW_OPERATOR == 1 )
			sa = new std::remove_reference<decltype(*sa)>::type[n_sa]; //    calloc( bwt_size, sizeof( bwt[0] ) ) );
	#else
			sa = reinterpret_cast<decltype(sa)>( calloc( n_sa, sizeof( sa[0] ) ) ); // FIX ME: Let us go away from the C-style memory allocation.
	#endif			
		sa[0] = -1; // set the first value of the suffix array to its special value.
		rxInputStream.read( reinterpret_cast<char *>( sa + 1 ), sizeof( sa[0] ) * (n_sa - 1) ); // read all suffix array values except the first one.
#else
		/* Take from: http://stackoverflow.com/questions/15138353/reading-the-binary-file-into-the-vector-of-unsigned-chars
		 * VERY IMPORTANT: Stop eating new lines in binary mode. (required for the back_inserter.)
		 */
		rxInputStream.unsetf( std::ios::skipws );
		
		/* Compute the expected suffix array size and reserve appropriate space with the vector.
		 */
		size_t uiExpectedSASize = (uiRefSeqLength + sa_intv) / sa_intv;
		
		/* This resize is quite expensive, because all elements are initialized with 0.
		 */
		sa.resize( uiExpectedSASize );

		/* Read the complete remaining file data 
		 */
#if 1	// SOLUTION 1
		sa[0] = -1; // special marker value ($FFFFFFFFFFFFFFFF with uint_64t)

		/* Read all suffix array values except the first one, which is a always a special one
		 */
		rxInputStream.read( reinterpret_cast<char *>( &sa[1] ), sizeof( uint64_t ) * (uiExpectedSASize - 1) ); 
		
		if ( rxInputStream.fail() )
		{	/* We should have no bad stream over here, or the input file did not comprise the expected number of bytes.
			 */
			throw std::runtime_error( std::string( "Unexpected bad after reading suffix array from stream. ") );
		} // if
#else	// SOLUTION 2 (FIX ME: Get a solution working the does not require the expensive resizing of the vector.)
		sa.push_back( -1 ); // set sa[0] to -1
		sa.insert( sa.end(), std::istream_iterator<uint64_t>( rxInputStream ), std::istream_iterator<uint64_t>() );
		std::copy( std::istream_iterator<uint64_t>( rxInputStream ), std::istream_iterator<uint64_t>( ), std::back_inserter( sa ) );
#endif // else of ( WORK_WITH_DEPRECTED_SA == 1 )

		/* Check whether we read exactly the required number of elements.
		 */
		if ( sa.size() != uiExpectedSASize )
		{
			throw std::runtime_error( std::string( "Reading suffix array from file system failed due to non matching expected size.") );
		} // if 
#endif
	} // method

public :
	/* Like bwt_occ4, but it computes the counters for two indices simultaneously. 
	 * If k and l belong to the same internal interval (the FM index has a block structure), bwt_2occ4 is more efficient than bwt_occ4.
	 * IMPORTANT: The indices are inclusive, i. e. we count 
	 * IMPORTANT: Requires k <= l
	 */
	void bwt_2occ4( t_bwtIndex k,										// first end index in BWT (k can be negative: (t_bwtIndex)(-1) )
					t_bwtIndex l,										// second end index in BWT (l can be negative: (t_bwtIndex)(-1) )
					bwt64bitCounter cntk[4], bwt64bitCounter cntl[4]	// (outputs) counter for k and l
				  )
	{	/* Adjusted versions of k and l, because $ is not in bwt
		*/
		t_bwtIndex _k = k - (k >= primary);
		t_bwtIndex _l = l - (l >= primary);
		
		if (   _l >> OCC_INTV_SHIFT != _k >> OCC_INTV_SHIFT // l and k belong to different intervals
			|| k == (t_bwtIndex)(-1)						// 
			|| l == (t_bwtIndex)(-1)						// 
		   ) 
		{	/* Interval decomposition step by two calls of bwt_occ4
			 */
			//// std::cout << "2 calls of bwt_occ4 with " << " k:" << k << " l:" << l << "\n";
			bwt_occ4( k, cntk );
			bwt_occ4( l, cntl );
		} // if
		else 
		{	/* The computation of bwt_occ4 but now for two indices in parallel.
			 */
			uint64_t x, y; // originally already 64 bit, but should work with a 32 bit value as well.
			uint32_t *p, tmp, *endk, *endl;

			/* Adjust k and l, because $ is not in bwt
			 */
			k -= (k >= primary);
			l -= (l >= primary);
			
			/* Step 1. (see bwt_occ4)
			 */
			p = bwt_occ_intv(k);
			
			//// std::cout << "k is " << k << "p is " << p << "\n";
			memcpy( cntk, p, 4 * sizeof(bwt64bitCounter) );
			
			//// for ( int i = 0; i < 4; ++i )
			//// 	std::cout << "cntk[" << i << "]=" << cntk[i] << "  ";
			//// std::cout << "\n";
			
			/* Step 2. (see bwt_occ4)
			 */
			p += sizeof(bwt64bitCounter); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
			// prepare cntk[]
			endk = p + ((k >> 4) - ((k&~OCC_INTV_MASK) >> 4));
			endl = p + ((l >> 4) - ((l&~OCC_INTV_MASK) >> 4));

			//// std::cout << "endk - p:" << endk - p << " endl - p:" << endl - p << "\n";
			
			for (x = 0; p < endk; ++p)
			{
				x += __occ_aux4(*p);
			} // for
			
			/* Correct due to precondition k <= l
			 */
			y = x;
			tmp = *p & ~((1U << ((~k & 15) << 1)) - 1);
			x += __occ_aux4(tmp) - (~k & 15);
			
			/* Continue summing up for l
			 */
			for (; p < endl; ++p)
			{
				y += __occ_aux4(*p);
			} // for
			tmp = *p & ~((1U << ((~l & 15) << 1)) - 1);
			y += __occ_aux4(tmp) - (~l & 15);
			
			memcpy(cntl, cntk, 4 * sizeof(bwt64bitCounter)); // could be moved to step 1, without changing semantics
			cntk[0] += x	   & 0xff ; cntl[0] += y       & 0xff;	// A
			cntk[1] += x >> 8  & 0xff ; cntl[1] += y >> 8  & 0xff;	// C
			cntk[2] += x >> 16 & 0xff ; cntl[2] += y >> 16 & 0xff;	// G
			cntk[3] += x >> 24        ; cntl[3] += y >> 24;			// T
		} // else
	} // method

	//markus
	const inline unsigned int getRefSeqLength() const { return uiRefSeqLength; }//function
	//end markus

	/* We keep the reference length private in order to avoid unexpected trouble.
	 * Delivers the length of the reference (pack) that belongs to the current FM-index.
	 */
	uint64_t getRefSeqLength( void )
	{
		return uiRefSeqLength;
	} // method

	/* Delivers the Position in the reference sequence T that belongs to the position k in the BWT.
	 * Uses the suffix array cache ...
	 */
	bwtint_t bwt_sa( bwtint_t uiBWTposition )
	{	/* Check uiBWTposition for out of range
		 */
		assert( ( uiBWTposition >= 0 ) && ( uiBWTposition < static_cast<bwtint_t>(uiRefSeqLength) ) ); // Out of range check for uiBWTposition

		bwtint_t sa = 0;
		const bwtint_t mask = sa_intv - 1; // this is why sa_intv must be a power of 2 
		
		/* We iterate so long until we meet an interval value of the cache.
		 */
		while (uiBWTposition & mask) 
		{
			++sa;
			uiBWTposition = bwt_invPsi( uiBWTposition );
		} // while
		
		/* Without setting bwt->sa[0] = -1 in bwt_cal_sa_step3, the following line should be
		 * changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) 
		 */
		auto uiReturnedPosition = sa + this->sa[uiBWTposition/sa_intv];
		assert( ( uiReturnedPosition >= 0 ) && ( uiReturnedPosition < static_cast<bwtint_t>(uiRefSeqLength) ) ); // Out of range for returned reference position
		
		return uiReturnedPosition;
	} // method

	/* Checks whether the files required for loading a pack does exist on the file system.
	 */
	static bool packExistsOnFileSystem( const std::string &rsPrefix )
	{
		//1 == use std o check for file existance
#if 0 
		return access( std::string(rsPrefix).append(".bwt").c_str(), F_OK ) != -1 &&
			access( std::string(rsPrefix).append(".sa").c_str(), F_OK ) != -1;
#else
		return	   boost::filesystem::exists( std::string(rsPrefix).append(".bwt") )
				&& boost::filesystem::exists( std::string(rsPrefix).append(".sa") );
#endif
	} // method
	

	/* Dump the current FM-Index to two separated files for BWT and SA.
	 */
	void vStoreFM_Index( const boost::filesystem::path &rxFileNamePrefix )
	{
		{	/* Save burrow wheeler transform
			 */
			std::ofstream xOutputStream( rxFileNamePrefix.string() + ".bwt", std::ios::binary | std::ios::trunc | std::ios::out );
			vSaveBWT( xOutputStream );
			xOutputStream.close();
		} // scope
		
		{	/* Save suffix array
			 */
			std::ofstream xOutputStream( rxFileNamePrefix.string() + ".sa", std::ios::binary | std::ios::trunc | std::ios::out );	
			vSaveSuffixArray( xOutputStream );
			xOutputStream.close();
		} // scope
	} // method


//markus

	/* Debug method.
	* Checks the correctness of a single interval.
	*/
	bool debugCheckInterval(SA_IndexInterval &rxInterval,
		const NucleotideSequence &rxQuerySeq,
		const BWACompatiblePackedNucleotideSequencesCollection &rxSequenceCollection,
		bool bOnRevFmIndex = false
		)
	{
		auto uiIntervalSize = rxInterval.x[2];
		if (uiIntervalSize > 200)
		{
			std::cout << "Skip inspection of interval, because it has " << rxInterval.x[2] << " many element and sequence length: " << rxInterval.sizeOfMatch() << std::endl;
			return true;
		} // if

		auto uiMatchSize = rxInterval.sizeOfMatch();
		auto uiQuerySeqPosition = rxInterval.startIndexOfMatchInQuery();
		std::cout << "before getSubsliceFromTo" << std::endl;

		auto xQuerySeqSlice = rxQuerySeq.fromTo(uiQuerySeqPosition, uiQuerySeqPosition + uiMatchSize);

		std::cout << "INSPECT INTERVAL with " << rxInterval.x[2] << " many element and sequence length: " << rxInterval.sizeOfMatch() << std::endl;
		/* Check all matches within the current interval.
		*/
		for (t_bwtIndex uiCounter = 0; uiCounter < uiIntervalSize; ++uiCounter)
		{
			std::cout << "SA index,  rxInterval.x[0]: " << rxInterval.x[0] << " counter: " << uiCounter << std::endl;
			auto uiRefSeqPosition = bwt_sa(rxInterval.x[0] + uiCounter);
			
			//by markus: using a fm index of the reversed reference in order to do forward extension
			if (bOnRevFmIndex)
				uiRefSeqPosition = rxSequenceCollection.uiUnpackedSizeForwardPlusReverse() - ( uiRefSeqPosition + uiMatchSize );

			std::cout << "Extract sequence" << std::endl;
			NucleotideSequence xSequence; // Will keep the required unpacked section of the reference sequence.
			rxSequenceCollection.vExtractSubsection(uiRefSeqPosition, uiRefSeqPosition + uiMatchSize, xSequence);
			auto xRefSeqSlice = xSequence.fullSequenceAsSlice();

			if (!(xRefSeqSlice.equals(xQuerySeqSlice)))
			{
				/* Something is wrong....
				* Dump some info about the problem.
				*/
				std::cout << "Data for wrong interval:" << std::endl;
				std::cout << "Size of interval: " << rxInterval.x[2] << std::endl;
				std::cout << uiRefSeqPosition << " " << uiMatchSize << " " << rxInterval.startIndexOfMatchInQuery() << std::endl;

				std::cout << "refSeq: " << *xRefSeqSlice.asSequenceOfACGT() << std::endl;
				std::cout << "querySeq: " << *xQuerySeqSlice.asSequenceOfACGT() << std::endl;

				//exit(0);

				return false;
			} // if
			else
			{
				std::cout << "confirmed match at " << uiRefSeqPosition << std::endl;
			}//else
		} // for (all matches)

		return true;
	}//debug function

//end markus

	/* Load an FM-Index previously stored by vStoreFM_Index.
	 */
	void vLoadFM_Index( const boost::filesystem::path &rxFileNamePrefix )
	{
		{	if ( !boost::filesystem::exists( rxFileNamePrefix.string() + ".bwt" ) )
			{
				throw std::runtime_error( "File opening error: " + rxFileNamePrefix.string() + ".bwt" );
			} // if

			std::ifstream xInputFiletream( rxFileNamePrefix.string() + ".bwt", std::ios::binary | std::ios::in );
			if ( xInputFiletream.fail() ) // check whether we could successfully open the stream
			{	
				throw std::runtime_error( "Opening of BWT of FM index failed: " + rxFileNamePrefix.string() + ".bwt" );
			} // if

			vRestoreBWT( xInputFiletream );
			xInputFiletream.close();
		} // scope
		
		/* The order is important over here, because during loading the suffix array we check with compatibility with the BWT.
		 */
		{	std::ifstream xInputFileStream( rxFileNamePrefix.string() + ".sa", std::ios::binary | std::ios::in );
			if ( xInputFileStream.fail() ) // check whether we could successfully open the stream
			{
				throw std::runtime_error( "Opening of suffix array of FM index failed: " + rxFileNamePrefix.string() + ".sa" );
			} // if
			vRestoreSuffixArray( xInputFileStream );
			xInputFileStream.close();
		} // scope
	} // method


	/* Debug function for comparing BWT.
	 * BWT can be build by several different functions. Here we can check for correctness.
	 * So long it does not compare the BWT sequence itself.
	 */
	bool operator== ( const FM_Index &rxOtherFM_Index )
	{
		std::string sErrorText = ( L2 != rxOtherFM_Index.L2 ) 
								 ?  "Different L2"
								 :	( primary != rxOtherFM_Index.primary ) 
									?  "Different primary"
									:	( uiRefSeqLength != uiRefSeqLength )
										? "Different seq_len"
#if ( confWORK_WITH_DEPRECTED_BWT == 1 )	
										:	( bwt_size != bwt_size )
#else
										:	( bwt != rxOtherFM_Index.bwt )
#endif
											?	"Different bwt_size"
											: ( sa != rxOtherFM_Index.sa ) 
												? "Different suffix arrays"
												: "";
		if ( sErrorText != "" )
		{
			BOOST_LOG_TRIVIAL( info ) << "BWT different: " << sErrorText;
		} // if
		
		return sErrorText == "";
	} // method

	/* Default constructor. (Initializes the fix count-table)
	 */
	FM_Index()
		: L2 ({ { 0, 0, 0, 0, 0, 0 } }),
#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
		  bwt( nullptr ),
#else
		  bwt(), // initialize the vector keeping the BWT
#endif
#if ( confWORK_WITH_DEPRECTED_SA == 1 )
		  sa( nullptr )
#else
		  sa() // initialize the vector keeping the suffix array
#endif
	{
		_vInitCountTable();
	} // default constructor

	/* FM-Index constructor. Builds a FM index on foundation of rxSequence. 
	 */
	FM_Index( const NucleotideSequence &rxSequence ) 
		: FM_Index() // call the default constructor
	{
		build_FM_Index( rxSequence );
	} // constructor

	/* FM-Index constructor. Builds a FM index on foundation of a given sequence collection. 
	 */
	FM_Index( const BWACompatiblePackedNucleotideSequencesCollection &rxSequenceCollection, // the pack for which we require a BWT
			  unsigned int uiAlgorithmSelection = 2 // 2 -> automatic algorithm selection
			)
		: FM_Index() // call the default constructor
	{
		build_FM_Index( rxSequenceCollection, uiAlgorithmSelection );
	} // constructor

	/* Destructor; responsilbe for deallocating memory.
	 */
	~FM_Index()
	{
#if ( confWORK_WITH_DEPRECTED_BWT == 1 )
		if ( bwt != nullptr )
		{	
 #if ( defBWT_USE_NEW_OPERATOR == 1 )
			delete[] bwt;
 #else
			free( bwt ); // FIX ME: Use unique pointer for bwt.
 #endif
		} // if
#endif
#if ( confWORK_WITH_DEPRECTED_SA == 1 )
		if ( sa != nullptr )
		{
	#if ( defBWT_USE_NEW_OPERATOR == 1 )
				delete[] sa;
	#else
				free( sa ); // FIX ME: Use unique pointer for sa.
	#endif
		} // if
#endif
	} // destructor
}; // class FM_Index


class FM_IndexContainer: public Container
{
public:
	std::shared_ptr<FM_Index> pIndex;

	FM_IndexContainer()
			:
		pIndex(new FM_Index())
	{}//constructor

	FM_IndexContainer(const FM_IndexContainer *pCpyFrom)
			:
		pIndex(pCpyFrom->pIndex)
	{}//copy constructor

	/*used to identify the FM_indexWrapper datatype in the aligner pipeline*/
    ContainerType getType(){return ContainerType::fM_index;}
	

	void vLoadFM_Index( std::string sPrefix  )
	{
		std::string sPath(sPrefix);
		pIndex->vLoadFM_Index(sPath);
	}//function

	static bool packExistsOnFileSystem(const std::string sPrefix )
	{
		return FM_Index::packExistsOnFileSystem(sPrefix);
	} // method
	

	/* wrap the function in oder to make it acessible to pyhton */
	void vStoreFM_Index(  std::string sPrefix  )
	{
		std::string sPath(sPrefix);
		pIndex->vStoreFM_Index(sPath);
	}//function

	std::shared_ptr<Container> copy()
    {
		return std::shared_ptr<Container>(new FM_IndexContainer(this));
	}//function
};//class

//function called in order to export this module
void exportFM_index();