/**
 * @file sw_sse_avx.h
 * @brief Implements the smith waterman algorithm.
 * @author Arne Kutzner
 */

#pragma once
/* IMPORTANT NOTE: 
 * For getting AVX2 code with g++ the compiler switch "-mavx2" has to be set.
 * AVX code (even if limited to 128bit) don't work on system having SSE support only.
 * WEB-Info regarding Intrinsics https://software.intel.com/sites/landingpage/IntrinsicsGuide/
 * Observations: With 128 bit: AVX is slightly slower than SSE
 */

#if _MSC_VER
	#define NOMINMAX
	#include <windows.h>
#endif

#include <vector>
#include <limits> 
// #include "sw_common.h"
#include "util/exception.h" // code throws aligner exceptions
#include "container/nucSeq.h" // sequence slices

/* Current Problem: Horizontal and vertical gap must be identical! */
#define GAP_EXT_HORIZONAL ( 1 )

/* Visual C++ and g++ indicate the setting of the AVX2 flag via the __AVX2__ symbol.
 */
#if (__AVX2__ == 1)
	#pragma message("Compilation of SIMD aligner using AVX2 intrinsics.")
	#define __USE_AVX2__ ( 1 )
#else
	#pragma message("Compilation of SIMD aligner using SSE intrinsics.")
	#define __USE_AVX2__ ( 0 )
#endif

#if (__USE_AVX2__ == 1)
	// AVX - 256 bit registers
	#ifdef __GNUC__
		#include <x86intrin.h> // AVX header
	#endif

	/* Shifting for AVX2 256 bit vectors:
	* See: https://stackoverflow.com/questions/20775005/8-bit-shift-operation-in-avx2-with-shifting-in-zeros
	*/
	template <unsigned int N> inline __m256i _mm256_shift_left(__m256i a)
	{
		__m256i mask = _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 3, 0));
		return _mm256_alignr_epi8(a, mask, 16 - N);
	}

	#define T__mXXXi __m256i
	#define _mmXXX_set1_epi32(x) _mm256_set1_epi32((x))
	#define _mmXXX_set1_epi16(x) _mm256_set1_epi16((x))
	#define _mmXXX_adds_epi16(x,y) _mm256_adds_epi16((x), (y))
	#define _mmXXX_cmpgt_epi16(x,y) _mm256_cmpgt_epi16((x), (y))
	#define _mmXXX_subs_epu16(x,y) _mm256_subs_epu16((x), (y))
	#define _mmXXX_max_epi16(x,y) _mm256_max_epi16((x), (y))
	#define _mmXXX_movemask_epi8(x) _mm256_movemask_epi8((x))
	#define _mm_store_siXXX(x,y) _mm256_store_si256((x), (y))
	#define _mm_load_siXXX(x) _mm256_load_si256((x))
	#define _mm_slli_siXXX(x,y) _mm256_shift_left<y>((x))

	/* Tricky extraction of the maximum among 16 parallel short-integers using 9 SIMD-statements.
	*/  
	#define __max_8(ret, xx) do { \
		(xx) = _mm256_max_epi16((xx), _mm256_permute2f128_si256((xx), (xx), 1)); \
		(xx) = _mm256_max_epi16((xx), _mm256_srli_si256((xx), 8)); \
		(xx) = _mm256_max_epi16((xx), _mm256_srli_si256((xx), 4)); \
		(xx) = _mm256_max_epi16((xx), _mm256_srli_si256((xx), 2)); \
		(ret) = _mm256_extract_epi16((xx), 0); \
	} while (0)
#else
	// SSE2 - 128 bit registers
	#ifdef __GNUC__
		#include <emmintrin.h> // SSE2 header
	#endif

	#define T__mXXXi __m128i
	#define _mmXXX_set1_epi32(x) _mm_set1_epi32((x))
	#define _mmXXX_set1_epi16(x) _mm_set1_epi16((x))
	#define _mmXXX_adds_epi16(x,y) _mm_adds_epi16((x), (y))
	#define _mmXXX_cmpgt_epi16(x,y) _mm_cmpgt_epi16((x), (y))
	#define _mmXXX_subs_epu16(x,y) _mm_subs_epu16((x), (y))
	#define _mmXXX_max_epi16(x,y) _mm_max_epi16((x), (y))
	#define _mmXXX_movemask_epi8(x) _mm_movemask_epi8((x))
	#define _mm_store_siXXX(x,y) _mm_store_si128((x), (y))
	#define _mm_load_siXXX(x) _mm_load_si128((x))
	#define _mm_slli_siXXX(x,y) _mm_slli_si128((x), (y))

	/* Tricky extraction of the maximum among 8 parallel short-integers using 7 SIMD-statements.
	*/
	#define __max_8(ret, xx) do { \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
		(ret) = _mm_extract_epi16((xx), 0); \
	} while (0)
#endif

// SIMD vector size as decided at comile time (SSE2 : 128 bit, AVX : 256 bit)
#define __mXXXi_SIZE sizeof(T__mXXXi)

/* WARNING: This parallel aligner works only for type int16_t (16 bit signed integer)
 * SSE2 recognition: https://stackoverflow.com/questions/6121792/how-to-check-if-a-cpu-supports-the-sse3-instruction-set
 */
template<class SCORE_TP,	// must be int16_t for 
		 class T_size_t		// unsigned interger type for iterators etc..
		>
struct SW_SIMD_Aligner 
{
	const T_size_t uiQueryLen;
	const T_size_t uiNumberOfBlocks;

	// const SimilarityMatrix<SCORE_TP, T_size_t> xSimilarityMatrix;

	const SmithWatermanParamaterSet<SCORE_TP> &pSWparameterSetRef;
	
	/* For __m128i and int16_t we have 8 values per block.
	 * For __m256i and int16_t we have 16 values per block.
	 * For __m128i and int32_t we have 4 values per block.
	 */
	static const int scoresPerBlock = sizeof(T__mXXXi) / sizeof(SCORE_TP);

	// Alloated memory must not be aligned!
	T__mXXXi *pStartOfReservedMemory, *pAlignedAllocatedMemoryRef;
	T__mXXXi *pH_VectorPreviousRow, *pH_VectorCurrentRow, *pE, *pH_VectorWithGlobalMaximum;

	/* Constructor
	 */
	SW_SIMD_Aligner( const NucSeq &rQuerySequence, // query sequence of alignment
					 const SmithWatermanParamaterSet<SCORE_TP> &SWparameterSet ) // alignment parameter
		: uiQueryLen( rQuerySequence.getSize() ),
		  uiNumberOfBlocks( ( rQuerySequence.getSize() + scoresPerBlock - 1 ) / scoresPerBlock ),
		  pSWparameterSetRef( SWparameterSet )
	{	
		if( sizeof(SCORE_TP) != 2 )
		{	// The current design works for int16_t only.
			throw AlignerException("Datatype for scores has wrong dimension (Must be int16_t)");
		} // if
		
		if ( SWparameterSet.iPossibleMaxScore(uiQueryLen) > std::numeric_limits<SCORE_TP>::max() )
		{	// The query might create a maximum beyond the max of SCORE_TP.
			throw AlignerException("Possible maximum for query exceeds maximum of scoring datatype");
		} // if

		/* Allocation memory for the query profile and the 4 auxiliary vectors.
		 * (We allocate this memory as a single block!)
		 */
		pStartOfReservedMemory = (T__mXXXi*)malloc(  (__mXXXi_SIZE - 1)	// space for alignment 
											       + sizeof(T__mXXXi) * uiNumberOfBlocks * ( SWparameterSet.uiAlphabetSize + 4 ) // +4 for H0, H1, E, MAX
												  ); 
		if( pStartOfReservedMemory == NULL )
		{
			throw AlignerException( "Memory allocation for SIMD alignment failed" );
		} // if
		
		/* Align with respect to the reserved memory for efficiency improvements.
		 * For T__mxxxi we need 16 byte aligned memory or things are inefficent.
		 */
		pAlignedAllocatedMemoryRef = (T__mXXXi*)(((T_size_t)pStartOfReservedMemory + (__mXXXi_SIZE - 1)) / __mXXXi_SIZE * __mXXXi_SIZE);

		/* Memory layout:
		 * 1. scoring profile for query - size : numberOfBlocks * SWparameterSet.uiAlphabetSize
		 * 2. H vector - size : numberOfBlocks
		 * 3. E vector - size : numberOfBlocks
		 * 4. H vector max - size : numberOfBlocks
		 */
		pH_VectorPreviousRow = pAlignedAllocatedMemoryRef + (uiNumberOfBlocks * SWparameterSet.uiAlphabetSize);
		pH_VectorCurrentRow = pH_VectorPreviousRow + uiNumberOfBlocks;
		pE = pH_VectorCurrentRow + uiNumberOfBlocks;
		pH_VectorWithGlobalMaximum = pE + uiNumberOfBlocks;

		initQueryProfile(rQuerySequence, SWparameterSet);
	} // constructor

	/* The destructor deallocates all reserved memory.
	 */
	~SW_SIMD_Aligner()
	{	// Aligned mem not used an longer ...
		if( pStartOfReservedMemory != NULL )
		{
			free( pStartOfReservedMemory );
		} // if
	} // destructor

	/* Write the query profile to memory.
	 * Similar technique to the serial approach.
	 */
	void initQueryProfile(const NucSeq &rQuerySequence,
						  const SmithWatermanParamaterSet<SCORE_TP> &SWparameterSet)
	{
		/* Initialize the Scoring Profile for the query as described in the paper.
		* The Scoring Profile starts at the beginning of the allocated memory.
		*/
		SimilarityMatrix<SCORE_TP> xSimilarityMatrix(SWparameterSet.iWeightMatch, SWparameterSet.iWeightMismatch, SWparameterSet.uiAlphabetSize);
		SCORE_TP *profileReference = (SCORE_TP*)pAlignedAllocatedMemoryRef;
		
		for (unsigned int alphabetIterator = 0; alphabetIterator < SWparameterSet.uiAlphabetSize; ++alphabetIterator) 
		{
			const SCORE_TP *scroringTableRow = xSimilarityMatrix.pSimilarityMatrixRef + (alphabetIterator * SWparameterSet.uiAlphabetSize);
			for( T_size_t i = 0; i < uiNumberOfBlocks; ++i ) 
			{
				for( T_size_t k = i; k < uiNumberOfBlocks * scoresPerBlock; k += uiNumberOfBlocks ) // p iterations
				{	/* We take std::numeric_limits<SCORE_TP>::min() for the "non existing columns".
					 * (Don't take 0 instead of min, because then we get faulty values.)
					 */
					*profileReference++ = (k >= uiQueryLen ? std::numeric_limits<SCORE_TP>::min() // non existing column
						: scroringTableRow[(rQuerySequence.getSequenceRefInHost())[k]]
						);
				} // for k - iteration over p, where p = 16 or 8 or 4
			} // for i - iteration over segment size
		} // for alphabetIterator
	} // method 
	
	/* The align method for SW-alignments - 16 bit version of the parallel SW implementation.
	 * The first gap costs -( _gapo + _gape) !
	 * The first row and column have both index 0 (important if a maximum is monitored in row 0)
	 * DatabaseSequence should be delivered in the already translated form.
	 * You have to deliver a fresh maxScoreProfiler
	 */
	
	/* Align for 16 bit sized scores.
	 * Returns maximum score.
	 */
	SCORE_TP align( const NucSeq &pReferenceSeq,
					 std::vector<T_size_t> &rvMaxScorePositions // vector will keep the positions of occurences of max score
#if (DO_CHECKS == 1)
					 , std::vector<T_scoring> &rvSwRowMaxima // vector collecting row maxima (for debugging)
#endif
					 ) 
	{
		auto pReference = pReferenceSeq.getSequenceRefInHost(); // uint_8t
		auto uiLenReference = pReferenceSeq.getSize();
		
		SCORE_TP iOverallMaxScore = 0;

		// Repeatedly used zero vector
		T__mXXXi zero = _mmXXX_set1_epi32(0);  
	
		/* Initialize the vectors for the gaps.
		 * SCORE_TP must be equal to int16_t or we are in trouble here
		 */
		T__mXXXi gapoePar8 = _mmXXX_set1_epi16( (int16_t)(pSWparameterSetRef.iGapOpen + pSWparameterSetRef.iGapExtend) ); 
		T__mXXXi gapePar8 = _mmXXX_set1_epi16( (int16_t)pSWparameterSetRef.iGapExtend ); 

		T__mXXXi gapoePar8F = _mmXXX_set1_epi16((int16_t)(pSWparameterSetRef.iGapOpen + GAP_EXT_HORIZONAL)); //// AK
		T__mXXXi gapePar8F = _mmXXX_set1_epi16((int16_t)(GAP_EXT_HORIZONAL)); //// AK
		
		/* We initialize the E, H0 and Hmax vectors with zero.
		 * Idea: here we could take something build in function.
		 */
		for( T_size_t uxIterator = 0; uxIterator < uiNumberOfBlocks; ++uxIterator ) 
		{
			_mm_store_siXXX(pE + uxIterator, zero); 
			
			/* In the first row we need 0 form the "previous row" !
			 */
			_mm_store_siXXX(pH_VectorPreviousRow + uxIterator, zero); 
			
			_mm_store_siXXX(pH_VectorWithGlobalMaximum + uxIterator, zero); 
		}

		/* The outer loop iterates over the matrix rows. (The database sequence)
		 */
		for( T_size_t uiRow = 0; uiRow < uiLenReference; ++uiRow ) 
		{
			T__mXXXi f = zero; 
			T__mXXXi maxima = zero; 

			/* Pick the profile row that belongs to the current symbol.
			 */
			T__mXXXi *profileReferenceForRow = pAlignedAllocatedMemoryRef + (pReference[uiRow] * uiNumberOfBlocks );

			T__mXXXi h = _mm_load_siXXX(pH_VectorPreviousRow + uiNumberOfBlocks - 1); 

			/* Shift the 128-bit value in a left by 2 bytes while shifting in zeros
			 * WARNING: 16 bit mechanics.
			 */
			h = _mm_slli_siXXX(h, 2); 
		
			/* The inner loop iterates over the blocks 
			 * Because of the segmentation it considers 8 values in parallel.
			 */
			for( T_size_t j = 0; j < uiNumberOfBlocks; ++j )
			{	/* SW cells are computed in the following order:
				 *   F(i,j+1) = max{H(i,j)-q, F(i,j) - r}
				 */
				/* Compute H(i,j) = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)} (parallel for 16/8/4 elements)
				 * Note that at the beginning, h = H(i - 1,j - 1)
				 */
				h = _mmXXX_adds_epi16( h, *(profileReferenceForRow++) ); 
				
				T__mXXXi e = _mm_load_siXXX( pE + j ); 
				h = _mmXXX_max_epi16( h, e ); 
				h = _mmXXX_max_epi16( h, f ); 
				
				/* At this point we have the correct h-values in h, regarding finding the maximum.
				 * However, f is not integrated here.
				 */
				maxima = _mmXXX_max_epi16( maxima, h ); 
				_mm_store_siXXX(pH_VectorCurrentRow + j, h); 
				
				/* We compute and store the <e>-values for the next row.
				 * E(i+1, j) = max{H(i, j) - g_init, E(i, j) - g_ext}
				 */
				h = _mmXXX_subs_epu16(h, gapoePar8); 
				e = _mmXXX_subs_epu16(e, gapePar8); 
				e = _mmXXX_max_epi16(e, h); 
				_mm_store_siXXX(pE + j, e); 
				
				/* We compute the <f>-values (left cells). These values can be wrong, we correct this in the lazy-F loop. 
				 * F(i, j+1) = max{H(i, j) - g_init, F(i, j) - g_ext}
				 */
				f = _mmXXX_subs_epu16(f, gapePar8F );  //// AK
				f = _mmXXX_max_epi16(f, h); 
				
				/* We load h for next iteration. Important is that we do this here and not at the beginning of the loop.
				 */
				h = _mm_load_siXXX(pH_VectorPreviousRow + j); 
			} // for j

			/* The lazy F-loop as described in the paper. This loop is responsible for fixing wrong f values.
			 * This loop can become quite expensive in specific situations.
			 * Observation: short queries (20 nt) create large k values. 
			 */
			for( int k = 0; k < __mXXXi_SIZE / 2; ++k )
			{
				f = _mm_slli_siXXX(f, 2); 
				for( T_size_t j = 0; j < uiNumberOfBlocks; ++j ) 
				{
					h = _mm_load_siXXX( pH_VectorCurrentRow + j ); 
					h = _mmXXX_max_epi16( h, f ); 
					_mm_store_siXXX( pH_VectorCurrentRow + j, h ); 
					
					h = _mmXXX_subs_epu16( h, gapoePar8F ); //// AK
					f = _mmXXX_subs_epu16( f, gapePar8F ); //// AK
					
					/* _mm_movemask_epi8 : Create mask from the most significant bit of each 8-bit element in a, and store the result in dst.
					 * Hint: movemask is simply a trick in order to finally get an integer that can be used for comparison.
					 * _mm_cmpgt_epi16 : Compare packed 16-bit integers in a and b for equality, and store the results in dst
					 */
					if( !_mmXXX_movemask_epi8( _mmXXX_cmpgt_epi16( f, h ) ) )  
					{
						goto exit_loop;
					}
				} // for j
			} // for k

			exit_loop:

			/* Get the max for row from the SIMD vector
			 */
			SCORE_TP iMaxScoreOfCurrentRow;
			__max_8( iMaxScoreOfCurrentRow, maxima );
			//// std::cout << iMaxScoreOfCurrentRow << std::endl;

#if (DO_CHECKS == 1)
			// Record computed row maximum 
			rvSwRowMaxima.push_back(iMaxScoreOfCurrentRow);
#endif
			
			/* Tracing of global maxima
			 */
			if( iMaxScoreOfCurrentRow >= iOverallMaxScore )
			{
				if( iMaxScoreOfCurrentRow > iOverallMaxScore )
				{	// fresh overall maximum detected 
					rvMaxScorePositions.clear();
					iOverallMaxScore = iMaxScoreOfCurrentRow;
				}
				rvMaxScorePositions.push_back(uiRow);
			} // if
			
			/* Swap the current and the previous H-vector.
			 */
			profileReferenceForRow = pH_VectorCurrentRow; 
			pH_VectorCurrentRow = pH_VectorPreviousRow; 
			pH_VectorPreviousRow = profileReferenceForRow;
		} // for uxIterator (row iteration)

		return iOverallMaxScore;
	} // end alignment method
}; // end struct