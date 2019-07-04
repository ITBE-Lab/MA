/**
 * @file fMIndex.h
 * @brief Implements the FM_index class.
 * @author Arne Kutzner
 */
#ifndef FM_INDEX_H
#define FM_INDEX_H

#include "container/bwt_large.h"
#include "container/is.h"
#include "container/pack.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <chrono> // time required for temporary filename construction
#include <cstdint>
#include <mutex>
#include <set>
#include <utility>
/// @endcond

namespace libMA
{
//@todo should have a version check built in

// extern Configuration xGlobalConfiguration;

typedef uint64_t bwt64bitCounter;
typedef int64_t t_bwtIndex; // IMPORTANT: We can have -1 in the context of occurrence counting.

/* requirement: (OCC_INTERVAL%16 == 0); please DO NOT change this line because some part of the code
 * assume OCC_INTERVAL=0x80 TO DO: Put these as constants into class FM_index.
 */
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL ( 1LL << OCC_INTV_SHIFT )
#define OCC_INTV_MASK ( OCC_INTERVAL - 1 )

/**
 * @brief A suffix array interval.
 * @details
 * The BWT interval class represents intervals with respect the suffix array!
 * See http://en.wikipedia.org/wiki/FM-index (intervals within the F-L matrix )
 * @ingroup container
 */
class SAInterval : public Container, public Interval<t_bwtIndex>
{
  private:
    t_bwtIndex startOfComplement;

  public:
    SAInterval( t_bwtIndex start, t_bwtIndex startOfComplement, t_bwtIndex size )
        : Interval( start, size ), startOfComplement( startOfComplement )
    {} // constructor

    SAInterval( ) : Interval( ), startOfComplement( 0 )
    {} // constructor

    SAInterval( const SAInterval& other ) : Interval( other ), startOfComplement( other.startOfComplement )
    {} // copy constructor

    // overload
    bool canCast( std::shared_ptr<Container> c ) const
    {
        return std::dynamic_pointer_cast<SAInterval>( c ) != nullptr;
    } // function

    // overload
    std::string getTypeName( ) const
    {
        return "SAInterval";
    } // function

    // overload
    std::shared_ptr<Container> getType( ) const
    {
        return std::shared_ptr<Container>( new SAInterval( ) );
    } // function

    /**
     * @brief Switch to the respective reverse complement.
     * @details
     * A SAInterval is the same size as the interval for the reverse complement sequence.
     * This fact can be used to implement forward extension.
     * @Note This only works when combined with FMD-indices.
     */
    SAInterval revComp( ) const
    {
        return SAInterval( startOfComplement, start( ), size( ) );
    } // function

    /**
     * @brief returns the start position of the SAInterval created using the reverse
     * complement of the sequence used for this interval.
     */
    t_bwtIndex startRevComp( ) const
    {
        return startOfComplement;
    } // function

    /*
     * @brief copys from another SAInterval.
     */
    inline SAInterval& operator=( const SAInterval& rxOther )
    {
        Interval::operator=( rxOther );
        assert( size( ) == rxOther.size( ) );
        startOfComplement = rxOther.startOfComplement;
        return *this;
    } // operator

    /*
     * @brief compares two Intervals.
     * @returns true if start and size are equal, false otherwise.
     */
    inline bool operator==( const SAInterval& rxOther )
    {
        return Interval::operator==( rxOther ) && startOfComplement == rxOther.startOfComplement;
    } // operator
}; // class ( SAInterval )

typedef unsigned char ubyte_t;


class SuffixArrayInterface : public Container
{
  public:
    virtual SAInterval EXPORTED extend_backward(
        // current interval
        const SAInterval& ik,
        // the character to extend with
        const uint8_t c );

    virtual SAInterval init_interval(
        // the character to init with
        const uint8_t c )
    {
        throw std::runtime_error( "unimplemented!" );
    } // method

    virtual uint64_t getRefSeqLength( void ) const
    {
        throw std::runtime_error( "unimplemented!" );
    } // method
}; // class

class SuffixArray : public SuffixArrayInterface
{}; // class

/**
 * @brief A Suffix array.
 * @details
 * FM-indices are central with the alignment process.
 * The original data-structure was part of the BWA-code.
 * @ingroup container
 */
class FMIndex : public SuffixArrayInterface
{
  public:
    /* C(), cumulative counts (part of the FM index)
     * Initialization happens in the constructor.
     */
    std::array<uint64_t, 6> L2; // uint64_t L2[6]; // 0, A, C, G, T, N

    /* S^{-1}(0), or the primary index of BWT.
     * BWT-index, position of $
     */
    t_bwtIndex primary;

    SAInterval EXPORTED getInterval( std::shared_ptr<NucSeq> pQuerySeq );

    t_bwtIndex EXPORTED get_ambiguity( std::shared_ptr<NucSeq> pQuerySeq );

    bool EXPORTED testSaInterval( std::shared_ptr<NucSeq> pQuerySeq, const Pack& rPack );

    bool EXPORTED test( const Pack& rPack, unsigned int uiNumTest );

  protected:
    typedef int64_t bwtint_t;

    /* Count-table for 4 values of 2 bit.
     * Table is initialized by _vInitCountTable
     */
    uint32_t _aCountTable[ 256 ];

    std::vector<uint32_t> bwt;

    /* Size of the reference sequence that was used for BWT construction.
     * Length of the reference sequence (forward plus reverse strand)
     */
    uint64_t uiRefSeqLength;

    // interval size (must be power of 2 and is normally 32, because of 64 bit / 2 == 32 bit)
    int sa_intv;

    /* suffix array
     */
    std::vector<bwtint_t> sa;

    /** Core BWT construction for short sequences.
     * Initializes the BWT using the given nucleotide sequence.
     * We have a serious problem with ambiguous bases over here, because the (compressed) BWT can't
     * represent them. WARNING: Do not pass sequences comprising ambiguous symbols (e.g. symbol
     * 'N').
     */
    void EXPORTED bwt_pac2bwt_step1( const NucSeq& fn_pac_arg );

/** Retrieve character at position k from the $-removed packed BWT without occurrence counter.
 * Used in the context of step 2 (bwt_bwtupdate_core_step2) for retrieving characters.
 * REMARK: Only correct in the case of OCC_INTERVAL==0x80
 */
#define bwt_B00( k ) ( bwt[ ( k ) >> 4 ] >> ( ( ~(k)&0xf ) << 1 ) & 3 )

    /** Injection of counting blocks into the BWT.
     * Inserts the required space for 4 uint_64 counters each 8 uint_32 blocks and initializes the
     * counters according to the counting. FIX ME: If we get exception over here, we have a memory
     * leak -> use unique poiters are vectors
     */
    void EXPORTED bwt_bwtupdate_core_step2( );

/** Retrieve character at position k from the $-removed packed BWT, which comprises occurrence
 * counter. (Note that bwt_t::bwt is not exactly the BWT string and therefore this macro is called
 * bwt_B0 instead of bwt_B.) REMARK: Only correct in the case of OCC_INTERVAL==0x80
 */
#define bwt_bwt( k ) ( bwt[ ( ( k ) >> 7 << 4 ) + sizeof( bwtint_t ) + ( ( (k)&0x7f ) >> 4 ) ] )
#define bwt_B0( k ) ( bwt_bwt( k ) >> ( ( ~(k)&0xf ) << 1 ) & 3 )

    /** The function Occ(c, k) is the number of occurrences of character c in the prefix L[1..k].
     * Ferragina and Manzini showed[1] that it is possible to compute Occ(c, k) in constant time.
     * ( Paolo Ferragina and Giovanni Manzini (2000). "Opportunistic Data Structures with
     * Applications". Proceedings of the 41st Annual Symposium on Foundations of Computer Science.
     * p.390. ) y contains 32 nucleotides encoded in 2 bits each. c is a single character A/C/G/T =
     * 0/1/2/3 Returned value is the number of occurrences of this character in y.
     */
    inline int __occ_aux( uint64_t y, int c )
    {
        // reduce nucleotide counting to bits counting
        y = ( ( c & 2 ) ? y : ~y ) >> 1 & ( ( c & 1 ) ? y : ~y ) & 0x5555555555555555ull;

        // count the number of 1s in y
        y = ( y & 0x3333333333333333ull ) + ( y >> 2 & 0x3333333333333333ull );

        return ( ( y + ( y >> 4 ) ) & 0xf0f0f0f0f0f0f0full ) * 0x101010101010101ull >> 56;
    } // method

    /** Counts the number of occurrences of the character c in the BWT until position k.
     * bwt_occ counts for a single character, while bwt_occ counts for all 4 characters (A/C/G/T)
     * simultaneously. k should be within [0..seqlen], the value of k shall be not equal to primary
     * (if you want to get the value at primary position deliver (bwtint_t)(-1)
     */
    inline bwtint_t bwt_occ( t_bwtIndex k, ubyte_t c )
    {
        bwtint_t n;
        uint32_t *p, *end;

        if( k == (t_bwtIndex)uiRefSeqLength )
        {
            return L2[ c + 1 ] - L2[ c ];
        } // if

        /* assertion: k != primary
         */
        if( k == ( bwtint_t )( -1 ) )
        {
            return 0;
        } // if

        k -= ( k >= primary ); // because $ is not in bwt

        /* retrieve Occ at k/OCC_INTERVAL
         */
        n = ( (bwtint_t*)( p = bwt_occ_intv( k ) ) )[ c ];
        p += sizeof( bwtint_t ); // jump to the start of the first BWT cell

        /* calculate Occ up to the last k/32
         */
        end = p + ( ( ( k >> 5 ) - ( ( k & ~OCC_INTV_MASK ) >> 5 ) ) << 1 );
        for( ; p < end; p += 2 )
        {
            n += __occ_aux( (uint64_t)p[ 0 ] << 32 | p[ 1 ], c );
        } // for

        // calculate Occ
        n += __occ_aux( ( (uint64_t)p[ 0 ] << 32 | p[ 1 ] ) & ~( ( 1ull << ( ( ~k & 31 ) << 1 ) ) - 1 ), c );
        if( c == 0 )
        {
            n -= ~k & 31; // corrected for the masked bits
        } // if

        return n;
    } // method

    /** Computes the inverse CSA. (Back to front projection in BWT)
     * @param k Position in BWT
     * @return Position of the character at k in F
     */
    inline bwtint_t bwt_invPsi( bwtint_t k )
    { /* Adapt, whether we are before or after the primary.
       */
        bwtint_t x = k - ( k > primary );

        /* Character in BWT at primary adjusted position
         */
        uint8_t c = bwt_B0( x );

        /* bwt_occ will adjust k one again!
         */
        x = L2[ c ] + bwt_occ( k, c );

        return ( k == primary ) ? 0 : x;
    } // method

    /** Small helper function used in the context of some assertion.
     */
    bool ispowerof2( unsigned int x )
    {
        return x && !( x & ( x - 1 ) );
    } // method

    /** Creation of the suffix array (SA is organized on the foundation of 64 bit counters)
     * @param intv Interval size for suffix array (must be a power of 2)
     */
    void EXPORTED bwt_cal_sa_step3( unsigned int intv );

    /** BWT construction for packs. Is is possible to choose among two different algorithms.
     * uiAlgorithmSelection : 0 manually selected algorithm for small inputs
     *                          1 manually selected algorithm for large inputs
     *                          2 automatic selection on foundation of input size
     */
    void EXPORTED build_FMIndex( const Pack& rxSequenceCollection, // the pack for which we compute a BWT
                                 unsigned int uiAlgorithmSelection = 2 // 2 -> automatic algorithm selection
    );

    /** Builds up a FMIndex for the given input sequence.
     * REMARK: The sequence should be free of ambiguous bases (character N).
     */
    void build_FMIndex( const NucSeq& fn_pac )
    {
        /* Construction of core BWT.
         */
        bwt_pac2bwt_step1( fn_pac );

        /* Step 2 and step 3 of FM-index creation
         */
        vPostProcessBWTAndCreateSA( );
    } // method

    /** Step 2 and step 3 of FM-index creation.
     */
    void vPostProcessBWTAndCreateSA( void )
    {
        /* Insertion of occurrence data into outcome of step 1.
         */
        bwt_bwtupdate_core_step2( );

        /* Computation of the suffix index cache. (Requires outcome of step 2)
         * We store the sa-value in intervals of 32 values.
         */
        bwt_cal_sa_step3( 32 );
    } // method

    /** Initializes the array _aCountTable.
     */
    void _vInitCountTable( )
    {
        for( unsigned int i = 0; i != 256; ++i )
        {
            uint32_t x = 0; // 4 byte x[3],x[2],x[1],x[0]

            /* Count the occurrences of the 4 2-bit symbols in i
             */
            for( unsigned int j = 0; j != 4; ++j )
            {
                x |= ( ( ( i & 3 ) == j ) + ( ( i >> 2 & 3 ) == j ) + ( ( i >> 4 & 3 ) == j ) + ( i >> 6 == j ) )
                     << ( j << 3 ); // (j << 3) is 0, 8, 16, 24
            } // for

            /* x[c] = number of c in i, where i is decomposed in groups of 2 bit and c in {0, 1, 2,
             * 3} Example: i=$FF = 11111111b, so x[0] = 0, x[1] = 0, x[2] = 0, x[3] = 4
             */
            _aCountTable[ i ] = x;
        } // for
    } // method

    /** Counts A,C,G,T (00, 01, 10, 11) in ux16symbols .
     * ux16symbols represents 16 nucleotides encoded using 2 bit each.
     * The returned uint32_t represents 4 unsigned 8 bit values.
     */
    inline uint32_t __occ_aux4( uint32_t ux16symbols )
    {
        return _aCountTable[ (ux16symbols)&0xff ] +
               _aCountTable[ ( ux16symbols ) >> 8 & 0xff ] // parallel add    of 4 bytes
               + _aCountTable[ ( ux16symbols ) >> 16 & 0xff ] // parallel add of 4 bytes
               + _aCountTable[ ( ux16symbols ) >> 24 ]; // parallel add of 4 bytes
    } // method

    /** The following two lines are ONLY correct when OCC_INTERVAL==0x80 (OCC_INTERVAL == 128, 2^7)
     * Delivers pointer to the first 32 bit-block of interval.
     * FIX ME: In the context of the BWT as vector this is not nice at all. (Work with references
     * over here)
     */
    inline uint32_t* bwt_occ_intv( t_bwtIndex k )
    {
        return &bwt[ 0 ] + ( ( k >> 7 ) << 4 ); // (k / 128) * 16  (we take 16, because we need twice the size)
    } // method

    /** k (input) is the position of some nucleotide within the BWT.
     * Represents count function within the BWT.
     * cnt[4] (output) delivers the number of A's, C's, G's and T's that are in front of nucleotide
     * number k. Special case k = -1: initialize cnt[4] with 0. See
     * http://en.wikipedia.org/wiki/FM-index for more info. (However, counting for k seems to start
     * with 0 instead of 1!)
     */
    inline void bwt_occ4( t_bwtIndex k, // position within the BWT (input)
                          bwt64bitCounter cnt[ 4 ] // output ( 4 counter )
    )
    {
        uint64_t x; // originally 64 bit, but should work with a 32 bit value as well.
        uint32_t *p, tmp, *end; // pointer into bwt.

        if( k == ( t_bwtIndex )( -1 ) )
        {
            std::memset( cnt, 0,
                         4 * sizeof( bwt64bitCounter ) ); // TO DO: Use C++ sizeof for arrays
            return;
        } // if

        /* Adjust k, because $ is not in bwt
         */
        k -= ( k >= primary );
        DEBUG_3( std::cout << "primary: " << ( k >= primary ) << " " << primary << std::endl; )

        /* Step 1.
         * Retrieve Occ at k/OCC_INTERVAL. (These are accumulated data, precomputed in the context
         * of the construction process) p points to a counter block start.
         */
        p = bwt_occ_intv( k );
        std::memcpy( cnt, p, 4 * sizeof( bwt64bitCounter ) ); // TO DO: Better C++ sizeof for arrays

        DEBUG_3(
            // std::cout << "bwt_size : " << bwt_size << " k: " << k << " " << (k < bwt_size) <<
            // "\n";
            std::cout << "cnt after copy " << cnt[ 0 ] << " " << cnt[ 1 ] << " " << cnt[ 2 ] << " " << cnt[ 3 ]
                      << "\n"; )

        /* Step 2. Jump (4(==sizeof(uint_32)) * 8 bytes(==sizeof(uint_64))) to the start of the
         * first BWT cell (nucleotide) (see bwt_bwtupdate_core) This line is dirty, works only if
         * the counter a 8 byte values.
         */
        p += sizeof( bwt64bitCounter );

        /* Pointer to the end (with respect to k) within the group
         */
        end = p + ( ( k >> 4 ) - ( ( k & ~OCC_INTV_MASK ) >> 4 ) ); // this is the end point of the following loop

        /* Add up nucleotides within the blocks up to k, using the parallel addition approach.
         * At the end, we have in x the counting for A, C, G, T
         */
        for( x = 0; p < end; ++p )
        {
            x += __occ_aux4( *p );
        } // for

        /* Inspect the 16 nucleotide block comprising nucleotide number k.
         */
        tmp = *p & ~( ( 1U << ( ( ~k & 15 ) << 1 ) ) - 1 );
        x += __occ_aux4( tmp ) - ( ~k & 15 );

        /* Add up the element counted for the final 128 nucleotides block
         */
        cnt[ 0 ] += x & 0xff; // A
        cnt[ 1 ] += x >> 8 & 0xff; // C
        cnt[ 2 ] += x >> 16 & 0xff; // G
        cnt[ 3 ] += x >> 24 & 0xff; // T

        DEBUG_3( std::cout << "GOT for k: " << k << " => " << cnt[ 0 ] << " " << cnt[ 1 ] << " " << cnt[ 2 ] << " "
                           << cnt[ 3 ] << "\n"; )
    } // method

    /** Writes the BWT of the FM-Index to the stream given as argument.
     * WARNING: The given stream should be opened with flag ios::binary.
     */
    void vSaveBWT( std::ostream& rxOutputStream )
    { /* Write the primary, 4 elements of L2 array as well as the BWT itself to the file system.
       */
        rxOutputStream.write( reinterpret_cast<const char*>( &primary ),
                              sizeof( primary ) ); // BWA code: err_fwrite( &bwt->primary, sizeof(bwtint_t), 1, fp );
        rxOutputStream.write( reinterpret_cast<const char*>( &L2[ 1 ] ),
                              sizeof( L2[ 0 ] ) * 4 ); //  BWA code: err_fwrite( bwt->L2 + 1, sizeof(bwtint_t), 4, fp );
        /* Write the vector to the file system.
         */
        rxOutputStream.write( reinterpret_cast<const char*>( &bwt[ 0 ] ), sizeof( bwt[ 0 ] ) * bwt.size( ) );
        rxOutputStream.flush( );
    } // method

    /** Writes the Suffix Array to the given stream.
     */
    void vSaveSuffixArray( std::ostream& rxOutputStream )
    {
        rxOutputStream.write( reinterpret_cast<const char*>( &primary ),
                              sizeof( primary ) ); // BWA code: err_fwrite( &bwt->primary, sizeof(bwtint_t), 1, fp );
        rxOutputStream.write( reinterpret_cast<const char*>( &L2[ 1 ] ),
                              sizeof( L2[ 0 ] ) * 4 ); //  BWA code: err_fwrite( bwt->L2 + 1, sizeof(bwtint_t), 4, fp );
        rxOutputStream.write( reinterpret_cast<const char*>( &sa_intv ),
                              sizeof( sa_intv ) ); // BWA code: err_fwrite( &bwt->sa_intv, sizeof(bwtint_t), 1, fp );
        rxOutputStream.write( reinterpret_cast<const char*>( &uiRefSeqLength ),
                              sizeof( uiRefSeqLength ) ); // BWA code: err_fwrite( &bwt->seq_len,
                                                          // sizeof(bwtint_t), 1, fp );
        /* Write the suffix array to an file. rxOutputStream must be in binary mode!
         * Taken from:
         * http://stackoverflow.com/questions/12372531/reading-and-writing-a-stdvector-into-a-file-correctly
         */
        rxOutputStream.write( reinterpret_cast<const char*>( &sa[ 1 ] ), ( sa.size( ) - 1 ) * sizeof( bwtint_t ) );
        //// The following apporach does not work: Why? (we have no iterator for uint_64) std::copy(
        /// sa.begin() + 1, sa.end(), std::ostreambuf_iterator<char>( rxOutputStream ) );
        rxOutputStream.flush( );
    } // method

    /** Reads the BWT from the given stream.
     * IMPROVEMENT: The original BWA design for BWT storage is quite crappy. Store the size of the
     * BWT as part of the data.
     */
    void vRestoreBWT( std::ifstream& rxInputStream )
    {
        /* Get the file size.
         * FIX ME: Move this code to a function of its own.
         */
        rxInputStream.seekg( 0, std::ios::end ); // go to the end of stream
        auto xFileEndPos = rxInputStream.tellg( ); // read position
        rxInputStream.seekg( 0, std::ios::beg ); // go to the begin of stream
        auto xFileStartPos = rxInputStream.tellg( ); // read position
        auto uiStreamLength = ( xFileEndPos - xFileStartPos ); // get the size of the stream

        if( rxInputStream.fail( ) ) // check whether we could successfully open the stream
        {
            throw std::runtime_error( "Opening of BWT of FM index failed: cannot seek." );
        } // if

        rxInputStream.seekg( 0, std::ios::beg ); // go back to the beginning of the stream.
        rxInputStream.read( reinterpret_cast<char*>( &primary ),
                            sizeof( primary ) ); // get the primary back into the data structure.
        rxInputStream.read( reinterpret_cast<char*>( &L2[ 1 ] ),
                            sizeof( L2[ 0 ] ) * 4 ); // get the L2 values back from the filesystem

        /* Compute the bwt size using the info about the stream size.
         */
        size_t uiExpectedBWTSize =
            ( uiStreamLength - ( sizeof( primary ) + ( sizeof( L2[ 0 ] ) * 4 ) ) ) / sizeof( bwt[ 0 ] );

        /* This resize is quite expensive, because all elements are initialized with 0.
         */
        bwt.resize( uiExpectedBWTSize );

        /* Restore vector by reading from file system.
         */
        rxInputStream.read( reinterpret_cast<char*>( &bwt[ 0 ] ), sizeof( bwt[ 0 ] ) * uiExpectedBWTSize );

        if( rxInputStream.fail( ) )
        { /* We should have no bad stream over here, or the input file did not comprise the expected
           * number of bytes.
           */
            throw std::runtime_error( std::string( "Unexpected fail after reading BWT from stream. " ) );
        } // if

        uiRefSeqLength = L2[ 4 ]; // restore seq_len;
        L2[ 0 ] = 0; // clean L2[0]
    } // method

    /** Restore suffix array by using the given input stream.
     */
    void vRestoreSuffixArray( std::istream& rxInputStream )
    {
        decltype( primary ) uiPrimaryReferenceValue;
        decltype( uiRefSeqLength ) uiSeqLenReferenceValue;

        rxInputStream.read( reinterpret_cast<char*>( &uiPrimaryReferenceValue ),
                            sizeof( primary ) ); // get the primary back into the data structure.
        if( primary != uiPrimaryReferenceValue )
        {
            throw std::runtime_error( "BWT and suffix array have different primary." );
        } // if

        char aSkipped[ sizeof( L2[ 0 ] ) * 4 ]; // small buffer for dummy reading
        rxInputStream.read( reinterpret_cast<char*>( aSkipped ),
                            sizeof( L2[ 0 ] ) * 4 ); // get the L2 values but ignore them.
        rxInputStream.read( reinterpret_cast<char*>( &sa_intv ),
                            sizeof( sa_intv ) ); // read suffix array interval size.
        rxInputStream.read( reinterpret_cast<char*>( &uiSeqLenReferenceValue ),
                            sizeof( uiSeqLenReferenceValue ) ); // read the sequence length for verification purposes
        if( this->uiRefSeqLength != uiSeqLenReferenceValue )
        {
            throw std::runtime_error( "SA-BWT inconsistency: suffix array has non matching sequence length stored." );
        } // if

        /* Take from:
         * http://stackoverflow.com/questions/15138353/reading-the-binary-file-into-the-vector-of-unsigned-chars
         * VERY IMPORTANT: Stop eating new lines in binary mode. (required for the back_inserter.)
         */
        rxInputStream.unsetf( std::ios::skipws );

        /* Compute the expected suffix array size and reserve appropriate space with the vector.
         */
        size_t uiExpectedSASize = ( uiRefSeqLength + sa_intv ) / sa_intv;

        /* This resize is quite expensive, because all elements are initialized with 0.
         */
        sa.resize( uiExpectedSASize );

        /* Read the complete remaining file data
         */
        sa[ 0 ] = -1; // special marker value ($FFFFFFFFFFFFFFFF with uint_64t)

        /* Read all suffix array values except the first one, which is a always a special one
         */
        rxInputStream.read( reinterpret_cast<char*>( &sa[ 1 ] ), sizeof( uint64_t ) * ( uiExpectedSASize - 1 ) );

        if( rxInputStream.fail( ) )
        { /* We should have no bad stream over here, or the input file did not comprise the expected
           * number of bytes.
           */
            throw std::runtime_error( std::string( "Unexpected bad after reading suffix array from stream. " ) );
        } // if

        /* Check whether we read exactly the required number of elements.
         */
        if( sa.size( ) != uiExpectedSASize )
        {
            throw std::runtime_error( std::string( "Reading suffix array from file system failed due to non matching "
                                                   "expected size." ) );
        } // if
    } // method

  public:
    /** Like bwt_occ4, but it computes the counters for two indices simultaneously.
     * If k and l belong to the same internal interval (the FM index has a block structure),
     * bwt_2occ4 is more efficient than bwt_occ4. IMPORTANT: The indices are inclusive, i. e. we
     * count IMPORTANT: Requires k <= l
     */
    inline void bwt_2occ4( t_bwtIndex k, // first end index in BWT (k can be negative: (t_bwtIndex)(-1) )
                           t_bwtIndex l, // second end index in BWT (l can be negative: (t_bwtIndex)(-1) )
                           bwt64bitCounter cntk[ 4 ], bwt64bitCounter cntl[ 4 ] // (outputs) counter for k and l
    )
    { /* Adjusted versions of k and l, because $ is not in bwt
       */
        t_bwtIndex _k = k - ( k >= primary );
        t_bwtIndex _l = l - ( l >= primary );

        if( _l >> OCC_INTV_SHIFT != _k >> OCC_INTV_SHIFT // l and k belong to different intervals
            || k == ( t_bwtIndex )( -1 ) //
            || l == ( t_bwtIndex )( -1 ) //
            || true )
        {
            /* Interval decomposition step by two calls of bwt_occ4
             */
            //// std::cout << "2 calls of bwt_occ4 with " << " k:" << k << " l:" << l << "\n";
            bwt_occ4( k, cntk );
            bwt_occ4( l, cntl );
        } // if
        else
        { /* The computation of bwt_occ4 but now for two indices in parallel.
           */
            uint64_t x,
                y; // originally already 64 bit, but should work with a 32 bit value as well.
            uint32_t *p, tmp, *endk, *endl;

            /* Adjust k and l, because $ is not in bwt
             */
            k -= ( k >= primary );
            l -= ( l >= primary );

            /* Step 1. (see bwt_occ4)
             */
            p = bwt_occ_intv( k );

            DEBUG_3( std::cout << "k is " << k << " p is " << p << "\n"; )
            memcpy( cntk, p, 4 * sizeof( bwt64bitCounter ) );

            DEBUG_3( for( int i = 0; i < 4; ++i ) std::cout << "cntk[" << i << "]=" << cntk[ i ] << "  ";
                     std::cout << "\n"; )

            /* Step 2. (see bwt_occ4)
             */
            p += sizeof( bwt64bitCounter ); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
            // prepare cntk[]
            endk = p + ( ( k >> 4 ) - ( ( k & ~OCC_INTV_MASK ) >> 4 ) );
            endl = p + ( ( l >> 4 ) - ( ( l & ~OCC_INTV_MASK ) >> 4 ) );

            DEBUG_3( std::cout << "endk - p:" << endk - p << " endl - p:" << endl - p << "\n"; )

            for( x = 0; p < endk; ++p )
            {
                x += __occ_aux4( *p );
            } // for

            /* Correct due to precondition k <= l
             */
            y = x;
            tmp = *p & ~( ( 1U << ( ( ~k & 15 ) << 1 ) ) - 1 );
            x += __occ_aux4( tmp ) - ( ~k & 15 );

            /* Continue summing up for l
             */
            for( ; p < endl; ++p )
            {
                y += __occ_aux4( *p );
            } // for
            tmp = *p & ~( ( 1U << ( ( ~l & 15 ) << 1 ) ) - 1 );
            y += __occ_aux4( tmp ) - ( ~l & 15 );

            memcpy( cntl, cntk,
                    4 * sizeof( bwt64bitCounter ) ); // could be moved to step 1, without changing semantics
            cntk[ 0 ] += x & 0xff;
            cntl[ 0 ] += y & 0xff; // A
            cntk[ 1 ] += x >> 8 & 0xff;
            cntl[ 1 ] += y >> 8 & 0xff; // C
            cntk[ 2 ] += x >> 16 & 0xff;
            cntl[ 2 ] += y >> 16 & 0xff; // G
            cntk[ 3 ] += x >> 24 & 0xff;
            cntl[ 3 ] += y >> 24 & 0xff; // T

        } // else
    } // method

    /**
     * @brief perform a backwards extension
     * @details
     * perform a backwards extension with the nucleotide c and the SAInterval ik
     * this also updates the position of the reverse complement interval
     */
    SAInterval EXPORTED extend_backward(
        // current interval
        const SAInterval& ik,
        // the character to extend with
        const uint8_t c );

    SAInterval init_interval(
        // the character to init with
        const uint8_t c )
    {
        return SAInterval( this->L2[ c ] + 1,
                           this->L2[ (int)NucSeq::nucleotideComplement( c ) ] + 1,
                           this->L2[ (int)c + 1 ] - this->L2[ (int)c ] );
    } // method

    /** We keep the reference length private in order to avoid unexpected trouble.
     * Delivers the length of the reference (pack) that belongs to the current FM-index.
     */
    uint64_t getRefSeqLength( void ) const
    {
        return uiRefSeqLength;
    } // method

    /** Delivers the Position in the reference .uiDeltasequence T that belongs to the position k in
     * the BWT. Uses the suffix array cache ...
     */
    bwtint_t bwt_sa( bwtint_t uiBWTposition )
    { /* Check uiBWTposition for out of range
       */
        assert( ( uiBWTposition >= 0 ) &&
                ( uiBWTposition <= static_cast<bwtint_t>( uiRefSeqLength ) ) ); // Out of range check for uiBWTposition

        bwtint_t sa = 0;
        const bwtint_t mask = sa_intv - 1; // this is why sa_intv must be a power of 2

        /* We iterate so long until we meet an interval value of the cache.
         */
        while( uiBWTposition & mask )
        {
            ++sa;
            uiBWTposition = bwt_invPsi( uiBWTposition );
        } // while

        /* Without setting bwt->sa[0] = -1 in bwt_cal_sa_step3, the following line should be
         * changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1)
         */
        auto uiReturnedPosition = sa + this->sa[ uiBWTposition / sa_intv ];
        assert( ( uiReturnedPosition >= 0 ) &&
                ( uiReturnedPosition <
                  static_cast<bwtint_t>( uiRefSeqLength ) ) ); // Out of range for returned reference position

        return uiReturnedPosition;
    } // method

    /** Checks whether the files required for loading a pack does exist on the file system.
     */
    static bool packExistsOnFileSystem( const std::string& rsPrefix )
    {
        // 1 == use std o check for file existance
        return fileExists( std::string( rsPrefix ).append( ".bwt" ) ) &&
               fileExists( std::string( rsPrefix ).append( ".sa" ) );
    } // method

    /** Dump the current FM-Index to two separated files for BWT and SA.
     */
    void vStoreFMIndex_boost( const std::string& rxFileNamePrefix )
    {
        { /* Save burrow wheeler transform
           */
            std::ofstream xOutputStream( rxFileNamePrefix + ".bwt",
                                         std::ios::binary | std::ios::trunc | std::ios::out );
            vSaveBWT( xOutputStream );
            xOutputStream.close( );
        } // scope

        { /* Save suffix array
           */
            std::ofstream xOutputStream( rxFileNamePrefix + ".sa", std::ios::binary | std::ios::trunc | std::ios::out );
            vSaveSuffixArray( xOutputStream );
            xOutputStream.close( );
        } // scope
    } // method

    /** wrap the vStoreFMIndex function in oder to make it acessible to pyhton */
    void vStoreFMIndex( const char* sPrefix )
    {
        std::string sPath( sPrefix );
        vStoreFMIndex_boost( sPath );
    } // function

    /** Load an FM-Index previously stored by vStoreFMIndex.
     */
    void vLoadFMIndex_boost( const std::string& rxFileNamePrefix )
    {
        {
            if( !fileExists( rxFileNamePrefix + ".bwt" ) )
            {
                throw std::runtime_error( "File opening error: " + rxFileNamePrefix + ".bwt" );
            } // if

            std::ifstream xInputFiletream( rxFileNamePrefix + ".bwt", std::ios::binary | std::ios::in );
            if( xInputFiletream.fail( ) ) // check whether we could successfully open the stream
            {
                throw std::runtime_error( "Opening of BWT of FM index failed: " + rxFileNamePrefix + ".bwt" );
            } // if

            vRestoreBWT( xInputFiletream );
            xInputFiletream.close( );
        } // scope

        /* The order is important over here, because during loading the suffix array we check with
         * compatibility with the BWT.
         */
        {
            std::ifstream xInputFileStream( rxFileNamePrefix + ".sa", std::ios::binary | std::ios::in );
            if( xInputFileStream.fail( ) ) // check whether we could successfully open the stream
            {
                throw std::runtime_error( "Opening of suffix array of FM index failed: " + rxFileNamePrefix + ".sa" );
            } // if
            vRestoreSuffixArray( xInputFileStream );
            xInputFileStream.close( );
        } // scope
    } // method

    void vLoadFMIndex( std::string sPrefix )
    {
        std::string sPath( sPrefix );
        vLoadFMIndex_boost( sPath );
    } // function

    /* Debug function for comparing BWT.
     * BWT can be build by several different functions. Here we can check for correctness.
     * So long it does not compare the BWT sequence itself.
     */
    bool operator==( const FMIndex& rxOtherFMIndex )
    {
        std::string sErrorText = ( L2 != rxOtherFMIndex.L2 )
                                     ? "Different L2"
                                     : ( primary != rxOtherFMIndex.primary )
                                           ? "Different primary"
                                           : ( uiRefSeqLength != rxOtherFMIndex.uiRefSeqLength )
                                                 ? "Different seq_len"
                                                 : ( bwt != rxOtherFMIndex.bwt )
                                                       ? "Different bwt_size"
                                                       : ( sa != rxOtherFMIndex.sa ) ? "Different suffix arrays" : "";
        if( sErrorText != "" )
        {
            std::cerr << "BWT different: " << sErrorText << std::endl;
        } // if

        return sErrorText == "";
    } // method

    // overload
    bool canCast( std::shared_ptr<Container> c ) const
    {
        return std::dynamic_pointer_cast<FMIndex>( c ) != nullptr;
    } // function

    // overload
    std::string getTypeName( ) const
    {
        return "FMIndex";
    } // function

    // overload
    std::shared_ptr<Container> getType( ) const
    {
        return std::shared_ptr<Container>( new FMIndex( ) );
    } // function

    /* Default constructor. (Initializes the fix count-table)
     */
    FMIndex( )
        : L2( {{0, 0, 0, 0, 0, 0}} ),
          bwt( ), // initialize the vector keeping the BWT
          sa( ) // initialize the vector keeping the suffix array
    {
        _vInitCountTable( );
    } // default constructor

    /* FM-Index constructor. Builds a FM index on foundation of rxSequence.
     */
    FMIndex( const NucSeq& rxSequence ) : FMIndex( ) // call the default constructor
    {
        build_FMIndex( rxSequence );
    } // constructor

    /* FM-Index constructor. Loads a fm index.
     */
    FMIndex( const std::string sFileName ) : FMIndex( ) // call the default constructor
    {
        vLoadFMIndex( sFileName );
    } // constructor

    /* FM-Index constructor. Builds a FM index on foundation of pxSequence.
     */
    FMIndex( const std::shared_ptr<NucSeq> pxSequence ) : FMIndex( ) // call the default constructor
    {
        build_FMIndex( *pxSequence );
    } // constructor

    /* FM-Index constructor. Builds a FM index on foundation of a given sequence collection.
     */
    FMIndex( const Pack& rxSequenceCollection, // the pack for which we require a BWT
             unsigned int uiAlgorithmSelection = 2 // 2 -> automatic algorithm selection
             )
        : FMIndex( ) // call the default constructor
    {
        build_FMIndex( rxSequenceCollection, uiAlgorithmSelection );
        DEBUG( if( !test( rxSequenceCollection, 10000 ) ) // test 10000 sequences
               {
                   std::cerr << "WARNING: suffix array test failed" << std::endl;
                   exit( 0 );
               } // if
        )
    } // constructor

    /* FM-Index constructor. Builds a FM index on foundation of a given sequence collection.
     */
    FMIndex(
        // the pack for which we require a BWT
        const std::shared_ptr<Pack>
            pxSequenceCollection )
        : FMIndex( ) // call the default constructor
    {
        build_FMIndex( *pxSequenceCollection, 2 );
        DEBUG( if( !test( *pxSequenceCollection, 10000 ) ) // test 10000 sequences
               {
                   std::cerr << "WARNING: suffix array test failed" << std::endl;
                   exit( 0 );
               } // if
        )
    } // constructor
}; // class FMIndex

} // namespace libMA


#ifdef WITH_PYTHON
/**
 * @brief function called in order to export this @ref Module "module"
 * @ingroup export
 */
void exportFM_index( py::module& rxPyModuleId );
#endif

#endif