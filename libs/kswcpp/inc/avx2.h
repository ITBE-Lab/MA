#pragma once

/* Include AVX2 header */
#if defined( __GNUC__ )
#include <x86intrin.h>
#else
#include <intrin.h>
#endif

#include <iostream>
#include <stdint.h>

#if defined Vc_CLANG || defined Vc_APPLECLANG
#define Vc_INTRINSIC_L inline
#define Vc_INTRINSIC_R __attribute__( ( always_inline ) )
#define Vc_INTRINSIC Vc_INTRINSIC_L Vc_INTRINSIC_R
#elif defined( __GNUC__ )
#define Vc_INTRINSIC_L inline
#define Vc_INTRINSIC_R __attribute__( ( __always_inline__, __artificial__ ) )
#define Vc_INTRINSIC Vc_INTRINSIC_L Vc_INTRINSIC_R
#elif defined( _MSC_VER )
#define Vc_INTRINSIC inline __forceinline
#else
#error "Unknown compiler. Don't know how to force inline ..."
#endif

/* TP_SELF must be a derived class of this base class.
 * Technique: https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
 */
template <typename TP_SELF> struct alignas( 32 ) Vec_m256i_base
{
    __m256i val;

    Vec_m256i_base( __m256i val )
    {
        this->val = val;
    } // constructor


    Vec_m256i_base( void )
    {} // default constructor


    Vc_INTRINSIC int8_t* pointerVal( void )
    {
        return static_cast<int8_t*>( static_cast<void*>( &this->val ) );
    } // method


    /* Reading via pure address related mapping.
     * WARNING: The physical organization of SIMD-vectors is different to the logic one.
     */
    Vc_INTRINSIC int8_t at( const size_t iIdx )
    {
        return this->pointerVal( )[ iIdx ];
    } // method


    void set( const size_t iIdx, const int8_t iVal )
    {
        this->pointerVal( )[ iIdx ] = iVal;
    } // method


    Vc_INTRINSIC const __m256i& unwrap( const __m256i& rhs ) const
    {
        return rhs;
    } // method


    Vc_INTRINSIC const __m256i& unwrap( const TP_SELF& rhs ) const
    {
        return rhs.val;
    } // method


    /* Assignment for derived types. */
    Vc_INTRINSIC TP_SELF& operator=( const TP_SELF& rhs )
    {
        /* rhs must wrap an aligned value, so a good compiler should go for an aligned assignment.
         * Assignment is equal to: _mm_store_si128( &(this->val), rhs.val );
         */
        this->val = rhs.val;
        return static_cast<TP_SELF&>( *this );
    } // method


    /* Set all bits of vector to zero */
    Vc_INTRINSIC TP_SELF& setzero( void )
    {
        this->val = _mm256_setzero_si256( );
        return static_cast<TP_SELF&>( *this );
    } // method


    Vc_INTRINSIC TP_SELF& moveAndZero32( int32_t iVal )
    {
        this->val = _mm256_insert_epi32( _mm256_setzero_si256( ), iVal, 0 );
        return static_cast<TP_SELF&>( *this );
    } // method


    /* Store this->value to unaligned address in memory */
    Vc_INTRINSIC const TP_SELF& store( void* addr ) const
    {
        _mm256_storeu_si256( static_cast<__m256i*>( addr ), this->val );
        return static_cast<const TP_SELF&>( *this );
    } // method


    /* Load this->val from unaligned address in memory */
    Vc_INTRINSIC const TP_SELF& load( void* addr )
    {
        this->val = _mm256_loadu_si256( static_cast<__m256i*>( addr ) );
        return static_cast<const TP_SELF&>( *this );
    } // method


    /* Number of bits that are set (value 1) in the 256 bit vector.
     * See: https://github.com/WojciechMula/sse-popcount
     */
    Vc_INTRINSIC size_t popcnt( void )
    {
#ifdef _MSC_VER
        return __popcnt64( _mm256_extract_epi64( this->val, 0 ) ) + __popcnt64( _mm256_extract_epi64( this->val, 1 ) ) +
               __popcnt64( _mm256_extract_epi64( this->val, 2 ) ) + __popcnt64( _mm256_extract_epi64( this->val, 3 ) );
#elif defined( __GNUC__ )
        return __builtin_popcountll( _mm256_extract_epi64( this->val, 0 ) ) +
               __builtin_popcountll( _mm256_extract_epi64( this->val, 1 ) ) +
               __builtin_popcountll( _mm256_extract_epi64( this->val, 2 ) ) +
               __builtin_popcountll( _mm256_extract_epi64( this->val, 3 ) );
#else
#error "No popcnt for this compiler so far ..."
#endif
    } // method


    /* Bitwise AND over the complete 256 bit vector */
    template <typename TP> Vc_INTRINSIC TP_SELF operator&( const TP& operand ) const
    {
        return TP_SELF( _mm256_and_si256( this->val, unwrap( operand ) ) );
    } // method


    /* Bitwise OR over the complete 256 bit vector */
    template <typename TP> Vc_INTRINSIC TP_SELF operator|( const TP& operand ) const
    {
        return TP_SELF( _mm256_or_si256( this->val, unwrap( operand ) ) );
    } // method


    /* Bitwise combination of and and complement over the 256 bit vector */
    template <typename TP> Vc_INTRINSIC TP_SELF andnot( const TP& operand ) const
    {
        return TP_SELF( _mm256_andnot_si256( this->val, unwrap( operand ) ) );
    } // method


    /* Byte-wise left shifting with respect to the complete vector.
     * https://github.com/blegal/AVX2_shift_and_rotate_si256/blob/master/src/avx2_func.hpp
     * If you use this with GCC, you need additionally a template indication!
     * See:
     * https://stackoverflow.com/questions/32300441/error-when-calling-an-integral-template-member-function-with-g-and-clang
     */
    template <unsigned int SHIFT> Vc_INTRINSIC TP_SELF shiftl( void ) const
    {
        if( SHIFT == 0 )
            return static_cast<const TP_SELF&>( *this );
        else if( SHIFT < 16 )
            return TP_SELF( _mm256_alignr_epi8(
                this->val, _mm256_permute2x128_si256( this->val, this->val, _MM_SHUFFLE( 0, 0, 2, 0 ) ),
                ( uint8_t )( 16 - SHIFT ) ) );
        else if( SHIFT == 16 )
            return TP_SELF( _mm256_permute2x128_si256( this->val, this->val, _MM_SHUFFLE( 0, 0, 2, 0 ) ) );
        else
            return TP_SELF(
                _mm256_slli_si256( _mm256_permute2x128_si256( this->val, this->val, _MM_SHUFFLE( 0, 0, 2, 0 ) ),
                                   ( uint8_t )( SHIFT - 16 ) ) );
    } // method


    /* Byte-wise right shifting with respect to the complete 256-bit vector * If you use this with GCC, you need
     * additionally a template indication! See:
     * https://stackoverflow.com/questions/32300441/error-when-calling-an-integral-template-member-function-with-g-and-clang
     */
    template <unsigned int SHIFT> Vc_INTRINSIC TP_SELF shiftr( void ) const
    {
        if( SHIFT == 0 )
            return static_cast<const TP_SELF&>( *this );
        else if( SHIFT < 16 )
            return TP_SELF(
                _mm256_alignr_epi8( _mm256_permute2x128_si256( this->val, this->val, _MM_SHUFFLE( 2, 0, 0, 1 ) ),
                                    this->val, ( uint8_t )( SHIFT ) ) );
        else if( SHIFT == 16 )
            return TP_SELF( _mm256_permute2x128_si256( this->val, this->val, _MM_SHUFFLE( 2, 0, 0, 1 ) ) );
        else
            return TP_SELF(
                _mm256_srli_si256( _mm256_permute2x128_si256( this->val, this->val, _MM_SHUFFLE( 2, 0, 0, 1 ) ),
                                   ( uint8_t )( SHIFT - 16 ) ) );
    } // method


    /* Vertical blend using a mask
     * According header: extern __m256i __cdecl _mm256_blendv_epi8(__m256i, __m256i, __m256i);
     */
    template <typename TP_VEC, typename TP_MASK>
    Vc_INTRINSIC TP_SELF blendv( const TP_VEC& vec, const TP_MASK& mask ) const
    {
        // return struct__m256i( _mm256_blendv_epi8( this->val, unwrap( vec ), _mm256_movemask_epi8( unwrap( mask ) ) )
        // );
        return TP_SELF( _mm256_blendv_epi8( this->val, unwrap( vec ), unwrap( mask ) ) );
    } // method
}; // struct


/* Vector of 32 signed 8 bit (int8_t) values */
struct alignas( 32 ) Vec_m256i_32_int8_t : public Vec_m256i_base<Vec_m256i_32_int8_t>
{
    static constexpr size_t SIZE = 32; // 32 elements of type int8_t

    Vec_m256i_32_int8_t( __m256i val ) : Vec_m256i_base<Vec_m256i_32_int8_t>( val )
    {} // constructor


    /* Default constructor does not initialize the vector! */
    Vec_m256i_32_int8_t( void )
    {} // default constructor

    /* Construct vector and initialize all elements with iVal */
    Vec_m256i_32_int8_t( int8_t iVal )
    {
        this->set1( iVal );
    } // constructor


    Vc_INTRINSIC Vec_m256i_32_int8_t& set1( int8_t iVal )
    {
        this->val = _mm256_set1_epi8( iVal );
        return static_cast<Vec_m256i_32_int8_t&>( *this );
    } // method


    /* Vertical maximum with respect to the argument */
    template <typename TP_VEC> Vc_INTRINSIC Vec_m256i_32_int8_t max( const TP_VEC& vec ) const
    {
        /* REQUIRES SSE 4.1*/
        return Vec_m256i_32_int8_t( _mm256_max_epi8( this->val, unwrap( vec ) ) );
    } // method


    /* Vertical maximum with respect to the argument */
    template <typename TP_VEC> Vc_INTRINSIC Vec_m256i_32_int8_t min( const TP_VEC& vec ) const
    {
        /* REQUIRES SSE 4.1*/
        return Vec_m256i_32_int8_t( _mm256_min_epi8( this->val, unwrap( vec ) ) );
    } // method


    /* Parallel vertical addition */
    template <typename TP> Vc_INTRINSIC Vec_m256i_32_int8_t operator+( const TP& operand ) const
    {
        return Vec_m256i_32_int8_t( _mm256_add_epi8( this->val, unwrap( operand ) ) );
    } // method


    /* Parallel vertical subtraction */
    template <typename TP> Vc_INTRINSIC Vec_m256i_32_int8_t operator-( const TP& operand ) const
    {
        return Vec_m256i_32_int8_t( _mm256_sub_epi8( this->val, unwrap( operand ) ) );
    } // method


    /* Parallel Comparison (greater than) of all vector elements */
    template <typename TP> Vc_INTRINSIC Vec_m256i_32_int8_t operator>( const TP& operand ) const
    {
        return Vec_m256i_32_int8_t( _mm256_cmpgt_epi8( this->val, unwrap( operand ) ) );
    } // method


    /* Parallel Comparison (equality) of all vector elements */
    template <typename TP> Vc_INTRINSIC Vec_m256i_32_int8_t operator==( const TP& operand ) const
    {
        return Vec_m256i_32_int8_t( _mm256_cmpeq_epi8( this->val, unwrap( operand ) ) );
    } // method
}; // struct


/* Vector of 8 signed 32 bit values */
struct alignas( 32 ) Vec_m256i_8_int32_t : public Vec_m256i_base<Vec_m256i_8_int32_t>
{
    static constexpr size_t SIZE = 8; // 8 times 32 bit (int32_t)
    typedef int32_t TP_ELEMENT;

    Vec_m256i_8_int32_t( __m256i val ) : Vec_m256i_base<Vec_m256i_8_int32_t>( val )
    {} // constructor


    /* Default constructor does not initialize the vector! */
    Vec_m256i_8_int32_t( void )
    {} // default constructor


    void dump( void ) const
    {
        std::cout << _mm256_extract_epi32( this->val, 0 ) << ", " << _mm256_extract_epi32( this->val, 1 ) << ", "
                  << _mm256_extract_epi32( this->val, 2 ) << ", " << _mm256_extract_epi32( this->val, 3 ) << ", "
                  << _mm256_extract_epi32( this->val, 4 ) << ", " << _mm256_extract_epi32( this->val, 5 ) << ", "
                  << _mm256_extract_epi32( this->val, 6 ) << ", " << _mm256_extract_epi32( this->val, 7 ) << std::endl;
    } // method


    /* Initializes all 4 vector elements to iVal */
    Vc_INTRINSIC Vec_m256i_8_int32_t& set1( const int32_t iVal )
    {
        this->val = _mm256_set1_epi32( iVal );
        return static_cast<Vec_m256i_8_int32_t&>( *this );
    } // method


    /* Initializes the 8 vector elements with the 8 arguments */
    Vc_INTRINSIC Vec_m256i_8_int32_t& set8( const int32_t iVal0, const int32_t iVal1, const int32_t iVal2,
                                            const int32_t iVal3, const int32_t iVal4, const int32_t iVal5,
                                            const int32_t iVal6, const int32_t iVal7 )
    {
        this->val = _mm256_setr_epi32( iVal0, iVal1, iVal2, iVal3, iVal4, iVal5, iVal6, iVal7 );
        return static_cast<Vec_m256i_8_int32_t&>( *this );
    } // method


    /* Signed expansion of an array of 8 int8_t values as int32_t into 256-bit vector.
     */
    static Vc_INTRINSIC Vec_m256i_8_int32_t expand_int8_t( const int8_t* pArr )
    { /* The below code is an efficient replacement for:
       * this->val = _mm_setr_epi32( pArr[0], pArr[1], pArr[2], pArr[3], pArr[4], pArr[5], pArr[6], pArr[7] );
       * The problem with the above line are the 8 (serial and so expensive) type castings from int8_t to int32_t
       * Observation: AVX2 on AMD-Ryzen is so slow that the latter code is slower than the above single line.
       */
        // Get 64 bit (8 times int8_t) from the pointer and store them in the lower lane
        const __m256i vLoadMask = _mm256_setr_epi32( 0xFFFFFFFF, 0xFFFFFFFF, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 );
        const __m256i v64Bit = _mm256_maskload_epi64( reinterpret_cast<const long long*>( pArr ), vLoadMask );
        // Distribution of 32bit segments across the two 128 bit lanes
        const __m256i vIdx = _mm256_setr_epi32( 0, 2, 3, 4, 1, 5, 6, 7 );
        const __m256i v32Bitx2InBothLanes = _mm256_permutevar8x32_epi32( v64Bit, vIdx );
        // Shuffle the bytes in both lanes to the correct positions
        const __m256i vShuffle = _mm256_setr_epi8( 0, 0x80U, 0x80U, 0x80U, 1, 0x80U, 0x80U, 0x80U, 2, 0x80U, 0x80U,
                                                   0x80U, 3, 0x80U, 0x80U, 0x80U, 0, 0x80U, 0x80U, 0x80U, 1, 0x80U,
                                                   0x80U, 0x80U, 2, 0x80U, 0x80U, 0x80U, 3, 0x80U, 0x80U, 0x80U );
        const __m256i vShuff8Bytes = _mm256_shuffle_epi8( v32Bitx2InBothLanes, vShuffle );
        // if greater 127, it was a negative value (255 = -1, ... , 128 = - 128)
        const __m256i vCompValues = _mm256_set1_epi32( INT8_MAX );
        // v_mask = 1 for the values greater 127, all other 0
        const __m256i vMask = _mm256_cmpgt_epi32( vShuff8Bytes, vCompValues );
        // expand negative values to 32 bit size
        const __m256i vNegValues = _mm256_or_si256( vShuff8Bytes, _mm256_set1_epi32( 0xFFFFFF00 ) );

        // blend positive and negative values using the mask
        return Vec_m256i_8_int32_t( _mm256_blendv_epi8( vShuff8Bytes, vNegValues, vMask ) );

        //- // blend positive and negative values using the mask
        //- this->val = _mm256_blendv_epi8( vShuff8Bytes, vNegValues, vMask );
        //-
        //- return static_cast<Vec_m256i_8_int32_t&>(*this);
    } // method


    /* Parallel addition of eight 32 bit values */
    template <typename TP> Vc_INTRINSIC Vec_m256i_8_int32_t operator+( const TP& operand ) const
    {
        return Vec_m256i_8_int32_t( _mm256_add_epi32( this->val, unwrap( operand ) ) );
    } // method


    /* Parallel Comparison (greater than) of eight 32 bit integers */
    template <typename TP> Vc_INTRINSIC Vec_m256i_8_int32_t operator>( const TP& operand ) const
    {
        return Vec_m256i_8_int32_t( _mm256_cmpgt_epi32( this->val, unwrap( operand ) ) );
    } // method


    /* Horizontal maximum of eight 32 bit values. */
    Vc_INTRINSIC int32_t hmax( void )
    {
        // std::cout << "A:" << std::endl; dump();
        const __m256i v32Bitx2InBothLanes = _mm256_permute4x64_epi64( this->val, _MM_SHUFFLE( 1, 0, 3, 2 ) );
        const __m256i max0 = _mm256_max_epi32( this->val, v32Bitx2InBothLanes );
        // std::cout << "B:" << std::endl; Vec_m256i_8_int32_t(max0).dump();
        const __m256i max1 = _mm256_shuffle_epi32( max0, _MM_SHUFFLE( 1, 0, 3, 2 ) ); // intra-lane shuffle
        // std::cout << "Ca:" << std::endl; Vec_m256i_8_int32_t(max1).dump();
        const __m256i max2 = _mm256_max_epi32( max0, max1 );
        // std::cout << "C:" << std::endl; Vec_m256i_8_int32_t(max2).dump();
        const __m256i max3 = _mm256_shuffle_epi32( max2, _MM_SHUFFLE( 0, 0, 0, 1 ) ); // intra-lane shuffle
        const __m256i max4 = _mm256_max_epi32( max2, max3 );
        // std::cout << "D:" << std::endl; Vec_m256i_8_int32_t(max4).dump();
        return _mm256_extract_epi32( max4, 0 );
    } // method
}; // struct


/* Vector of 16 signed 16 bit values */
struct alignas( 32 ) Vec_m256i_16_int16_t : public Vec_m256i_base<Vec_m256i_16_int16_t>
{
    static constexpr size_t SIZE = 16; // 16 times 16 bit (int16_t)
    typedef int16_t TP_ELEMENT;

    Vec_m256i_16_int16_t( __m256i val ) : Vec_m256i_base<Vec_m256i_16_int16_t>( val )
    {} // constructor


    /* Default constructor does not initialize the vector! */
    Vec_m256i_16_int16_t( void )
    {} // default constructor


    /* Dump the vector for debugging purposes */
    void dump( void )
    {
        std::cout << _mm256_extract_epi16( this->val, 0 ) << ", " << _mm256_extract_epi16( this->val, 1 ) << ", "
                  << _mm256_extract_epi16( this->val, 2 ) << ", " << _mm256_extract_epi16( this->val, 3 ) << ", "
                  << _mm256_extract_epi16( this->val, 4 ) << ", " << _mm256_extract_epi16( this->val, 5 ) << ", "
                  << _mm256_extract_epi16( this->val, 6 ) << ", " << _mm256_extract_epi16( this->val, 7 ) << ","
                  << _mm256_extract_epi16( this->val, 8 ) << ", " << _mm256_extract_epi16( this->val, 9 ) << ", "
                  << _mm256_extract_epi16( this->val, 10 ) << ", " << _mm256_extract_epi16( this->val, 11 ) << ", "
                  << _mm256_extract_epi16( this->val, 12 ) << ", " << _mm256_extract_epi16( this->val, 13 ) << ", "
                  << _mm256_extract_epi16( this->val, 14 ) << ", " << _mm256_extract_epi16( this->val, 15 )
                  << std::endl;
    } // method


    /* Initializes all 4 vector elements to iVal */
    Vc_INTRINSIC Vec_m256i_16_int16_t& set1( const int16_t iVal )
    {
        this->val = _mm256_set1_epi16( iVal );
        return *this;
    } // method


    /* Signed expansion of an array of 8 int8_t values as int32_t into 256-bit vector.
     */
    static Vc_INTRINSIC Vec_m256i_16_int16_t expand_int8_t( const int8_t* pArr )
    { /* The below code is an efficient replacement (on Intel processors) for:
       * this->val = _mm256_setr_epi16( pArr[0], pArr[1], pArr[2], pArr[3], pArr[4], pArr[5], pArr[6], pArr[7],
                                        pArr[8], pArr[9], pArr[10], pArr[11], pArr[12], pArr[13], pArr[14], pArr[15]);
       * The problem with the above line are the 16 (serial and so expensive) type castings from int16_t to int32_t
       */
        // Get 64 bit (8 times int8_t) from the pointer and store them in the lower lane
        const __m256i vLoadMask =
            _mm256_setr_epi32( 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x0, 0x0, 0x0, 0x0 );
        const __m256i v128Bit = _mm256_maskload_epi64( reinterpret_cast<const long long*>( pArr ), vLoadMask );
        // Distribution of 16bit segments across the two 128 bit lanes
        const __m256i vIdx = _mm256_setr_epi32( 0, 1, 4, 5, 2, 3, 6, 7 );
        const __m256i v32Bitx2InBothLanes = _mm256_permutevar8x32_epi32( v128Bit, vIdx );
        // Shuffle the bytes in both lanes to the correct positions
        const __m256i vShuffle =
            _mm256_setr_epi8( 0, 0x80U, 1, 0x80U, 2, 0x80U, 3, 0x80U, 4, 0x80U, 5, 0x80U, 6, 0x80U, 7, 0x80U, 0, 0x80U,
                              1, 0x80U, 2, 0x80U, 3, 0x80U, 4, 0x80U, 5, 0x80U, 6, 0x80U, 7, 0x80U );
        const __m256i vShuff8Bytes = _mm256_shuffle_epi8( v32Bitx2InBothLanes, vShuffle );
        // if greater 127, it was a negative value (255 = -1, ... , 128 = - 128)
        const __m256i vCompValues = _mm256_set1_epi16( INT8_MAX );
        // vMask = 1 for the values greater 127, all other 0
        const __m256i vMask = _mm256_cmpgt_epi16( vShuff8Bytes, vCompValues );
        // expand negative values to 32 bit size
        const __m256i vNegValues = _mm256_or_si256( vShuff8Bytes, _mm256_set1_epi16( 0xFF00U ) );

        // blend positive and negative values using the mask
        return Vec_m256i_16_int16_t( _mm256_blendv_epi8( vShuff8Bytes, vNegValues, vMask ) );

        //- // blend positive and negative values using the mask
        //- this->val = _mm256_blendv_epi8( vShuff8Bytes, vNegValues, vMask );
        //-
        //- return *this;
    } // method


    /* Parallel addition of eight 32 bit values */
    template <typename TP> Vc_INTRINSIC Vec_m256i_16_int16_t operator+( const TP& operand ) const
    {
        return Vec_m256i_16_int16_t( _mm256_add_epi16( this->val, unwrap( operand ) ) );
    } // method


    /* Parallel Comparison (greater than) of eight 32 bit integers */
    template <typename TP> Vc_INTRINSIC Vec_m256i_16_int16_t operator>( const TP& operand ) const
    {
        return Vec_m256i_16_int16_t( _mm256_cmpgt_epi16( this->val, unwrap( operand ) ) );
    } // method


    /* Horizontal maximum of sixteen 16 bit values. */
    Vc_INTRINSIC int16_t hmax( void )
    {
        const __m256i v32Bitx2InBothLanes = _mm256_permute4x64_epi64( this->val, _MM_SHUFFLE( 1, 0, 3, 2 ) );
        const __m256i max0 = _mm256_max_epi16( this->val, v32Bitx2InBothLanes );
        const __m256i max1 = _mm256_shuffle_epi32( max0, _MM_SHUFFLE( 1, 0, 3, 2 ) ); // intra-lane shuffle
        const __m256i max2 = _mm256_max_epi16( max0, max1 );
        const __m256i max3 = _mm256_shuffle_epi32( max2, _MM_SHUFFLE( 0, 0, 0, 1 ) ); // intra-lane shuffle
        const __m256i max4 = _mm256_max_epi16( max2, max3 );
        const __m256i max5 = _mm256_shufflelo_epi16( max4, _MM_SHUFFLE( 0, 0, 0, 1 ) );
        const __m256i max6 = _mm256_max_epi16( max4, max5 );
        return _mm256_extract_epi16( max6, 0 );
    } // method
}; // struct