/* 128 bit SIMD vector for Intel architectures.
 https://github.com/samyvilar/dyn_perf/blob/master/sse2.h
 https://github.com/svn2github/framewave/blob/master/branch/FRAMEWAVE_1.3/Framewave/sdk/SSEPlus/include/emulation/SSEPlus_emulation_SSE2.h
 https://gitlab.eurecom.fr/oai/openairinterface5g/blob/d8ea9ea70664867b80acccf33616ec0083039933/openair1/PHY/sse_intrin.h
 */
#pragma once

/* Include appropriate header */
#include <emmintrin.h>
#include <iostream>
#include <smmintrin.h>
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
 * INSTR: SSE2: 20
 *		  SSE3: 30
 *		  SSSE3 : 31
 *		  SSE4.1: 41
 *        SSE4.2: 42
 */
template <int INSTR, typename TP_SELF> struct alignas( 16 ) Vec_m128i_base
{
#if defined( __SSE4_1__ )
    /* Avoid contradictions with with respect to defines and INSTR */
    static_assert( INSTR <= 41 );
#elif defined( __SSSE3__ )
    /* Avoid contradictions with with respect to defines and INSTR */
    static_assert( INSTR <= 31 );
#endif

    __m128i val; // wrapped 128bit vector

    Vec_m128i_base( void )
    {} // default constructor


    Vec_m128i_base( __m128i val )
    {
        this->val = val;
    } // constructor


    Vc_INTRINSIC int8_t* pointerVal( void )
    {
        return static_cast<int8_t*>( static_cast<void*>( &this->val ) );
    } // method

    /* Pure address related mapped.
     * Delivers byte at position via address mapping.
     */
    Vc_INTRINSIC int8_t at( const size_t iIdx )
    {
        return this->pointerVal( )[ iIdx ];
    } // method


    void set( const size_t iIdx, const int8_t iVal )
    {
        this->pointerVal( )[ iIdx ] = iVal;
    } // method


    /* Argument unwrapping */
    Vc_INTRINSIC const __m128i& unwrap( const __m128i& rhs ) const
    {
        return rhs;
    } // method


    /* Argument unwrapping */
    Vc_INTRINSIC const __m128i& unwrap( const TP_SELF& rhs ) const
    {
        return rhs.val;
    } // method


    /* Assignment for derived types.
     */
    Vc_INTRINSIC TP_SELF& operator=( const TP_SELF& rhs )
    {
        /* rhs must wrap an aligned value, so a good compiler should go for an aligned
         * assignment. Assignment is equal to: _mm_store_si128( &(this->val), rhs.val );
         */
        this->val = rhs.val;
        return static_cast<TP_SELF&>( *this );
    } // method


    Vc_INTRINSIC TP_SELF& setzero( void )
    {
        this->val = _mm_setzero_si128( );
        return static_cast<TP_SELF&>( *this );
    } // method


    Vc_INTRINSIC TP_SELF& moveAndZero32( int32_t iVal )
    {
        this->val = _mm_cvtsi32_si128( iVal );
        return static_cast<TP_SELF&>( *this );
    } // method


    /* See:
     * https://stackoverflow.com/questions/17354971/fast-counting-the-number-of-set-bits-in-m128i-register
     * and
     * https://github.com/WojciechMula/sse-popcount
     * Support of POPCNT is indicated via the CPUID.01H:ECX.POPCNT[Bit 23] flag
     */
    Vc_INTRINSIC size_t popcnt( void )
    {
        const __m128i n_hi = _mm_unpackhi_epi64( this->val, this->val );
#ifdef _MSC_VER
        // SSE 4.2
        return __popcnt64( _mm_cvtsi128_si64( this->val ) ) + __popcnt64( _mm_cvtsi128_si64( n_hi ) );
#elif defined( __GNUC__ )
        // For -msse4.2 GCC will use the POPCNT processor instruction
        return __builtin_popcountll( _mm_cvtsi128_si64( this->val ) ) +
               __builtin_popcountll( _mm_cvtsi128_si64( n_hi ) );
#else
#error "No hardware support for popcnt available"
#endif
    } // method


    /* Store this->value to unaligned address in memory */
    Vc_INTRINSIC const TP_SELF& store( void* addr ) const
    {
        _mm_storeu_si128( static_cast<__m128i*>( addr ), this->val );
        return static_cast<const TP_SELF&>( *this );
    } // method


    /* Load this->val from unaligned address in memory */
    Vc_INTRINSIC const TP_SELF& load( void* addr )
    {
        this->val = _mm_loadu_si128( static_cast<__m128i*>( addr ) );
        return static_cast<const TP_SELF&>( *this );
    } // method


    /* Bitwise AND over the complete 128 bit vector */
    template <typename TP> Vc_INTRINSIC TP_SELF operator&( const TP& operand ) const
    {
        return TP_SELF( _mm_and_si128( this->val, unwrap( operand ) ) );
    } // method


    /* Bitwise OR over the complete 128 bit vector */
    template <typename TP> Vc_INTRINSIC TP_SELF operator|( const TP& operand ) const
    {
        return TP_SELF( _mm_or_si128( this->val, unwrap( operand ) ) );
    } // method


    /* Bitwise combination of and and complement over the complete 128 bit vector */
    template <typename TP> Vc_INTRINSIC TP_SELF andnot( const TP& operand ) const
    {
        return TP_SELF( _mm_andnot_si128( this->val, unwrap( operand ) ) );
    } // method


    /* Byte-wise shift left.
     * If you use this with GCC, you need additionally a template indication!
     * See:
     * https://stackoverflow.com/questions/32300441/error-when-calling-an-integral-template-member-function-with-g-and-clang
     */
    template <unsigned int SHIFT> Vc_INTRINSIC TP_SELF shiftl( void ) const
    {
        return TP_SELF( _mm_slli_si128( this->val, SHIFT ) );
    } // method


    /* Byte-wise shift right
     * If you use this with GCC, you need additionally a template indication!
     * See:
     * https://stackoverflow.com/questions/32300441/error-when-calling-an-integral-template-member-function-with-g-and-clang
     */
    template <unsigned int SHIFT> Vc_INTRINSIC TP_SELF shiftr( void ) const
    {
        return TP_SELF( _mm_srli_si128( this->val, SHIFT ) );
    } // method


    /* Vertical blend using mask */
    template <typename TP_VEC, typename TP_MASK>
    Vc_INTRINSIC TP_SELF blendv( const TP_VEC& vec, const TP_MASK& mask ) const
    {
        if( INSTR >= 41 ) // compile time decision
        { /* SSE4.1 only */
#if( !defined( __GNUC__ ) || defined( __SSE4_1__ ) )
            return TP_SELF( _mm_blendv_epi8( this->val, unwrap( vec ), unwrap( mask ) ) );
#endif
        } // if
        else
        { /* Plain SSE2 */
            return TP_SELF( _mm_or_si128( _mm_and_si128( unwrap( mask ), unwrap( vec ) ),
                                          _mm_andnot_si128( unwrap( mask ), this->val ) ) );
        } // else
    } // method
}; // struct


/* SIMD Vector of 16 int8_t (signed 8 bit) values */
template <int INSTR> struct alignas( 16 ) Vec_m128i_16_int8_t : public Vec_m128i_base<INSTR, Vec_m128i_16_int8_t<INSTR>>
{
    static constexpr size_t SIZE = 16; // 16 elements of type int8_t
    typedef int8_t TP_ELEMENT;

    typedef Vec_m128i_base<INSTR, Vec_m128i_16_int8_t<INSTR>> TP_SUPER;

    Vec_m128i_16_int8_t( __m128i val ) : Vec_m128i_base<INSTR, Vec_m128i_16_int8_t<INSTR>>( val )
    {} // constructor


    /* Default constructor does not initialize the vector! */
    Vec_m128i_16_int8_t( void )
    {} // default constructor


    Vec_m128i_16_int8_t( int8_t iVal )
    {
        this->set1( iVal );
    } // constructor


    Vc_INTRINSIC Vec_m128i_16_int8_t& set1( int8_t iVal )
    {
        this->val = _mm_set1_epi8( iVal );
        return static_cast<Vec_m128i_16_int8_t&>( *this );
    } // method


    /* Vertical maximum with respect to the argument */
    template <typename TP_VEC> Vc_INTRINSIC Vec_m128i_16_int8_t max( const TP_VEC& vec ) const
    {
        if( INSTR >= 41 ) // compile time decision
        { // SSE4.1 only
#if( !defined( __GNUC__ ) || defined( __SSE4_1__ ) )
            return Vec_m128i_16_int8_t( _mm_max_epi8( this->val, TP_SUPER::unwrap( vec ) ) );
#endif
        } // if
        else
        { // Plain SSE2
            __m128i mask = _mm_cmpgt_epi8( this->val, TP_SUPER::unwrap( vec ) );
            return Vec_m128i_16_int8_t(
                _mm_or_si128( _mm_and_si128( mask, this->val ), _mm_andnot_si128( mask, TP_SUPER::unwrap( vec ) ) ) );
        } // else
    } // method


    /* Vertical maximum with respect to the argument */
    template <typename TP_VEC> Vc_INTRINSIC Vec_m128i_16_int8_t min( const TP_VEC& vec ) const
    {
        if( INSTR >= 41 ) // compile time decision
        { // SSE4.1 only
#if( !defined( __GNUC__ ) || defined( __SSE4_1__ ) )
            return Vec_m128i_16_int8_t( _mm_min_epi8( this->val, TP_SUPER::unwrap( vec ) ) );
#endif
        } // if
        else
        { // Plain SSE2
            __m128i mask = _mm_cmplt_epi8( this->val, TP_SUPER::unwrap( vec ) );
            return Vec_m128i_16_int8_t(
                _mm_or_si128( _mm_and_si128( mask, this->val ), _mm_andnot_si128( mask, TP_SUPER::unwrap( vec ) ) ) );
        } // else
    } // method


    /* Parallel addition of four 8 bit values */
    template <typename TP> Vc_INTRINSIC Vec_m128i_16_int8_t operator+( const TP& operand ) const
    {
        return Vec_m128i_16_int8_t( _mm_add_epi8( this->val, TP_SUPER::unwrap( operand ) ) );
    } // method


    /* Parallel substraction of four 8 bit values */
    template <typename TP> Vc_INTRINSIC Vec_m128i_16_int8_t operator-( const TP& operand ) const
    {
        return Vec_m128i_16_int8_t( _mm_sub_epi8( this->val, TP_SUPER::unwrap( operand ) ) );
    } // method


    /* Parallel Comparison (greater than) of all vector elements */
    template <typename TP> Vc_INTRINSIC Vec_m128i_16_int8_t operator>( const TP& operand ) const
    {
        return Vec_m128i_16_int8_t( _mm_cmpgt_epi8( this->val, TP_SUPER::unwrap( operand ) ) );
    } // method


    /* Parallel Comparison (equality) of all vector elements */
    template <typename TP> Vc_INTRINSIC Vec_m128i_16_int8_t operator==( const TP& operand ) const
    {
        return Vec_m128i_16_int8_t( _mm_cmpeq_epi8( this->val, TP_SUPER::unwrap( operand ) ) );
    } // method
}; // struct


/* Vector of 8 signed 16 bit values */
template <int INSTR> struct alignas( 16 ) Vec_m128i_8_int16_t : public Vec_m128i_base<INSTR, Vec_m128i_8_int16_t<INSTR>>
{
    static constexpr size_t SIZE = 8; // 8 times 16 bit
    typedef int16_t TP_ELEMENT;

    typedef Vec_m128i_base<INSTR, Vec_m128i_8_int16_t<INSTR>> TP_SUPER;

    Vec_m128i_8_int16_t( __m128i val ) : TP_SUPER( val )
    {} // constructor


    /* Default constructor does not initialize the vector! */
    Vec_m128i_8_int16_t( void )
    {} // default constructor


    /* Dump the vector std::cout */
    void dump( void )
    {
        std::cout << (int16_t)_mm_extract_epi16( this->val, 0 ) << ", " << (int16_t)_mm_extract_epi16( this->val, 1 )
                  << ", " << (int16_t)_mm_extract_epi16( this->val, 2 ) << ", "
                  << (int16_t)_mm_extract_epi16( this->val, 3 ) << ", " << (int16_t)_mm_extract_epi16( this->val, 4 )
                  << ", " << (int16_t)_mm_extract_epi16( this->val, 5 ) << ", "
                  << (int16_t)_mm_extract_epi16( this->val, 6 ) << ", " << (int16_t)_mm_extract_epi16( this->val, 7 )
                  << ", " << std::endl;
    } // method


    /* Initializes all 8 vector elements to iVal */
    Vc_INTRINSIC Vec_m128i_8_int16_t& set1( const int16_t iVal )
    {
        this->val = _mm_set1_epi16( iVal );
        return static_cast<Vec_m128i_8_int16_t&>( *this );
    } // method


    /* Parallel addition of four 32 bit values */
    template <typename TP> Vc_INTRINSIC Vec_m128i_8_int16_t operator+( const TP& operand ) const
    {
        return Vec_m128i_8_int16_t( _mm_add_epi16( this->val, TP_SUPER::unwrap( operand ) ) );
    } // method


    /* Parallel Comparison (greater than) of four 32 bit integers */
    template <typename TP> Vc_INTRINSIC Vec_m128i_8_int16_t operator>( const TP& operand ) const
    {
        return Vec_m128i_8_int16_t( _mm_cmpgt_epi16( this->val, TP_SUPER::unwrap( operand ) ) );
    } // method


    /* Signed expansion of an array of 4 int8_t values.
     */
    static Vc_INTRINSIC Vec_m128i_8_int16_t expand_int8_t( const int8_t* pArr )
    {
#if( !defined( __GNUC__ ) || defined( __SSSE3__ ) )
        if( INSTR >= 31 ) // compile time decision
        { /* The below code is fast for architectures with a fast _mm_shuffle_epi8
           * implementation; requires SSSE3.
           */
            // Load the 8 int8_t values via loading a single 64 bit float value
            __m128i vZero = _mm_setzero_si128( );
            const __m128i v8Bytes =
                _mm_castpd_si128( _mm_loadl_pd( _mm_castsi128_pd( vZero ), reinterpret_cast<double const*>( pArr ) ) );
            // Shuffle the bytes to the appropriate positions with respect to int16_t
            const __m128i vShuffle =
                _mm_setr_epi8( INT8_C( 0 ), 0x80U, INT8_C( 1 ), 0x80U, INT8_C( 2 ), 0x80U, INT8_C( 3 ), 0x80U,
                               INT8_C( 4 ), 0x80U, INT8_C( 5 ), 0x80U, INT8_C( 6 ), 0x80U, INT8_C( 7 ), 0x80U );
            const __m128i vShuff8Bytes = _mm_shuffle_epi8( v8Bytes, vShuffle );
            // if greater 127, it was a negative value (255 = -1, ... , 128 = - 128)
            const __m128i vCompValues = _mm_set1_epi16( INT8_MAX );
            // v_mask = 1 for the values greater 127, all other 0
            const __m128i vMask = _mm_cmpgt_epi16( vShuff8Bytes, vCompValues );
            // expand negative values to 32 bit size
            const __m128i vNegValues = _mm_or_si128( vShuff8Bytes, _mm_set1_epi16( (int16_t)0xFF00 ) );
            // blend positive and negative values using the mask
            return Vec_m128i_8_int16_t( vShuff8Bytes ).blendv( vNegValues, vMask );
        } // if
        else
#endif
        { /* Plain SSE2 */
            return Vec_m128i_8_int16_t( _mm_setr_epi16( pArr[ 0 ], pArr[ 1 ], pArr[ 2 ], pArr[ 3 ], pArr[ 4 ],
                                                        pArr[ 5 ], pArr[ 6 ], pArr[ 7 ] ) );
        } // else
    } // method


    /* Horizontal maximum of the eight 16 bit values.
     */
    Vc_INTRINSIC int16_t hmax( void )
    {
        __m128i max1 = _mm_shuffle_epi32( this->val, _MM_SHUFFLE( 0, 0, 3, 2 ) );
        __m128i max2 = _mm_max_epi16( this->val, max1 );
        __m128i max3 = _mm_shuffle_epi32( max2, _MM_SHUFFLE( 0, 0, 0, 1 ) );
        __m128i max4 = _mm_max_epi16( max2, max3 );
        __m128i max5 = _mm_shufflelo_epi16( max4, _MM_SHUFFLE( 0, 0, 0, 1 ) );
        __m128i max6 = _mm_max_epi16( max4, max5 );
        return (int16_t)_mm_extract_epi16( max6, 0 );
    } // method
}; // struct


/* Vector of 4 signed 32 bit values */
template <int INSTR> struct alignas( 16 ) Vec_m128i_4_int32_t : public Vec_m128i_base<INSTR, Vec_m128i_4_int32_t<INSTR>>
{
    static constexpr size_t SIZE = 4; // 4 times 32 bit (int32_t)
    typedef int32_t TP_ELEMENT;

    typedef Vec_m128i_base<INSTR, Vec_m128i_4_int32_t<INSTR>> TP_SUPER;

    Vec_m128i_4_int32_t( __m128i val ) : TP_SUPER( val )
    {} // constructor


    /* Default constructor does not initialize the vector! */
    Vec_m128i_4_int32_t( void )
    {} // default constructor


    /* For debuging ... */
    void dump( void )
    {
        std::cout << _mm_extract_epi32( this->val, 0 ) << ", " << _mm_extract_epi32( this->val, 1 ) << ", "
                  << _mm_extract_epi32( this->val, 2 ) << ", " << _mm_extract_epi32( this->val, 3 ) << std::endl;
    } // method


    /* Getter, delivers value */
    Vc_INTRINSIC TP_ELEMENT operator[]( const size_t iIdx )
    {
        return _mm_extract_epi32( this->val, iIdx );
    } // method


    /* Initializes all 4 vector elements to iVal */
    Vc_INTRINSIC Vec_m128i_4_int32_t& set1( const int32_t iVal )
    {
        this->val = _mm_set1_epi32( iVal );
        return static_cast<Vec_m128i_4_int32_t&>( *this );
    } // method


    /* Initializes the 4 vector elements with the 4 arguments */
    Vc_INTRINSIC Vec_m128i_4_int32_t& set4( const int32_t iVal0, const int32_t iVal1, const int32_t iVal2,
                                            const int32_t iVal3 )
    {
        this->val = _mm_setr_epi32( iVal0, iVal1, iVal2, iVal3 );
        return static_cast<Vec_m128i_4_int32_t&>( *this );
    } // method


    /* Signed expansion of an array of 4 int8_t values.
     */
    static Vc_INTRINSIC Vec_m128i_4_int32_t expand_int8_t( const int8_t* pArr )
    {
#if( !defined( __GNUC__ ) || defined( __SSSE3__ ) )
        if( INSTR >= 31 ) // compile time decision
        { /* The below code is fast for architectures with a fast _mm_shuffle_epi8
           * implementation; requires SSSE3.
           * Load the 4 int8_t values via loading a single 32 bit integer */
            const __m128i v4Bytes = _mm_cvtsi32_si128( *( reinterpret_cast<const int32_t*>( pArr ) ) );
            /* Shuffle the bytes to the correct positions */
            const __m128i vShuffle = _mm_setr_epi8( 0, 0x80U, 0x80U, 0x80U, 1, 0x80U, 0x80U, 0x80U, 2, 0x80U, 0x80U,
                                                    0x80U, 3, 0x80U, 0x80U, 0x80U );
            const __m128i vShuff4Bytes = _mm_shuffle_epi8( v4Bytes, vShuffle );
            /* if greater 127, it was a negative value (255 = -1, ... , 128 = - 128) */
            const __m128i vCompValues = _mm_setr_epi32( 0x7F, 0x7F, 0x7F, 0x7F );
            /* v_mask = 1 for the values greater 127, all other 0 */
            const __m128i vMask = _mm_cmpgt_epi32( vShuff4Bytes, vCompValues );
            /* expand negative values to 32 bit size */
            const __m128i vNegValues =
                _mm_or_si128( vShuff4Bytes, _mm_setr_epi32( 0xFFFFFF00, 0xFFFFFF00, 0xFFFFFF00, 0xFFFFFF00 ) );

            /* blend positive and negative values using the mask */
            return Vec_m128i_4_int32_t( vShuff4Bytes ).blendv( vNegValues, vMask );
        } // if
        else
#endif
        { /* Plain SSE2 */
            return Vec_m128i_4_int32_t( _mm_setr_epi32( pArr[ 0 ], pArr[ 1 ], pArr[ 2 ], pArr[ 3 ] ) );
        } // else
    } // method


    /* Parallel addition of four 32 bit values */
    template <typename TP> Vc_INTRINSIC Vec_m128i_4_int32_t operator+( const TP& operand ) const
    {
        return Vec_m128i_4_int32_t( _mm_add_epi32( this->val, TP_SUPER::unwrap( operand ) ) );
    } // method


    /* Parallel Comparison (greater than) of four 32 bit integers */
    template <typename TP> Vc_INTRINSIC Vec_m128i_4_int32_t operator>( const TP& operand ) const
    {
        return Vec_m128i_4_int32_t( _mm_cmpgt_epi32( this->val, TP_SUPER::unwrap( operand ) ) );
    } // method


    /* Parallel max of four 32 bit integers */
    template <typename TP> Vc_INTRINSIC Vec_m128i_4_int32_t max( const TP& vec ) const
    {
        if( INSTR >= 41 ) // compile time decision
        { // SSE4.1 only
#if( !defined( __GNUC__ ) || defined( __SSE4_1__ ) )
            return Vec_m128i_4_int32_t( _mm_max_epi32( this->val, TP_SUPER::unwrap( vec ) ) );
#endif
        } // if
        else
        { // Plain SSE2
            __m128i mask = _mm_cmpgt_epi32( this->val, TP_SUPER::unwrap( vec ) );
            // ssp_logical_bitwise_select_SSE2
            return Vec_m128i_4_int32_t(
                _mm_or_si128( _mm_and_si128( mask, this->val ), _mm_andnot_si128( mask, TP_SUPER::unwrap( vec ) ) ) );
        } // else
    } // method


    /* Horizontal maximum of the four 32 bit values.
     * Found at:
     * https://stackoverflow.com/questions/9877700/getting-max-value-in-a-m128i-vector-with-sse
     */
    Vc_INTRINSIC int32_t hmax( void )
    {
        auto max1 = _mm_shuffle_epi32( this->val, _MM_SHUFFLE( 0, 0, 3, 2 ) );
        auto max2 = this->max( max1 ); // _mm_max_epi32( this->val, max1 );
        auto max3 = _mm_shuffle_epi32( max2.val, _MM_SHUFFLE( 0, 0, 0, 1 ) );
        auto max4 = max2.max( max3 ); //  _mm_max_epi32( max2, max3 );
        return _mm_cvtsi128_si32( max4.val );
    } // method
}; // struct


/* Check sizes of all vector types */
// static_assert( sizeof( Vec_m128i_16_int8_t<0> ) == sizeof( __m128i ) );
// static_assert( sizeof( Vec_m128i_8_int16_t<0> ) == sizeof( __m128i ) );
// static_assert( sizeof( Vec_m128i_4_int32_t<0> ) == sizeof( __m128i ) );
