/* This is header delivers a class for SIMD vector simulation.
 * The idea is to rely on this class if there is no support of SIMD vectors on a specific platform.
 * The code is not optimized in its current form.
 * TODO: Compatibility with float types.
 */

#pragma once
#include <stdint.h>

#define USE_C_MEM

#ifdef USE_C_MEM
#include <cstring> // for std::memcpy
#endif

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


/* WARNING: TP_VEC_ELE must be an integer-type!
 */
template <size_t VEC_SIZE, typename TP_VEC_ELE> struct alignas( 16 ) Vec_scalar
{
    static constexpr size_t SIZE = VEC_SIZE; // Number of elements of the vector

    typedef Vec_scalar<VEC_SIZE, TP_VEC_ELE> TP_SELF;

    TP_VEC_ELE val[ VEC_SIZE ]; // val

    static constexpr size_t BYTE_SIZE = VEC_SIZE * sizeof( TP_VEC_ELE );


    /* Default creates new array */
    Vec_scalar( void ) // :
                       // val( new int[SIZE] )
    {} // constructor


    /* Initializing constructor */
    Vec_scalar( const TP_VEC_ELE* pInitValues )
    {
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            this->val[ uiItr ] = pInitValues[ uiItr ];
    } // constructor

    Vec_scalar( int8_t iVal )
    {
        this->set1( iVal );
    } // constructor


    Vc_INTRINSIC int8_t* pointerVal( void )
    {
        return reinterpret_cast<int8_t*>( this->val );
    } // method


    Vc_INTRINSIC int8_t byteGet( const size_t iIdx )
    {
        return reinterpret_cast<int8_t*>( this->val )[ iIdx ];
    } // method


    Vc_INTRINSIC void byteSet( const size_t iIdx, const int8_t iVal )
    {
        reinterpret_cast<int8_t*>( this->val )[ iIdx ] = iVal;
    } // method


    /* Setter, delivers reference */
    Vc_INTRINSIC TP_VEC_ELE& operator[]( const size_t iIdx )
    {
        return this->val[ iIdx ];
    } // method


    /* Getter, delivers value */
    Vc_INTRINSIC TP_VEC_ELE operator[]( const size_t iIdx ) const
    {
        return this->val[ iIdx ];
    } // method


    Vc_INTRINSIC TP_SELF& setzero( void )
    { // we could use memset instead
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            this->val[ uiItr ] = 0;
        return *this;
    } // method


    /* Set the first element in the vector to iVal.
     * TO DO: arg to int32_t
     */
    Vc_INTRINSIC TP_SELF& moveAndZero32( int8_t iVal )
    {
        for( size_t uiItr = 1; uiItr < SIZE; uiItr++ )
            this->val[ uiItr ] = 0;

        this->val[ 0 ] = iVal;
        return *this;
    } // method


    /* Copy assignment */
    Vc_INTRINSIC TP_SELF& operator=( const TP_SELF& rhs )
    {
        // we could use memcpy instead
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            this->val[ uiItr ] = rhs.val[ uiItr ];
        return *this;
    } // method


    /* Store the vector on the given address */
    Vc_INTRINSIC const TP_SELF& store( void* addr ) const
    {
#ifdef USE_C_MEM
        std::memcpy( addr, this->val, BYTE_SIZE );
#else
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            ( reinterpret_cast<TP_VEC_ELE*>( addr ) )[ uiItr ] = this->val[ uiItr ];
#endif
        return *this;
    } // method


    /* Store the vector on the given address */
    Vc_INTRINSIC TP_SELF& load( void* addr )
    {
        // we could use memcpy instead
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            this->val[ uiItr ] = ( reinterpret_cast<TP_VEC_ELE*>( addr ) )[ uiItr ];
        return *this;
    } // method


    /* Number of bits in the vector vector.
     * Here we assume a 64-bit architecture.
     * See: https://github.com/WojciechMula/sse-popcount/blob/master/popcnt-bit-parallel-scalar.cpp
     */
    Vc_INTRINSIC size_t popcnt( void ) //// BYTE_SIZE
    {
        static const uint8_t aLookup8bit[ 256 ] = {
            /* 0 */ 0,  /* 1 */ 1,  /* 2 */ 1,  /* 3 */ 2,
            /* 4 */ 1,  /* 5 */ 2,  /* 6 */ 2,  /* 7 */ 3,
            /* 8 */ 1,  /* 9 */ 2,  /* a */ 2,  /* b */ 3,
            /* c */ 2,  /* d */ 3,  /* e */ 3,  /* f */ 4,
            /* 10 */ 1, /* 11 */ 2, /* 12 */ 2, /* 13 */ 3,
            /* 14 */ 2, /* 15 */ 3, /* 16 */ 3, /* 17 */ 4,
            /* 18 */ 2, /* 19 */ 3, /* 1a */ 3, /* 1b */ 4,
            /* 1c */ 3, /* 1d */ 4, /* 1e */ 4, /* 1f */ 5,
            /* 20 */ 1, /* 21 */ 2, /* 22 */ 2, /* 23 */ 3,
            /* 24 */ 2, /* 25 */ 3, /* 26 */ 3, /* 27 */ 4,
            /* 28 */ 2, /* 29 */ 3, /* 2a */ 3, /* 2b */ 4,
            /* 2c */ 3, /* 2d */ 4, /* 2e */ 4, /* 2f */ 5,
            /* 30 */ 2, /* 31 */ 3, /* 32 */ 3, /* 33 */ 4,
            /* 34 */ 3, /* 35 */ 4, /* 36 */ 4, /* 37 */ 5,
            /* 38 */ 3, /* 39 */ 4, /* 3a */ 4, /* 3b */ 5,
            /* 3c */ 4, /* 3d */ 5, /* 3e */ 5, /* 3f */ 6,
            /* 40 */ 1, /* 41 */ 2, /* 42 */ 2, /* 43 */ 3,
            /* 44 */ 2, /* 45 */ 3, /* 46 */ 3, /* 47 */ 4,
            /* 48 */ 2, /* 49 */ 3, /* 4a */ 3, /* 4b */ 4,
            /* 4c */ 3, /* 4d */ 4, /* 4e */ 4, /* 4f */ 5,
            /* 50 */ 2, /* 51 */ 3, /* 52 */ 3, /* 53 */ 4,
            /* 54 */ 3, /* 55 */ 4, /* 56 */ 4, /* 57 */ 5,
            /* 58 */ 3, /* 59 */ 4, /* 5a */ 4, /* 5b */ 5,
            /* 5c */ 4, /* 5d */ 5, /* 5e */ 5, /* 5f */ 6,
            /* 60 */ 2, /* 61 */ 3, /* 62 */ 3, /* 63 */ 4,
            /* 64 */ 3, /* 65 */ 4, /* 66 */ 4, /* 67 */ 5,
            /* 68 */ 3, /* 69 */ 4, /* 6a */ 4, /* 6b */ 5,
            /* 6c */ 4, /* 6d */ 5, /* 6e */ 5, /* 6f */ 6,
            /* 70 */ 3, /* 71 */ 4, /* 72 */ 4, /* 73 */ 5,
            /* 74 */ 4, /* 75 */ 5, /* 76 */ 5, /* 77 */ 6,
            /* 78 */ 4, /* 79 */ 5, /* 7a */ 5, /* 7b */ 6,
            /* 7c */ 5, /* 7d */ 6, /* 7e */ 6, /* 7f */ 7,
            /* 80 */ 1, /* 81 */ 2, /* 82 */ 2, /* 83 */ 3,
            /* 84 */ 2, /* 85 */ 3, /* 86 */ 3, /* 87 */ 4,
            /* 88 */ 2, /* 89 */ 3, /* 8a */ 3, /* 8b */ 4,
            /* 8c */ 3, /* 8d */ 4, /* 8e */ 4, /* 8f */ 5,
            /* 90 */ 2, /* 91 */ 3, /* 92 */ 3, /* 93 */ 4,
            /* 94 */ 3, /* 95 */ 4, /* 96 */ 4, /* 97 */ 5,
            /* 98 */ 3, /* 99 */ 4, /* 9a */ 4, /* 9b */ 5,
            /* 9c */ 4, /* 9d */ 5, /* 9e */ 5, /* 9f */ 6,
            /* a0 */ 2, /* a1 */ 3, /* a2 */ 3, /* a3 */ 4,
            /* a4 */ 3, /* a5 */ 4, /* a6 */ 4, /* a7 */ 5,
            /* a8 */ 3, /* a9 */ 4, /* aa */ 4, /* ab */ 5,
            /* ac */ 4, /* ad */ 5, /* ae */ 5, /* af */ 6,
            /* b0 */ 3, /* b1 */ 4, /* b2 */ 4, /* b3 */ 5,
            /* b4 */ 4, /* b5 */ 5, /* b6 */ 5, /* b7 */ 6,
            /* b8 */ 4, /* b9 */ 5, /* ba */ 5, /* bb */ 6,
            /* bc */ 5, /* bd */ 6, /* be */ 6, /* bf */ 7,
            /* c0 */ 2, /* c1 */ 3, /* c2 */ 3, /* c3 */ 4,
            /* c4 */ 3, /* c5 */ 4, /* c6 */ 4, /* c7 */ 5,
            /* c8 */ 3, /* c9 */ 4, /* ca */ 4, /* cb */ 5,
            /* cc */ 4, /* cd */ 5, /* ce */ 5, /* cf */ 6,
            /* d0 */ 3, /* d1 */ 4, /* d2 */ 4, /* d3 */ 5,
            /* d4 */ 4, /* d5 */ 5, /* d6 */ 5, /* d7 */ 6,
            /* d8 */ 4, /* d9 */ 5, /* da */ 5, /* db */ 6,
            /* dc */ 5, /* dd */ 6, /* de */ 6, /* df */ 7,
            /* e0 */ 3, /* e1 */ 4, /* e2 */ 4, /* e3 */ 5,
            /* e4 */ 4, /* e5 */ 5, /* e6 */ 5, /* e7 */ 6,
            /* e8 */ 4, /* e9 */ 5, /* ea */ 5, /* eb */ 6,
            /* ec */ 5, /* ed */ 6, /* ee */ 6, /* ef */ 7,
            /* f0 */ 4, /* f1 */ 5, /* f2 */ 5, /* f3 */ 6,
            /* f4 */ 5, /* f5 */ 6, /* f6 */ 6, /* f7 */ 7,
            /* f8 */ 5, /* f9 */ 6, /* fa */ 6, /* fb */ 7,
            /* fc */ 6, /* fd */ 7, /* fe */ 7, /* ff */ 8}; // lookup table

        uint8_t* pData = reinterpret_cast<uint8_t*>( this->val );
        size_t uiResult = 0;
        size_t uiItr = 0;

        while( uiItr + 4 * 8 <= BYTE_SIZE )
        {
            uint64_t uiPartial = 0; // packed_byte
#define ITER                                                                                            \
    {                                                                                                   \
        const uint64_t t1 = *reinterpret_cast<const uint64_t*>( pData + uiItr );                        \
        const uint64_t t2 = ( t1 & 0x5555555555555555llu ) + ( ( t1 >> 1 ) & 0x5555555555555555llu );   \
        const uint64_t t3 = ( t2 & 0x3333333333333333llu ) + ( ( t2 >> 2 ) & 0x3333333333333333llu );   \
        const uint64_t t4 = ( t3 & 0x0f0f0f0f0f0f0f0fllu ) + ( ( t3 >> 4 ) & 0x0f0f0f0f0f0f0f0fllu );   \
        uiPartial += t4;                                                                                \
        uiItr += 8;                                                                                     \
    }
            ITER ITER ITER ITER
#undef ITER
                const uint64_t t5 = ( uiPartial & 0x00ff00ff00ff00ffllu ) +
                                    ( ( uiPartial >> 8 ) & 0x00ff00ff00ff00ffllu );
            const uint64_t t6 =
                ( t5 & 0x0000ffff0000ffffllu ) + ( ( t5 >> 16 ) & 0x0000ffff0000ffffllu );
            const uint64_t t7 =
                ( t6 & 0x00000000ffffffffllu ) + ( ( t6 >> 32 ) & 0x00000000ffffffffllu );
            uiResult += t7;
        } // while

        for( ; uiItr < BYTE_SIZE; uiItr++ )
            uiResult += aLookup8bit[ pData[ uiItr ] ];

        return uiResult;
    } // method


    /* Logical it is a right shifting, but low level its a left shifting */
    template <unsigned int SHIFT> Vc_INTRINSIC TP_SELF shiftr( void ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < BYTE_SIZE - SHIFT; uiItr++ )
            reinterpret_cast<char*>( aTmp.val )[ uiItr ] =
                reinterpret_cast<const char*>( this->val )[ uiItr + SHIFT ];
        for( size_t uiItr = BYTE_SIZE - SHIFT; uiItr < BYTE_SIZE; uiItr++ )
            reinterpret_cast<char*>( aTmp.val )[ uiItr ] = 0;

        return aTmp;
    } // method


    /* Logical it is a left shifting, but low level its a right shifting */
    template <unsigned int SHIFT> Vc_INTRINSIC TP_SELF shiftl( void ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < BYTE_SIZE - SHIFT; uiItr++ )
            reinterpret_cast<char*>( aTmp.val )[ uiItr + SHIFT ] =
                reinterpret_cast<const char*>( this->val )[ uiItr ];
        for( size_t uiItr = 0; uiItr < SHIFT; uiItr++ )
            reinterpret_cast<char*>( aTmp.val )[ uiItr ] = 0;

        return aTmp;
    } // method


    /* Vertical blend using a mask (simulates _mm_blendv_epi8) */
    template <typename TP_VEC, typename TP_MASK>
    Vc_INTRINSIC TP_SELF blendv( const TP_VEC& vec, const TP_MASK& mask ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < BYTE_SIZE; uiItr++ )
            reinterpret_cast<char*>( aTmp.val )[ uiItr ] =
                ( reinterpret_cast<const char*>( mask.val )[ uiItr ] & 0x80 )
                    ? reinterpret_cast<const char*>( vec.val )[ uiItr ]
                    : reinterpret_cast<const char*>( this->val )[ uiItr ];

        return aTmp;
    } // method


    /* Bitwise OR over the complete 128 bit vector */
    template <typename TP_VEC> Vc_INTRINSIC TP_SELF operator|( const TP_VEC& operand ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            aTmp.val[ uiItr ] = this->val[ uiItr ] | operand.val[ uiItr ];

        return aTmp;
    } // method


    /* Bitwise AND over the complete 128 bit vector */
    template <typename TP_VEC> Vc_INTRINSIC TP_SELF operator&( const TP_VEC& operand ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            aTmp.val[ uiItr ] = this->val[ uiItr ] & operand.val[ uiItr ];

        return aTmp;
    } // method

	/* Bitwise AND over the complete 128 bit vector */
    template <typename TP_VEC> Vc_INTRINSIC TP_SELF andnot( const TP_VEC& operand ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            aTmp.val[ uiItr ] = !(this->val[ uiItr ] & operand.val[ uiItr ]);

        return aTmp;
    } // method


    /* Initialize all vector elements with iVal */
    Vc_INTRINSIC TP_SELF& set1( TP_VEC_ELE iVal )
    {
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            this->val[ uiItr ] = iVal;
        return *this;
    } // method


    /* Vertical maximum with respect to the argument */
    Vc_INTRINSIC TP_SELF max( const TP_SELF& vec ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            aTmp.val[ uiItr ] =
                this->val[ uiItr ] > vec.val[ uiItr ] ? this->val[ uiItr ] : vec.val[ uiItr ];
        return aTmp;
    }; // method


    /* Vertical minimum with respect to the argument */
    Vc_INTRINSIC TP_SELF min( const TP_SELF& vec ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            aTmp.val[ uiItr ] =
                this->val[ uiItr ] < vec.val[ uiItr ] ? this->val[ uiItr ] : vec.val[ uiItr ];
        return aTmp;
    }; // method


    Vc_INTRINSIC TP_SELF operator+( const TP_SELF& operand ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            aTmp.val[ uiItr ] = this->val[ uiItr ] + operand.val[ uiItr ];
        return aTmp;
    } // method


    Vc_INTRINSIC TP_SELF operator-( const TP_SELF& operand ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            aTmp.val[ uiItr ] = this->val[ uiItr ] - operand.val[ uiItr ];
        return aTmp;
    } // method


    /* Parallel Comparison (equality) of all vector elements.
     * WARNING: Code is not compatible with floating point types.
     */
    Vc_INTRINSIC TP_SELF operator==( const TP_SELF& operand ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            aTmp.val[ uiItr ] =
                this->val[ uiItr ] == operand.val[ uiItr ] ? static_cast<TP_VEC_ELE>( -1 ) : 0;
        return aTmp;
    } // method


    /* Parallel Comparison (greater than) of all vector elements
     * WARNING: Code is not compatible with floating point types.
     */
    Vc_INTRINSIC TP_SELF operator>( const TP_SELF& operand ) const
    {
        TP_SELF aTmp;
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            aTmp.val[ uiItr ] =
                this->val[ uiItr ] > operand.val[ uiItr ] ? static_cast<TP_VEC_ELE>( -1 ) : 0;
        return aTmp;
    } // method


    /* Signed expansion of an array of 4 int8_t values.
     */
    Vc_INTRINSIC TP_SELF& expand_int8_t( const int8_t* pArr )
    {
        for( size_t uiItr = 0; uiItr < SIZE; uiItr++ )
            this->val[ uiItr ] = static_cast<TP_VEC_ELE>( pArr[ uiItr ] );
        return *this;
    } // method


    /* Horizontal maximum of the eight 16 bit values.
     */
    Vc_INTRINSIC TP_VEC_ELE hmax( void )
    {
        TP_VEC_ELE max = this->val[ 0 ];
        for( size_t uiItr = 1; uiItr < SIZE; uiItr++ )
            if( this->val[ uiItr ] > max )
                max = this->val[ uiItr ];
        return max;
    } // method
}; // struct
