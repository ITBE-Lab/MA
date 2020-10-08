#pragma once
#include <cassert>
#include <cstring>
#include <iostream>
#include <memory>
#include <unordered_map>

#if defined( __GNUC__ )
#include <stdlib.h>
#endif

/* Found at https://stackoverflow.com/questions/8456236/how-is-a-vectors-data-aligned
 * std::vector<T, AlignmentAllocator<T, 16> > bla;
 */
template <typename T, std::size_t N = 16> class AlignmentAllocator
{
  public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef T* pointer;
    typedef const T* const_pointer;

    typedef T& reference;
    typedef const T& const_reference;

  public:
    inline AlignmentAllocator( ) throw( )
    {}

    template <typename T2> inline AlignmentAllocator( const AlignmentAllocator<T2, N>& ) throw( )
    {}

    inline ~AlignmentAllocator( ) throw( )
    {}

    inline pointer adress( reference r )
    {
        return &r;
    }

    inline const_pointer adress( const_reference r ) const
    {
        return &r;
    }

    inline pointer allocate( size_type n )
    {
#if defined( _MSC_VER )
        return (pointer)_aligned_malloc( n * sizeof( value_type ), N );
#elif defined( __GNUC__ )
        std::cout << "ALLOCATE MEMORY" << std::endl;
        return (pointer)aligned_alloc( N, n * sizeof( value_type ) );
#else
#error "Aligned memory allocation unknown for this compiler"
#endif
    }

    inline void deallocate( pointer p, size_type )
    {
#if defined( _MSC_VER )
        _aligned_free( p );
#elif defined( __GNUC__ )
        free( p );
#else
#error "Aligned memory allocation unknown for this compiler"
#endif
    }

    inline void construct( pointer p, const value_type& wert )
    {
        new( p ) value_type( wert );
    }

    inline void destroy( pointer p )
    {
        p->~value_type( );
    }

    inline size_type max_size( ) const throw( )
    {
        return size_type( -1 ) / sizeof( value_type );
    }

    template <typename T2> struct rebind
    {
        typedef AlignmentAllocator<T2, N> other;
    };

    bool operator!=( const AlignmentAllocator<T, N>& other ) const
    {
        return !( *this == other );
    }

    // Returns true if and only if storage allocated from *this
    // can be deallocated from other, and vice versa.
    // Always returns true for stateless allocators.
    bool operator==( const AlignmentAllocator<T, N>& other ) const
    {
        return true;
    }
};

#ifdef USE_VECTOR
template <typename SIMD_VEC_TP>
#endif
class AlignedMemoryManager
{
  public:
    uint8_t* pMemMat = NULL;
    uint8_t* pMemMatAligned = NULL;
    size_t uiCapacityMemMat = 0;

    void* pMemH = NULL;
    size_t uiCapacityMemH = 0;

    AlignedMemoryManager( )
    {
        //// std::cout << "Make Memory Manager" << std::endl;
    } // default constructor

#ifdef USE_VECTOR
    std::vector<SIMD_VEC_TP, AlignmentAllocator<SIMD_VEC_TP, sizeof( SIMD_VEC_TP )>> vs_u;
    std::vector<SIMD_VEC_TP, AlignmentAllocator<SIMD_VEC_TP, sizeof( SIMD_VEC_TP )>> vs_v;
    std::vector<SIMD_VEC_TP, AlignmentAllocator<SIMD_VEC_TP, sizeof( SIMD_VEC_TP )>> vs_x;
    std::vector<SIMD_VEC_TP, AlignmentAllocator<SIMD_VEC_TP, sizeof( SIMD_VEC_TP )>> vs_y;
    std::vector<SIMD_VEC_TP, AlignmentAllocator<SIMD_VEC_TP, sizeof( SIMD_VEC_TP )>> vs_x2;
    std::vector<SIMD_VEC_TP, AlignmentAllocator<SIMD_VEC_TP, sizeof( SIMD_VEC_TP )>> vs_y2;
    std::vector<SIMD_VEC_TP, AlignmentAllocator<SIMD_VEC_TP, sizeof( SIMD_VEC_TP )>> vs_s;
#endif

    /* Delivers and clears aligned memory.
     * Alignment is with respect to the type __mXXXi.
     * Allocates on demand.
     */
    template <typename __mXXXi>
#ifdef USE_VECTOR
    __mXXXi* reserveMemMatrix( size_t uiRequestedSize, size_t uiVecSize )
#else
    __mXXXi* reserveMemMatrix( size_t uiRequestedSize )
#endif
    {
        if( ( uiRequestedSize > this->uiCapacityMemMat ) || ( this->pMemMat == NULL ) )
        {
            //// std::cout << "Allocate memory " << uiRequestedSize << " " << this->uiCapacityMemMat <<std::endl;
            if( this->pMemMat != NULL )
                free( this->pMemMat );

            this->pMemMat = (uint8_t*)malloc( uiRequestedSize * sizeof( __mXXXi ) + sizeof( __mXXXi ) );
            this->pMemMatAligned = (uint8_t*)( ( ( size_t )( this->pMemMat ) + ( sizeof( __mXXXi ) - 1 ) ) /
                                               sizeof( __mXXXi ) * sizeof( __mXXXi ) );

            // this->pMemMat = (uint8_t*)aligned_alloc( sizeof(__mXXXi), uiRequestedSize * sizeof( __mXXXi ) );
            // this->pMemMatAligned = this->pMemMat;
            this->uiCapacityMemMat = uiRequestedSize;
        } // if
        memset( this->pMemMatAligned, 0, uiRequestedSize * sizeof( __mXXXi ) );

#ifdef USE_VECTOR
        /* Set all vectors to appropriate size */
        vs_u.resize( uiVecSize );
        vs_v.resize( uiVecSize );
        vs_x.resize( uiVecSize );
        vs_y.resize( uiVecSize );
        vs_x2.resize( uiVecSize );
        vs_y2.resize( uiVecSize );
        vs_s.resize( uiVecSize );

        /* Set all vectors to zero */
        for( auto uiIndex = 0; uiIndex < uiVecSize; uiIndex++ )
        {
            vs_u[ uiIndex ].setzero( );
            vs_v[ uiIndex ].setzero( );
            vs_x[ uiIndex ].setzero( );
            vs_y[ uiIndex ].setzero( );
            vs_x2[ uiIndex ].setzero( );
            vs_y2[ uiIndex ].setzero( );
            vs_s[ uiIndex ].setzero( );
        } // for
#endif

        return (__mXXXi*)( this->pMemMatAligned );
    } // method

    template <typename T_Scoring> T_Scoring* reserveMemH_Vector( size_t uiRequestedSize )
    {
        if( ( uiRequestedSize > this->uiCapacityMemH ) || ( this->pMemH == NULL ) )
        {
            if( this->pMemH != NULL )
                free( this->pMemH );

            this->pMemH = (void*)malloc( uiRequestedSize * sizeof( T_Scoring ) );
            this->uiCapacityMemH = uiRequestedSize;
        } // if

        return (T_Scoring*)( this->pMemH );
    } // method

    ~AlignedMemoryManager( )
    {
        if( this->pMemMat != NULL )
            free( this->pMemMat );

        if( this->pMemH != NULL )
            free( this->pMemH );
    } // destructor
}; // class

template <typename TP_DIFF_VEC, int64_t CHUNK_SIZE, size_t HASH_TABLE_GB_MIN_SIZE> class CIGARMemoryManager
{
    std::unordered_map<int64_t, std::array<TP_DIFF_VEC, CHUNK_SIZE>> xData;
    void* mem2 = nullptr;
    TP_DIFF_VEC* p_pstruct = nullptr;
    int64_t iPROffset = 0;
    TP_DIFF_VEC xZero;
    std::array<TP_DIFF_VEC, CHUNK_SIZE> xZeroArray;

  public:
    CIGARMemoryManager( )
    {} // constructor

    void resize( size_t uiSize )
    {
        if( uiSize * sizeof( TP_DIFF_VEC ) <= HASH_TABLE_GB_MIN_SIZE * 1073741824 )
        {
            if( p_pstruct != nullptr ) // safety check
                free( mem2 );
            // std::cout << "allocating " << uiSize * sizeof( TP_DIFF_VEC ) << " / " << HASH_TABLE_GB_MIN_SIZE *
            // 1073741824
            //          << " bytes." << std::endl;
            // std::cout << "that is " << uiSize * sizeof( TP_DIFF_VEC ) / 1073741824.0 << " / " <<
            // HASH_TABLE_GB_MIN_SIZE
            //          << " gigabytes." << std::endl;
            mem2 = malloc( uiSize * sizeof( TP_DIFF_VEC ) );
            p_pstruct = (TP_DIFF_VEC*)( ( ( (size_t)mem2 + ( sizeof( TP_DIFF_VEC ) - 1 ) ) / sizeof( TP_DIFF_VEC ) ) *
                                        sizeof( TP_DIFF_VEC ) );
        } // if
        else
        {
            // std::cout << "using hash table because: " << uiSize * sizeof( TP_DIFF_VEC ) << " / "
            //          << HASH_TABLE_GB_MIN_SIZE * 1073741824 << " bytes are required." << std::endl;
            // std::cout << "that is " << uiSize * sizeof( TP_DIFF_VEC ) / 1073741824.0 << " / " <<
            // HASH_TABLE_GB_MIN_SIZE
            //          << " gigabytes." << std::endl;

            // zero initialize zero value
            xZero.setzero( );
            // zero initialize zero array
            for( size_t uiI = 0; uiI < xZeroArray.size( ); uiI++ )
                xZeroArray[ uiI ] = xZero;
        } // else
    } // method

    void set( size_t iI, TP_DIFF_VEC xVal )
    {
        if( p_pstruct == nullptr )
        {
            for( size_t uiX = 0; uiX < sizeof( TP_DIFF_VEC ); uiX++ )
                if( xVal.at( uiX ) != 0 )
                {
                    // insert zero array (this triggers a copy; however this is necessary) if no array exists for that
                    // chunk.
                    // Then sets the appropriate element in the array.
                    // We need to zero initialize arrays manually as their constructors do not guarantee that
                    // and because we ignore inserting zero elements.
                    xData.try_emplace( ( iPROffset + iI ) / CHUNK_SIZE, xZeroArray )
                        .first->second[ ( iPROffset + iI ) % CHUNK_SIZE ] = xVal;
                    return;
                } // if
        } // if
        else
            p_pstruct[ iPROffset + iI ] = xVal;
    }

    TP_DIFF_VEC& get( size_t iI )
    {
        if( p_pstruct == nullptr )
        {
            auto xIt = xData.find( iI / CHUNK_SIZE );
            if( xIt == xData.end( ) )
                return xZero;
            return xIt->second[ iI % CHUNK_SIZE ];
        } // if
        else
            return p_pstruct[ iI ];
    }

    void setPROffset( int64_t iI )
    {
        iPROffset = iI;
    } // method

    template <typename TP_RET_TYPE> TP_RET_TYPE accessAs( size_t uiIndex )
    {
        if( p_pstruct == nullptr )
        {
            assert( sizeof( TP_RET_TYPE ) <= sizeof( TP_DIFF_VEC ) );
            size_t iI = uiIndex / ( sizeof( TP_DIFF_VEC ) / sizeof( TP_RET_TYPE ) );
            size_t uiModulo = uiIndex % ( sizeof( TP_DIFF_VEC ) / sizeof( TP_RET_TYPE ) );
            return ( (TP_RET_TYPE*)&this->get( iI ) )[ uiModulo ];
        } // if
        else
            return ( (TP_RET_TYPE*)p_pstruct )[ uiIndex ];
    }

    ~CIGARMemoryManager( )
    {
        if( p_pstruct != nullptr )
            free( mem2 );
    }
}; // class

//// void ksw_extd2_old( void *km,
//// 					int qlen,
//// 					const uint8_t *query,
//// 					int tlen, // reference length
//// 					const uint8_t *target, // reference
//// 					int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2,
//// 					int w, // bandwidth
//// 					int zdrop, //
//// 					int end_bonus,
//// 					int flag,
//// 					ksw_extz_t *ez,
//// 					AlignedMemoryManager &rxMemManager
//// 				  );
