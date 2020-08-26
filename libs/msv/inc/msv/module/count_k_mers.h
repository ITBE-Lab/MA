/**
 * @file count_k_mers.h
 * @brief set of containers and modules to count k_mers/minimizers
 * @author Markus Schmidt
 */
#pragma once

#include "ma/container/minimizer_index.h"
#include "ma/container/nucSeq.h"
#include "ma/container/seed.h"
#include "ms/module/module.h"
#include <unordered_map>

using namespace libMS;
using namespace libMA;

// custom specialization of std::hash can be injected in namespace std
namespace std
{
template <size_t I> struct hash<std::array<uint8_t, I>>
{
    std::size_t operator( )( std::array<uint8_t, I> const& xBuff ) const noexcept
    {
        size_t uiRet = 0;
        for( size_t uiI = 0; uiI < I; uiI++ )
            uiRet = std::hash<uint8_t>{}( xBuff[ uiI ] ) ^ ( uiRet << 1 );
        return uiRet;
    }
};
} // namespace std

namespace libMSV
{

/**
 * @brief k-mer counting datastructure
 * @details
 * Performs thread-SAVE counting of k-mers.
 * Minimizes locking by using chunks and spin locks.
 */
template <nucSeqIndex NUM_CHUNKS, nucSeqIndex K> class __KMerCounter : public Container
{
    // make sure there is no empty space in the uint8_t arrays
    static_assert( ( K - NUM_CHUNKS ) % 4 == 0 );
    // make sure NUM_CHUNKS nt actually fit in the uint32_t index that is used
    static_assert( NUM_CHUNKS <= sizeof( uint32_t ) * 2 );
    class Chunk
    {
        using ARR_T = std::array<uint8_t, ( K - NUM_CHUNKS ) / 4>;
        std::unordered_map<ARR_T, size_t> xCountMap;
#define SPIN_LOCK 1
#if SPIN_LOCK
        std::atomic_flag xLock = ATOMIC_FLAG_INIT;
#else
        std::mutex xLock;
#endif
      public:
        Chunk( ) : xCountMap( )
        {} // constructor

        Chunk( const Chunk& rOther ) = delete;

        void inc( const NucSeq& xSeq, size_t uiCnt )
        {
#if SPIN_LOCK
            while( xLock.test_and_set( std::memory_order_acquire ) )
                ; // spin to acquire lock
#else
            std::lock_guard<std::mutex> xGuard( xLock );
#endif


            xCountMap[ toArray( xSeq ) ] += uiCnt;


#if SPIN_LOCK
            xLock.clear( std::memory_order_release ); // release lock
#endif
        } // method


        /* Set the value at position uiPosition in the packed array vArr.
         * Works only for the virtual forward strand.
         */
        inline void vSetNucleotideOnPos( ARR_T& vArr, const uint64_t uiPosition, const uint8_t uiValue )
        { /* We expect a correct position, when we come here.
           */
            // clear the two bits
            vArr[ ( size_t )( uiPosition >> 2 ) ] &= ~( 3 << ( ( ~uiPosition & 3UL ) << 1 ) );
            // add the new value
            vArr[ ( size_t )( uiPosition >> 2 ) ] |= uiValue << ( ( ~uiPosition & 3UL ) << 1 );
        } // inline method

        /* Get the value at position uiPosition in the unpacked sequence.
         * Works only for the virtual forward strand.
         */
        inline uint8_t getNucleotideOnPos( ARR_T& vArr, const uint64_t uiPosition )
        { /* We expect a correct position, when we come here.
           */
            return vArr[ ( size_t )( uiPosition >> 2 ) ] >> ( ( ~uiPosition & 3UL ) << 1 ) & 3;
        } // inline method

        ARR_T toArray( const NucSeq& xSeq )
        {
            ARR_T vRet{};
            for( size_t uiI = 0; uiI < K - NUM_CHUNKS; uiI++ )
                vSetNucleotideOnPos( vRet, uiI, xSeq[ uiI + NUM_CHUNKS ] );
            return vRet;
        } // method

        NucSeq fromArray( ARR_T vArr, uint32_t uiIdx )
        {
            NucSeq xRet;
            xRet.resize( K );
            for( size_t uiI = 0; uiI < NUM_CHUNKS; uiI++ )
            {
                xRet[ NUM_CHUNKS - ( uiI + 1 ) ] = uiIdx & 3;
                uiIdx >>= 2;
            } // for

            for( size_t uiI = 0; uiI < K - NUM_CHUNKS; uiI++ )
                xRet[ uiI + NUM_CHUNKS ] = getNucleotideOnPos( vArr, uiI );
            return xRet;
        } // method

        size_t get( const NucSeq& xSeq )
        {
            auto xIt = xCountMap.find( toArray( xSeq ) );
            if( xIt != xCountMap.end( ) )
                return xIt->second;
            return 0;
        } // method

        template <typename F> void iterate( F&& fDo, uint32_t uiIdx )
        {
            for( auto& xPair : xCountMap )
                fDo( fromArray( xPair.first, uiIdx ), xPair.second );
        } // method
    }; // class
  public:
    // array size should be 4 ^ NUM_CHUNKS which is 2 ^ (NUM_CHUNKS * 2)...
    std::array<Chunk, 2 << ( NUM_CHUNKS * 2 )> vChunks{};
    const nucSeqIndex uiW;

    __KMerCounter( nucSeqIndex uiW ) : uiW( uiW )
    {} // constructor

    __KMerCounter( const __KMerCounter& rOther ) = delete;

    __KMerCounter operator=( const __KMerCounter& rOther ) = delete;

    static uint32_t getChunkId( const NucSeq& xSection )
    {
        uint32_t uiRet = 0;
        for( size_t uiI = 0; uiI < NUM_CHUNKS && uiI < xSection.length( ); uiI++ )
            uiRet = ( uiRet << 2 ) ^ ( xSection[ uiI ] < 4 ? xSection[ uiI ] : 0 );
        assert( uiRet < ( 2 << ( NUM_CHUNKS * 2 ) ) );
        return uiRet;
    } // method

    template <typename F>
    static bool toKMers( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiFrom, nucSeqIndex uiTo, nucSeqIndex uiW, F&& fDo )
    {
        for( nucSeqIndex uiI = uiFrom; uiI <= uiTo - K; uiI += uiW )
        {
            NucSeq xSection( *pSeq, uiI, uiI + K );
            if( !fDo( xSection ) )
                return false;
        } // for
        return true;
    } // method

    template <typename F> bool toKMers( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiFrom, nucSeqIndex uiTo, F&& fDo )
    {
        return toKMers( pSeq, 0, pSeq->length( ), uiW, fDo );
    } // method

    template <typename F> bool toKMers( std::shared_ptr<NucSeq> pSeq, F&& fDo )
    {
        return toKMers( pSeq, 0, pSeq->length( ), fDo );
    } // method

    void addKMer( const NucSeq& xSeq, size_t uiCnt )
    {
        if( xSeq.length( ) != K )
            return;
        for( size_t uiI = 0; uiI < K; uiI++ )
            if( xSeq[ uiI ] > 3 )
                return;
        vChunks[ getChunkId( xSeq ) ].inc( xSeq, uiCnt );
    } // method

    void addSequence( std::shared_ptr<NucSeq> pSeq )
    {
        toKMers( pSeq, [&]( const NucSeq& xSeq ) {
            addKMer( xSeq, 1 );
            return true;
        } );
        pSeq->vReverseAll( );
        pSeq->vSwitchAllBasePairsToComplement( );
        toKMers( pSeq, [&]( const NucSeq& xSeq ) {
            addKMer( xSeq, 1 );
            return true;
        } );
        pSeq->vReverseAll( );
        pSeq->vSwitchAllBasePairsToComplement( );
    } // method

    bool isUnique( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiFrom, nucSeqIndex uiTo, nucSeqIndex uiMaxOcc )
    {
        for( size_t uiI = 0; uiI < pSeq->length( ); uiI++ )
            if( ( *pSeq )[ uiI ] > 3 )
                return false;
        return toKMers( pSeq, uiFrom, uiTo,
                        [&]( const NucSeq& xSeq ) { return vChunks[ getChunkId( xSeq ) ].get( xSeq ) <= uiMaxOcc; } );
    } // method

    bool isUnique( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiMaxOcc )
    {
        return isUnique( pSeq, 0, pSeq->length( ), uiMaxOcc );
    } // method

    bool isUnique( std::shared_ptr<NucSeq> pSeq )
    {
        return isUnique( pSeq, 0 );
    } // method

    template <typename F> void iterate( F&& fDo )
    {
        for( uint32_t uiI = 0; uiI < vChunks.size( ); uiI++ )
            vChunks[ uiI ].iterate( fDo, uiI );
    } // method
}; // class

using KMerCounter = __KMerCounter<6, 18>;

class GetKMerCounter : public Module<KMerCounter, false>
{
    const nucSeqIndex uiW;

  public:
    GetKMerCounter( const ParameterSetManager& rParameters, nucSeqIndex uiW ) : uiW( uiW )
    {} // constructor

    std::shared_ptr<KMerCounter> execute( )
    {
        return std::make_shared<KMerCounter>( uiW );
    } // method
}; // class

class KMerCounterModule : public Module<NucSeq, false, NucSeq, KMerCounter>
{
  public:
    KMerCounterModule( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<NucSeq> execute( std::shared_ptr<NucSeq> pSeq, std::shared_ptr<KMerCounter> pCounter )
    {
        pCounter->addSequence( pSeq );
        return pSeq;
    } // method
}; // class

class KMerCountFilterModule : public Module<Seeds, false, NucSeq, Seeds, KMerCounter>
{
  public:
    const nucSeqIndex uiMaxOcc;

    KMerCountFilterModule( const ParameterSetManager& rParameters, nucSeqIndex uiMaxOcc ) : uiMaxOcc( uiMaxOcc )
    {} // constructor

    std::shared_ptr<Seeds> execute( std::shared_ptr<NucSeq> pSeq, std::shared_ptr<Seeds> pSeeds,
                                    std::shared_ptr<KMerCounter> pCounter )
    {
        auto pRet = std::make_shared<Seeds>( );
        pRet->reserve( pSeeds->size( ) );
        for( auto& xSeed : *pSeeds )
            if( pCounter->isUnique( pSeq, xSeed.start( ), xSeed.end( ), uiMaxOcc ) )
                pRet->push_back( xSeed );
        return pRet;
    } // method
}; // class


/**
 * @brief k-mer counting datastructure
 * @details
 * Performs thread-SAVE counting of k-mers.
 * Minimizes locking by using chunks and spin locks.
 */
template <nucSeqIndex NUM_CHUNKS_BITS, nucSeqIndex K, typename hash_t> class __HashCounter : public Container
{
    static_assert( ( K - NUM_CHUNKS_BITS ) <= sizeof( hash_t ) * 8 );
    // make sure NUM_CHUNKS nt actually fit in the uint64_t index that is used
    static_assert( NUM_CHUNKS_BITS <= sizeof( uint64_t ) * 8 );
    class Chunk
    {
        std::unordered_map<hash_t, size_t> xCountMap;
#define SPIN_LOCK 1
#if SPIN_LOCK
        std::atomic_flag xLock = ATOMIC_FLAG_INIT;
#else
        std::mutex xLock;
#endif
      public:
        Chunk( ) : xCountMap( )
        {} // constructor

        Chunk( const Chunk& rOther ) = delete;

        void inc( hash_t uiHash, size_t uiCnt )
        {
#if SPIN_LOCK
            while( xLock.test_and_set( std::memory_order_acquire ) )
                ; // spin to acquire lock
#else
            std::lock_guard<std::mutex> xGuard( xLock );
#endif


            xCountMap[ uiHash ] += uiCnt;


#if SPIN_LOCK
            xLock.clear( std::memory_order_release ); // release lock
#endif
        } // method

        size_t get( hash_t uiHash )
        {
            auto xIt = xCountMap.find( uiHash );
            if( xIt != xCountMap.end( ) )
                return xIt->second;
            return 0;
        } // method

        template <typename F> void iterate( F&& fDo, uint64_t uiChunkBits )
        {
            for( auto& xPair : xCountMap )
                fDo( ( ( (uint64_t)xPair.first ) << NUM_CHUNKS_BITS ) | uiChunkBits, xPair.second );
        } // method
    }; // class
  public:
    // array size should be 2 ^ NUM_CHUNKS_BITS
    std::array<Chunk, 2 << NUM_CHUNKS_BITS> vChunks{};

    __HashCounter( )
    {} // constructor

    __HashCounter( const __HashCounter& rOther ) = delete;

    __HashCounter operator=( const __HashCounter& rOther ) = delete;

    size_t chunkId( uint64_t uiHash )
    {
        // clear all but the lower NUM_CHUNKS_BITS bits of uiHash
        return ( size_t )( ( uiHash << ( 8 * sizeof( uint64_t ) - NUM_CHUNKS_BITS ) ) >>
                           ( 8 * sizeof( uint64_t ) - NUM_CHUNKS_BITS ) );
    }

    hash_t chunkKey( uint64_t uiHash )
    {
        assert( ( uiHash >> NUM_CHUNKS_BITS ) < ( 2ull << K ) );
        // cut off the lower NUM_CHUNKS_BITS bits of uiHash
        return ( hash_t )( uiHash >> NUM_CHUNKS_BITS );
    }

    void addHash( uint64_t uiHash, size_t uiCnt )
    {
        vChunks[ chunkId( uiHash ) ].inc( chunkKey( uiHash ), uiCnt );
    } // method
    void addHash( uint64_t uiHash )
    {
        addHash( uiHash, 1 );
    } // method

    size_t get( uint64_t uiHash )
    {
        return vChunks[ chunkId( uiHash ) ].get( chunkKey( uiHash ) );
    } // method

    bool isUnique( uint64_t uiHash, nucSeqIndex uiMaxOcc )
    {
        return get( uiHash ) <= uiMaxOcc;
    } // method

    bool isUnique( uint64_t uiHash )
    {
        return isUnique( uiHash, 0 );
    } // method

    template <typename F> void iterate( F&& fDo )
    {
        for( uint64_t uiI = 0; uiI < vChunks.size( ); uiI++ )
            vChunks[ (size_t)uiI ].iterate( fDo, uiI );
    } // method
}; // class

using HashCounter = __HashCounter<16, 56, uint64_t>;


class MMCounterModule : public Module<NucSeq, false, NucSeq, HashCounter>
{
  public:
    nucSeqIndex k, w;
    MMCounterModule( const ParameterSetManager& rParameters )
        : k( rParameters.getSelected( )->xMinimizerK->get( ) ), w( rParameters.getSelected( )->xMinimizerW->get( ) )
    {} // constructor

    std::shared_ptr<NucSeq> execute( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<HashCounter> pCounter )
    {
        pQuery->vTranslateToCharacterForm( );
        const char* sSeq = (const char*)pQuery->pxSequenceRef;
        const int iSize = (int)pQuery->length( );
        for( auto uiHash : minimizer::Index::_getHash( sSeq, iSize, k, w ) )
            pCounter->addHash( uiHash );
        pQuery->vTranslateToNumericForm( );
        return pQuery;
    } // method
}; // class

class MMFilteredSeeding : public Module<Seeds, false, minimizer::Index, NucSeq, Pack, HashCounter>
{
  public:
    const nucSeqIndex uiMaxOcc;
    bool bRectangular;

    MMFilteredSeeding( const ParameterSetManager& rParameters )
        : uiMaxOcc( rParameters.getSelected( )->xMMFilterMaxOcc->get( ) ),
         bRectangular( rParameters.getSelected( )->xRectangularSoc->get( ) )
    {} // constructor

    std::shared_ptr<Seeds> execute( std::shared_ptr<minimizer::Index> pMMIndex, std::shared_ptr<NucSeq> pQuery,
                                    std::shared_ptr<Pack> pPack, std::shared_ptr<HashCounter> pCounter )
    {
        pQuery->vTranslateToCharacterForm( );
        const char* sSeq = (const char*)pQuery->pxSequenceRef;
        const int iSize = (int)pQuery->length( );
        // pack arguments for c function
        auto xPtr = std::make_pair( pCounter, uiMaxOcc );
        auto pRet = pMMIndex->seed_one(
            sSeq, iSize, bRectangular, pPack,
            // lambda function filters query minimizers before they are turned to seeds via hash table lookup
            [/*cannot capture since lambda needs to be passed to c as function pointer*/]( mm128_t* a, size_t& n,
                                                                                           void* pArg ) {
                // unpack arguments from c function
                auto pPair = static_cast<std::pair<std::shared_ptr<HashCounter>, nucSeqIndex>*>( pArg );
                size_t uiI = 0;
                while( uiI < n )
                {
                    if( !pPair->first->isUnique( minimizer::Index::_getHash( a[ uiI ] ), pPair->second ) )
                        a[ uiI ] = a[ --n ];
                    else
                        uiI++;
                } // while
            },
            // c function can only take void* as arguments
            static_cast<void*>( &xPtr ) );
        pQuery->vTranslateToNumericForm( );
        return pRet;
    } // method

    static std::vector<size_t> for_( const Seed& rS, std::shared_ptr<minimizer::Index> pMMIndex,
                                     std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pPack,
                                     std::shared_ptr<HashCounter> pCounter, bool bRectangular )
    {
        assert( rS.end( ) <= pQuery->length( ) );
        pQuery->vTranslateToCharacterForm( );
        const char* sSeq = (const char*)( pQuery->pxSequenceRef + rS.start( ) );
        const int iSize = (int)rS.size( );
        std::vector<size_t> vRet;
        // pack arguments for c function
        auto xPtr = std::make_pair( pCounter, &vRet );
        auto pRet = pMMIndex->seed_one(
            sSeq, iSize, bRectangular, pPack,
            // lambda function filters query minimizers before they are turned to seeds via hash table lookup
            [/*cannot capture since lambda needs to be passed to c as function pointer*/]( mm128_t* a, size_t& n,
                                                                                           void* pArg ) {
                // unpack arguments from c function
                auto pPair = static_cast<std::pair<std::shared_ptr<HashCounter>, std::vector<size_t>*>*>( pArg );
                for( size_t uiI = 0; uiI < n; uiI++ )
                    pPair->second->push_back( pPair->first->get( minimizer::Index::_getHash( a[ uiI ] ) ) );
            },
            // c function can only take void* as arguments
            static_cast<void*>( &xPtr ) );
        pQuery->vTranslateToNumericForm( );
        return vRet;
    } // method

    static size_t getMinCount( const Seed& rS, std::shared_ptr<minimizer::Index> pMMIndex,
                               std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pPack,
                               std::shared_ptr<HashCounter> pCounter, bool bRectangular )
    {
        auto vRet = for_( rS, pMMIndex, pQuery, pPack, pCounter, bRectangular );
        if( vRet.size( ) == 0 )
            return 0;
        size_t uiRet = vRet[ 0 ];
        for( auto x : vRet )
            uiRet = std::min( uiRet, x );
        return uiRet;
    } // method
    static size_t getMaxCount( const Seed& rS, std::shared_ptr<minimizer::Index> pMMIndex,
                               std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pPack,
                               std::shared_ptr<HashCounter> pCounter, bool bRectangular )
    {
        auto vRet = for_( rS, pMMIndex, pQuery, pPack, pCounter, bRectangular );
        if( vRet.size( ) == 0 )
            return 0;
        size_t uiRet = vRet[ 0 ];
        for( auto x : vRet )
            uiRet = std::max( uiRet, x );
        return uiRet;
    } // method
}; // class

} // namespace libMSV

#ifdef WITH_PYTHON
/**
 * @brief exports the modules to python
 * @ingroup export
 */
void exportCountKMers( libMS::SubmoduleOrganizer& xOrganizer );
#endif