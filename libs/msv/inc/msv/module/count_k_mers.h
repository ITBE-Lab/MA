/**
 * @file count_k_mers.h
 * @brief set of containers and modules to count k_mers/minimizers
 * @author Markus Schmidt
 */
#pragma once

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
 * Tries to minimize locking by using chunks
 */
template <nucSeqIndex NUM_CHUNKS, nucSeqIndex K> class __KMerCounter : public Container
{
    class Chunk
    {
        using ARR_T = std::array<uint8_t, ( K - NUM_CHUNKS ) / 2>;
        std::unordered_map<ARR_T, size_t> xCountMap;
#define SPIN_LOCK 1
#if SPIN_LOCK
        std::atomic_flag xLock = ATOMIC_FLAG_INIT;
#else
        std::mutex xLock;
#endif

      public:
        Chunk( )
        {} // constructor
        Chunk( const Chunk& rOther ) = delete;

        void inc( const NucSeq& xSeq, size_t uiCnt )
        {
#if SPIN_LOCK
            while( xLock.test_and_set( std::memory_order_acquire ) ) // acquire lock
                ; // spin
#else
            std::lock_guard<std::mutex> xGuard( xLock );
#endif


            xCountMap[ toArray( xSeq ) ] += uiCnt;


#if SPIN_LOCK
            xLock.clear( std::memory_order_release ); // release lock
#endif
        } // method

        uint8_t comp( const NucSeq& xSeq, size_t uiI )
        {
            uint8_t uiA = xSeq[ 2 * uiI + NUM_CHUNKS ];
            if( uiA >= 4 )
                uiA = 0;
            uint8_t uiB = xSeq[ 2 * uiI + NUM_CHUNKS ];
            if( uiB >= 4 )
                uiB = 0;
            return ( uiA << 2 ) ^ uiB;
        }

        ARR_T toArray( const NucSeq& xSeq )
        {
            ARR_T vRet;
            for( size_t uiI = 0; uiI < ( K - NUM_CHUNKS ) / 2; uiI++ )
                vRet[ uiI ] = comp( xSeq, uiI );
            return vRet;
        }

        NucSeq fromArray( ARR_T vArr, size_t uiIdx )
        {
            NucSeq xRet;
            xRet.resize( K );
            for( size_t uiI = 0; uiI < NUM_CHUNKS; uiI++ )
            {
                xRet[ uiI ] = uiIdx ^ 3;
                uiIdx >>= 2;
            } // for

            for( size_t uiI = 0; uiI < K - NUM_CHUNKS; uiI += 2 )
            {
                xRet[ uiI + NUM_CHUNKS ] = vArr[ uiI ] ^ 3;
                xRet[ uiI + NUM_CHUNKS + 1 ] = ( vArr[ uiI ] >> 2 ) ^ 3;
            } // for
            return xRet;
        }

        size_t get( const NucSeq& xSeq )
        {
            auto xIt = xCountMap.find( toArray( xSeq ) );
            if( xIt != xCountMap.end( ) )
                return xIt->second;
            return 0;
        } // method

        template <typename F> void iterate( F&& fDo, size_t uiIdx )
        {
            for( auto& xPair : xCountMap )
                fDo( fromArray( xPair.first, uiIdx ), xPair.second );
        } // method
    }; // class
  public:
    // array size should be 4 ^ NUM_CHUNKS which is 2 ^ (NUM_CHUNKS * 2)...
    std::array<Chunk, 2 << ( NUM_CHUNKS * 2 )> vChunks;
    const nucSeqIndex uiW;

    __KMerCounter( nucSeqIndex uiW ) : uiW( uiW )
    {} // constructor

    __KMerCounter( const __KMerCounter& rOther ) = delete;

    __KMerCounter operator=( const __KMerCounter& rOther ) = delete;

    static size_t getChunkId( const NucSeq& xSection )
    {
        size_t uiRet = 0;
        for( size_t uiI = 0; uiI < NUM_CHUNKS && uiI < xSection.length( ); uiI++ )
            uiRet = ( uiRet << 2 ) ^ ( xSection[ uiI ] < 4 ? xSection[ uiI ] : 0 );
        return uiRet;
    }

    template <typename F>
    static bool toKMers( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiFrom, nucSeqIndex uiTo, nucSeqIndex uiW, F&& fDo )
    {
        for( nucSeqIndex uiI = uiFrom; uiI < uiTo - K; uiI += uiW )
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
        vChunks[ getChunkId( xSeq ) ].inc( xSeq, uiCnt );
    }

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
        return toKMers( pSeq, uiFrom, uiTo,
                        [&]( const NucSeq& xSeq ) { return vChunks[ getChunkId( xSeq ) ].get( xSeq ) <= uiMaxOcc; } );
    } // method

    bool isUnique( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiMaxOcc )
    {
        return isUnique( pSeq, 0, pSeq->length( ), uiMaxOcc );
    } // method

    template <typename F> void iterate( F&& fDo )
    {
        for( size_t uiI = 0; uiI < vChunks.size(); uiI++ )
            vChunks[ uiI ].iterate( fDo, uiI );
    } // method
}; // class

using KMerCounter = __KMerCounter<8, 18>;

class GetKMerCounter : public Module<KMerCounter, false>
{
    const nucSeqIndex uiK;
    const nucSeqIndex uiW;

  public:
    GetKMerCounter( const ParameterSetManager& rParameters, nucSeqIndex uiK, nucSeqIndex uiW ) : uiK( uiK ), uiW( uiW )
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

} // namespace libMSV

#ifdef WITH_PYTHON
/**
 * @brief exports the modules to python
 * @ingroup export
 */
void exportCountKMers( libMS::SubmoduleOrganizer& xOrganizer );
#endif