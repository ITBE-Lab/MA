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

namespace libMSV
{

class ExclusiveReadWrite
{
    const int iThreadFree = 0;
    const int iBlockedState = -1;
    std::atomic<int> iState;

  public:
    ExclusiveReadWrite( ) : iState( iThreadFree )
    {}

    template <typename F> void save_read( F&& fDo )
    {
        size_t uiTries = 0;
        while( true )
        {
            uiTries++;
            if( uiTries % 10 == 0 )
                std::this_thread::sleep_for( std::chrono::nanoseconds( uiTries / 10 ) );
            auto iLocalState = iState.load( );
            if( iLocalState != iBlockedState && iState.compare_exchange_weak( iLocalState, iLocalState + 1 ) )
                break;
        } // while
        fDo( );
        // atomic decrement
        iState--;
    } // method

    template <typename F> void save_write( F&& fDo )
    {
        // change iState from 0 to -1.
        // 0 means no thread is currently accessing
        // -1 means any other threads are blocked from acessing
        size_t uiTries = 0;
        while( true )
        {
            uiTries++;
            if( uiTries % 15 == 0 )
                std::this_thread::sleep_for( std::chrono::nanoseconds( uiTries / 15 ) );

            auto iExpectedState = iThreadFree;
            if( iState.compare_exchange_weak( iExpectedState, iBlockedState ) )
                break;
        } // while

        fDo( );
        iState = iThreadFree;
    } // method
}; // class

/**
 * @brief k-mer counting datastructure
 * @details
 * Performs lock-FREE thread-SAVE counting of k-mers, where any k-mer sequence is only ever stored once.
 * Hence it is memory efficient and multithreadable.
 * Datastructure is lock-free as long a no new k-mer sequence is inserted.
 * Meaning at the beginning of the insertion process this will lock threads regularly, but as soon as most possible
 * k-mers are inserted into here it will act lock-free.
 */
class KMerCounter : public Container, private ExclusiveReadWrite
{
    /**
     * @brief wrapper for std::atomic
     * @details
     * We need the copy contructor in order to stick Count in an unordered_map.
     * Whenever the map increases in size elements might get copied:
     * During this copy operation the atomic operations are broken, however, using ExclusiveReadWrite we
     * make sure that inserting a new element into xCountMap is only done if there is only one thread working
     * on the datastructure.
     * Read operations on xCountMap do never trigger move/copy operations, hence simultaneous increasing of the count
     * are valid.
     */
    class Count
    {
      public:
        std::atomic<size_t> uiCnt;

        Count( ) : uiCnt( 0 )
        {} // constructor

        Count( const Count& rOther ) : uiCnt( rOther.uiCnt.load( ) )
        {} // copy constructor
    }; // class
  public:
    std::unordered_map<std::shared_ptr<CompressedNucSeq>, Count> xCountMap;
    const nucSeqIndex uiK;
    const nucSeqIndex uiW;
    std::atomic<size_t> uiNumReads;

    KMerCounter( nucSeqIndex uiK, nucSeqIndex uiW ) : uiK( uiK ), uiW( uiW ), uiNumReads( 0 )
    {
        xCountMap.reserve( 100000 );
    } // constructor

    template <typename F>
    static bool toKMers( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiFrom, nucSeqIndex uiTo, nucSeqIndex uiK,
                         nucSeqIndex uiW, F&& fDo )
    {
        for( nucSeqIndex uiI = uiFrom; uiI < uiTo - uiK; uiI += uiW )
        {
            NucSeq xSection( *pSeq, uiI, uiI + uiK );
            auto pComp = std::make_shared<CompressedNucSeq>( );
            pComp->compress( xSection );
            if( !fDo( pComp ) )
                return false;
        } // for
        return true;
    } // method

    template <typename F> bool toKMers( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiFrom, nucSeqIndex uiTo, F&& fDo )
    {
        return toKMers( pSeq, 0, pSeq->length( ), uiK, uiW, fDo );
    } // method

    template <typename F> bool toKMers( std::shared_ptr<NucSeq> pSeq, F&& fDo )
    {
        return toKMers( pSeq, 0, pSeq->length( ), fDo );
    } // method

    void addKMer( const std::shared_ptr<CompressedNucSeq>& pComp, size_t uiCnt )
    {
        bool bExists;
        save_read( [&]( ) { bExists = xCountMap.count( pComp ) > 0; } );
        if( bExists )
            save_read( [&]( ) { xCountMap[ pComp ].uiCnt += uiCnt; } );
        else
            save_write( [&]( ) { xCountMap[ pComp ].uiCnt += uiCnt; } );
    }

    void addSequence( std::shared_ptr<NucSeq> pSeq )
    {
        uiNumReads++;
        std::vector<std::shared_ptr<CompressedNucSeq>> vAll;
        vAll.reserve( pSeq->length( ) * 2 );
        toKMers( pSeq, [&]( std::shared_ptr<CompressedNucSeq> pComp ) {
            vAll.push_back( pComp );
            return true;
        } );
        pSeq->vReverseAll( );
        pSeq->vSwitchAllBasePairsToComplement( );
        toKMers( pSeq, [&]( std::shared_ptr<CompressedNucSeq> pComp ) {
            vAll.push_back( pComp );
            return true;
        } );
        pSeq->vReverseAll( );
        pSeq->vSwitchAllBasePairsToComplement( );


        std::vector<std::shared_ptr<CompressedNucSeq>> vExisting;
        vExisting.reserve( pSeq->length( ) * 2 );
        std::vector<std::shared_ptr<CompressedNucSeq>> vNew;
        vNew.reserve( pSeq->length( ) * 2 );
        save_read( [&]( ) {
            for( auto pComp : vAll )
                if( xCountMap.count( pComp ) > 0 )
                    vExisting.push_back( pComp );
                else
                    vNew.push_back( pComp );
            for( auto pComp : vExisting )
                xCountMap.at( pComp ).uiCnt++;
        } );

        if( vNew.size( ) > 0 )
            save_write( [&]( ) {
                for( auto pComp : vNew )
                    xCountMap[ pComp ].uiCnt++;
            } );
    } // method

    void merge( std::shared_ptr<KMerCounter> pOther )
    {
        for( auto& xEle : pOther->xCountMap )
            addKMer( xEle.first, xEle.second.uiCnt.load( ) );
    } // method

    void clear( )
    {
        save_write( [&]( ) {
            xCountMap.clear( );
            uiNumReads = 0;
        } );
    } // method

    bool isUnique( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiFrom, nucSeqIndex uiTo, nucSeqIndex uiMaxOcc )
    {
        return toKMers( pSeq, uiFrom, uiTo, [&]( std::shared_ptr<CompressedNucSeq> pComp ) {
            // make sure we do not accidentally insert an elemtent here (insert is not thread save)
            size_t uiCnt = 0;
            save_read( [&]( ) {
                auto xIt = xCountMap.find( pComp );
                if( xIt != xCountMap.end( ) )
                    uiCnt = xIt->second.uiCnt.load( );
            } );
            return uiCnt <= uiMaxOcc;
        } );
    } // method

    /// @note thread save
    bool isUnique( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiMaxOcc )
    {
        return isUnique( pSeq, 0, pSeq->length( ), uiMaxOcc );
    } // method

}; // class

class GetKMerCounter : public Module<KMerCounter, false>
{
    const nucSeqIndex uiK;
    const nucSeqIndex uiW;

  public:
    GetKMerCounter( const ParameterSetManager& rParameters, nucSeqIndex uiK, nucSeqIndex uiW ) : uiK( uiK ), uiW( uiW )
    {} // constructor

    std::shared_ptr<KMerCounter> execute( )
    {
        return std::make_shared<KMerCounter>( uiK, uiW );
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