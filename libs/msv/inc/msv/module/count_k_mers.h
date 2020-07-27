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

class KMerCounter : public Container
{
  public:
    std::unordered_map<ByteBuffer, size_t> xCountMap;
    const nucSeqIndex uiK;
    const nucSeqIndex uiW;
    size_t uiNumReads = 0;

    KMerCounter( nucSeqIndex uiK, nucSeqIndex uiW ) : uiK( uiK ), uiW( uiW )
    {} // constructor

    /// @note thread save
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

    /// @note thread save
    template <typename F> bool toKMers( std::shared_ptr<NucSeq> pSeq, F&& fDo )
    {
        return toKMers( pSeq, 0, pSeq->length( ), fDo );
    } // method

    /// @note NOT thread save
    void addSequence( std::shared_ptr<NucSeq> pSeq )
    {
        uiNumReads++;
        toKMers( pSeq, [&]( std::shared_ptr<CompressedNucSeq> pComp ) {
            xCountMap[ pComp->xCompSeqBuf ] += 1;
            return true;
        } );
        pSeq->vReverseAll( );
        pSeq->vSwitchAllBasePairsToComplement( );
        toKMers( pSeq, [&]( std::shared_ptr<CompressedNucSeq> pComp ) {
            xCountMap[ pComp->xCompSeqBuf ] += 1;
            return true;
        } );
        pSeq->vReverseAll( );
        pSeq->vSwitchAllBasePairsToComplement( );
    } // method

    /// @note NOT thread save
    void merge( std::shared_ptr<KMerCounter> pOther )
    {
        for( auto& xEle : pOther->xCountMap )
            xCountMap[ xEle.first ] += xEle.second;
    } // method

    /// @note NOT thread save
    void clear( )
    {
        xCountMap.clear( );
        uiNumReads = 0;
    } // method

    /// @note thread save
    bool isUnique( std::shared_ptr<NucSeq> pSeq, nucSeqIndex uiFrom, nucSeqIndex uiTo, nucSeqIndex uiMaxOcc )
    {
        return toKMers( pSeq, uiFrom, uiTo, [&]( std::shared_ptr<CompressedNucSeq> pComp ) {
            // make sure we do not accidentally insert an elemtent here (insert is not thread save)
            auto xIt = xCountMap.find( pComp->xCompSeqBuf );
            if( xIt == xCountMap.end( ) )
                return false;
            return xIt->second <= uiMaxOcc;
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