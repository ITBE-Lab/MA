/**
 * @file computeCoverage.h
 * @brief Computes the coverage for SV-calls
 * @author Markus Schmidt
 */
#pragma once

#include "container/segment.h"
#include "container/svDb.h"
#include "container/svJump.h"
#include "contrib/intervalTree/IntervalTree.h"
#include "module/module.h"

namespace libMA
{

/**
 * @brief Computes the coverage for SV-calls
 * @details
 * Uses an interval tree internally
 * @ingroup module
 */
class ComputeCoverage : public Module<Container, false, SuffixArrayInterface, NucSeq, SegmentVector>
{
    class SvCallWrapper
    {
      public:
        SvCall xCall;
        std::array<std::vector<Seed>, 2> avCoverageAnalysis;
        std::mutex xMutex;

        SvCallWrapper( SvCall xCall ) : xCall( xCall ), avCoverageAnalysis( ), xMutex( )
        {} // constructor

        SvCallWrapper( SvCallWrapper& rOther ) = delete;
        SvCallWrapper& operator=( SvCallWrapper& rOther ) = delete;

        inline bool checkInterval( nucSeqIndex uiStart, nucSeqIndex uiEnd, Seed& rS, int64_t iReadId, bool bIsFrom,
                                   const size_t uiMaxDistance, const size_t uiAllowedOverlap )
        {
            if( rS.start_ref( ) > uiEnd + uiMaxDistance || rS.end_ref( ) + uiMaxDistance < uiStart )
                return false;

            bool bIsLeft = rS.start_ref( ) + uiAllowedOverlap < uiStart;
            bool bIsRight = rS.end_ref( ) > uiEnd + uiAllowedOverlap;


            if( bIsLeft )
                avCoverageAnalysis[ xCall.bSwitchStrand && !bIsFrom ? 1 : 0 ].push_back( rS );
            if( bIsRight )
                avCoverageAnalysis[ xCall.bSwitchStrand && !bIsFrom ? 0 : 1 ].push_back( rS );

            return true;
        } // method

        inline void addSeed( Seed& rS, std::shared_ptr<NucSeq> pQuery, const size_t uiMaxDistance,
                             const size_t uiAllowedOverlap )
        {
            std::lock_guard<std::mutex> xGuard( xMutex );
            size_t uiNumFound = 0;

            if( this->checkInterval( this->xCall.uiFromStart, this->xCall.uiFromStart + this->xCall.uiFromSize, rS,
                                     pQuery->iId, true, uiMaxDistance, uiAllowedOverlap ) )
                uiNumFound += 1;
            if( this->checkInterval( this->xCall.uiToStart, this->xCall.uiToStart + this->xCall.uiToSize, rS,
                                     pQuery->iId, false, uiMaxDistance, uiAllowedOverlap ) )
                uiNumFound += 1;
            // if this did neither overlap the start nor the end interval something is wrong wiht the compute_coverage
            // function below
            assert( uiNumFound > 0 );
        } // method

        inline void compute_coverage( )
        {
            xCall.uiCoverage = 1;
            for( size_t uiIdx : {0, 1} )
                if( avCoverageAnalysis[ uiIdx ].size( ) > 0 )
                {
                    // sort seeds by reference position
                    std::sort( avCoverageAnalysis[ uiIdx ].begin( ), avCoverageAnalysis[ uiIdx ].end( ),
                               []( Seed& rA, Seed& rB ) { return rA.start_ref( ) < rB.start_ref( ); } );

                    std::vector<std::pair<size_t, nucSeqIndex>> vCoverageList;

                    std::vector<Seed>::iterator itStart = avCoverageAnalysis[ uiIdx ].begin( );
                    std::vector<Seed>::iterator itEnd = avCoverageAnalysis[ uiIdx ].begin( );
                    size_t uiCoverageCnt = 0;
                    size_t uiCoverageSum = 0;
                    nucSeqIndex uiLastPos = itStart->start_ref( );

                    while( itStart != avCoverageAnalysis[ uiIdx ].end( ) ||
                           itEnd != avCoverageAnalysis[ uiIdx ].end( ) )
                    {
                        if( itStart != avCoverageAnalysis[ uiIdx ].end( ) && itStart->start_ref( ) < itEnd->end_ref( ) )
                        {
                            if( itStart->start_ref( ) > uiLastPos )
                            {
                                vCoverageList.emplace_back( uiCoverageCnt, uiLastPos - itStart->start_ref( ) );
                                uiCoverageSum += uiLastPos - itStart->start_ref( );
                            } // if
                            uiCoverageCnt++;
                            uiLastPos = itStart->start_ref( );
                            itStart++;
                        } // if
                        else
                        {
                            if( itEnd->end_ref( ) > uiLastPos )
                            {
                                vCoverageList.emplace_back( uiCoverageCnt, uiLastPos - itEnd->end_ref( ) );
                                uiCoverageSum += uiLastPos - itEnd->end_ref( );
                            } // if
                            uiCoverageCnt--;
                            uiLastPos = itEnd->end_ref( );

                            itEnd++;
                        } // else
                    } // while

                    std::sort( vCoverageList.begin( ), vCoverageList.end( ) );
                    size_t uiIdx2 = 0;
                    while( uiIdx2 + 1 < vCoverageList.size( ) && uiCoverageSum > vCoverageList[ uiIdx2 ].second * 2 )
                        uiCoverageSum -= vCoverageList[ uiIdx2++ ].second * 2;

                    assert( vCoverageList[ uiIdx2 ].first < 1000 );
                    if( vCoverageList[ uiIdx2 ].first > xCall.uiCoverage )
                        xCall.uiCoverage = (uint32_t)vCoverageList[ uiIdx2 ].first;
                } // if
        } // method
    }; // class

  public:
    std::shared_ptr<SV_DB> pDb;
    const unsigned int uiMaxAmbiguity;
    ///
    const size_t uiMinSeedSize;
    const size_t uiMaxDistance = 50;
    const size_t uiAllowedOverlap = 3;

    std::shared_ptr<interval_tree::IntervalTree<nucSeqIndex, std::shared_ptr<SvCallWrapper>>> pIntervalTree;
    std::vector<std::shared_ptr<SvCallWrapper>> vCalls;

    /**
     * @brief Initialize a BinarySeeding Module
     * @details
     * if bLrExtension is True our extension scheme is used,
     * otherwise the extension scheme by Li et Al. is used.
     * Our approach is faster and computes seeds of higher quality.
     * However Li et Al.s approach will increase the overall accuracy of the alignment
     * for short queries.
     */
    ComputeCoverage( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iCallerId,
                     int64_t iSequencerId )
        : pDb( pDb ),
          uiMaxAmbiguity( rParameters.getSelected( )->xMaxAmbiguitySv->get( ) ),
          uiMinSeedSize( rParameters.getSelected( )->xMinSeedSizeSV->get( ) )
    {
        SvCallsFromDb xCallGetter( rParameters, pDb, iCallerId );
        std::vector<interval_tree::Interval<nucSeqIndex, std::shared_ptr<SvCallWrapper>>> vIntervalVec;
        while( xCallGetter.hasNext( ) )
        {
            auto xCall = xCallGetter.next( );
            vCalls.push_back( std::make_shared<SvCallWrapper>( xCall ) );
            vIntervalVec.emplace_back( xCall.uiFromStart, xCall.uiFromStart + xCall.uiFromSize, vCalls.back( ) );
            vIntervalVec.emplace_back( xCall.uiToStart, xCall.uiToStart + xCall.uiToSize, vCalls.back( ) );
        } // while
        pIntervalTree = std::make_shared<interval_tree::IntervalTree<nucSeqIndex, std::shared_ptr<SvCallWrapper>>>(
            std::move( vIntervalVec ) );
    } // constructor

    ~ComputeCoverage( )
    {
        for( std::shared_ptr<SvCallWrapper> pWrap : vCalls )
        {
            pWrap->compute_coverage( );
            pDb->updateCoverage( pWrap->xCall );
        } // for
    } // deconstructor

    virtual std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<SuffixArrayInterface> pFM_index,
                                                         std::shared_ptr<NucSeq>
                                                             pQuerySeq,
                                                         std::shared_ptr<SegmentVector>
                                                             pSeeds )
    {
        pSeeds->forEachSeed(
            *std::dynamic_pointer_cast<FMIndex>( pFM_index ), pQuerySeq->length( ), uiMaxAmbiguity, uiMinSeedSize, true,
            [&]( Seed xS ) //
            { //
                pIntervalTree->visit_overlapping(
                    uiMaxDistance > xS.start_ref( ) ? 0 : xS.start_ref( ) - uiMaxDistance,
                    xS.end_ref( ) + uiMaxDistance,
                    [&]( const interval_tree::Interval<nucSeqIndex, std::shared_ptr<SvCallWrapper>>& rInterval ) //
                    { //
                        rInterval.value->addSeed( xS, pQuerySeq, uiMaxDistance, uiAllowedOverlap );
                    } // lambda
                ); // visit_overlapping function call
                return true;
            } // lambda
        ); // forEachSeed function call

        return std::make_shared<Container>( );
    } // method

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief exports the Segmentation @ref Module "module" to python.
 * @ingroup export
 */
void exportComputeCoverage( py::module& rxPyModuleId );
#endif