/**
 * @file sweepSvJumps.h
 * @author Markus Schmidt
 */
#pragma once

#include "container/squeezedVector.h"
#include "container/svDb.h"
#include "module/module.h"

namespace libMA
{

class GenomeSection : public Container, public Interval<int64_t>
{
    using Interval<int64_t>::Interval;
}; // class

/**
 * @brief generates evenly spaced intervals over the length of the given pack
 * @details used for parallel implementation of the complete bipartite subgraph (CBSG) sweep.
 * last segment will most likely extend over the end of the genome.
 */
class GenomeSectionFactory : public Module<GenomeSection, true>
{
  public:
    int64_t iRefSize;
    int64_t iSectionSize;
    int64_t iCurrStart;
    /**
     * @brief
     * @details
     */
    GenomeSectionFactory( const ParameterSetManager& rParameters, std::shared_ptr<Pack> pPack )
        : iRefSize( (int64_t)pPack->uiStartOfReverseStrand( ) ),
          // compute the number of genome sections so that there are 100 sections for each thread
          // 50 * 2 because forw & rev. further the sections should all be at least 10000nt long
          // otherwise we do so much extra work with re overlapping part between sections,
          // that parallel execution is not worth it.
          iSectionSize( std::max( iRefSize / ( int64_t )( rParameters.getNumThreads( ) * 50 ), (int64_t)100000 ) ),
          iCurrStart( 0 )
    {} // constructor

    virtual std::shared_ptr<GenomeSection> EXPORTED execute( )
    {
        // setFinished( );
        // return std::make_shared<GenomeSection>( 0, std::numeric_limits<int64_t>::max( ) - 10000 );

        std::shared_ptr<GenomeSection> pRet;
        if( iCurrStart % 2 == 0 ) // forward strand
            pRet = std::make_shared<GenomeSection>( ( iCurrStart / 2 ) * iSectionSize, iSectionSize );
        else // reverse strand
            pRet = std::make_shared<GenomeSection>(
                ( iCurrStart / 2 ) * iSectionSize + std::numeric_limits<int64_t>::max( ) / (int64_t)2, iSectionSize );

        iCurrStart++;
        if( ( iCurrStart / 2 ) * iSectionSize >= iRefSize )
            setFinished( );
        return pRet;
    } // method

    virtual bool requiresLock( ) const
    {
        return true;
    } // method
}; // class


class CompleteBipartiteSubgraphClusterVector : public Container
{
  public:
    std::vector<std::shared_ptr<SvCall>> vContent;
}; // class

/**
 * @brief saves all computed clusters in the database
 * @details
 */
class SvCallSink : public Module<Container, false, CompleteBipartiteSubgraphClusterVector>
{
  public:
    std::shared_ptr<SV_DB::SvCallInserter> pInserter;
    std::mutex xLock;
    /**
     * @brief
     * @details
     */
    SvCallSink( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDB, std::string rsSvCallerName,
                std::string rsSvCallerDesc, int64_t uiJumpRunId )
        : pInserter( std::make_shared<SV_DB::SvCallInserter>( pDB, rsSvCallerName, rsSvCallerDesc, uiJumpRunId ) )
    {} // constructor

    virtual std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<CompleteBipartiteSubgraphClusterVector> pVec )
    {
        std::lock_guard<std::mutex> xGuard( xLock );
        for( auto pCall : pVec->vContent )
            pInserter->insertCall( *pCall );
        return std::make_shared<Container>( );
    } // method
}; // class

/**
 * @brief
 * @details
 */
class CompleteBipartiteSubgraphSweep : public Module<CompleteBipartiteSubgraphClusterVector, false, GenomeSection>
{
  public:
    const ParameterSetManager& rParameters;
    std::shared_ptr<SV_DB> pSvDb;
    std::shared_ptr<Pack> pPack;
    int64_t iSvCallerRunId;
    int64_t iMaxFuzziness;
    nucSeqIndex uiGenomeSize;
    size_t uiSqueezeFactor;
    size_t uiCenterStripUp;
    size_t uiCenterStripDown;
    double dEstimateCoverageFactor;
    std::vector<double> vEstimatedCoverageList;
    /**
     * @brief
     * @details
     */
    CompleteBipartiteSubgraphSweep( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pSvDb,
                                    std::shared_ptr<Pack> pPack, int64_t iSvCallerRunId )
        : rParameters( rParameters ),
          pSvDb( pSvDb ),
          pPack( pPack ),
          iSvCallerRunId( iSvCallerRunId ),
          // @todo this does not consider tail edges (those should be limited in size and then here we should use the
          // max of both limits)
          // also this should be the maximal cluster width not the maximal CBSG width
          iMaxFuzziness( (int64_t)rParameters.getSelected( )->xJumpH->get( ) ),
          uiGenomeSize( pPack->uiStartOfReverseStrand( ) ),
          uiSqueezeFactor( 5000 ),
          uiCenterStripUp( 5000 ),
          uiCenterStripDown( 1000 ),
          dEstimateCoverageFactor( 0.1 ),
          vEstimatedCoverageList( pSvDb->pContigCovTable->getEstimatedCoverageList( iSvCallerRunId, pPack ) )
    {} // constructor

    virtual std::shared_ptr<CompleteBipartiteSubgraphClusterVector>
        EXPORTED execute( std::shared_ptr<GenomeSection> pSection )
    {
        SortedSvJumpFromSql xEdges(
            rParameters, pSvDb, iSvCallerRunId,
            // make sure we overlap the start of the next interval, so that clusters that span over two intervals
            // are being collected. -> for this we just keep going after the end of the interval
            pSection->iStart > iMaxFuzziness ? pSection->iStart - iMaxFuzziness : 0,
            pSection->iSize + iMaxFuzziness * 5 );

        nucSeqIndex uiForwStrandStart = (nucSeqIndex)pSection->start( );
        nucSeqIndex uiForwStrandEnd = (nucSeqIndex)pSection->end( );
        if( pSection->iStart >= std::numeric_limits<int64_t>::max( ) / (int64_t)2 )
        {
            uiForwStrandStart = ( nucSeqIndex )( pSection->start() - std::numeric_limits<int64_t>::max( ) / (int64_t)2 );
            uiForwStrandEnd = ( nucSeqIndex )( pSection->end( ) - std::numeric_limits<int64_t>::max( ) / (int64_t)2 );
        } // if

        SqueezedVector<std::shared_ptr<SvCall>> xPointerVec( uiGenomeSize, uiSqueezeFactor, uiCenterStripUp,
                                                             uiCenterStripDown );

        auto pRet = std::make_shared<CompleteBipartiteSubgraphClusterVector>( );

#if DEBUG_LEVEL > 0
        std::set<int64_t> xVisitedStart;
        std::vector<std::shared_ptr<SvCall>> xActiveClusters;
#endif

        while( xEdges.hasNextStart( ) || xEdges.hasNextEnd( ) )
        {
            if( xEdges.nextStartIsSmaller( ) )
            {
                auto pEdge = xEdges.getNextStart( );
#if DEBUG_LEVEL > 0
                xVisitedStart.insert( pEdge->iId );
#endif
                auto pNewCluster = std::make_shared<SvCall>( pEdge );

                size_t uiStart = xPointerVec.to_physical_coord( pNewCluster->uiFromStart + pNewCluster->uiFromSize - 1,
                                                                pNewCluster->uiToStart );
                size_t uiEnd = xPointerVec.to_physical_coord( pNewCluster->uiFromStart,
                                                              pNewCluster->uiToStart + pNewCluster->uiToSize - 1 );
                // set the clusters y coodinate to the physical coords (we won't use the actual coords anyways)
                // this is necessary, since we need to work with these coords when joining clusters
                pNewCluster->uiToStart = uiStart;
                pNewCluster->uiToSize = uiEnd - uiStart;
                assert( uiEnd >= uiStart );

                // join with all covered clusters; make sure that we don't join the same cluster twice
                std::shared_ptr<SvCall> pLastJoined;
                for( size_t uiI = uiStart; uiI <= uiEnd; uiI++ )
                    if( xPointerVec.get( )[ uiI ] != nullptr && pLastJoined != xPointerVec.get( )[ uiI ] )
                    {
                        pLastJoined = xPointerVec.get( )[ uiI ];
#if DEBUG_LEVEL > 0
                        assert( pLastJoined->uiOpenEdges > 0 );
                        xActiveClusters.erase( std::remove_if( xActiveClusters.begin( ), xActiveClusters.end( ),
                                                               [&]( auto pX ) { return pX == pLastJoined; } ),
                                               xActiveClusters.end( ) );
#endif
                        pNewCluster->join( *pLastJoined );
                    } // if

#if DEBUG_LEVEL > 0
                xActiveClusters.push_back( pNewCluster );
#endif
                // redirect all covered pointers to the new cluster
                for( size_t uiI = pNewCluster->uiToStart; uiI <= pNewCluster->uiToStart + pNewCluster->uiToSize; uiI++ )
                {
#if DEBUG_LEVEL > 0
                    if( xPointerVec.get( )[ uiI ] != nullptr )
                        for( int64_t iId : xPointerVec.get( )[ uiI ]->vSupportingJumpIds )
                            assert( std::find( pNewCluster->vSupportingJumpIds.begin( ),
                                               pNewCluster->vSupportingJumpIds.end( ),
                                               iId ) != pNewCluster->vSupportingJumpIds.end( ) );
#endif
                    xPointerVec.get( )[ uiI ] = pNewCluster;
                } // for
            } // if
            else
            {
                auto pEndJump = xEdges.getNextEnd( );
#if DEBUG_LEVEL > 0
                assert( xVisitedStart.count( pEndJump->iId ) != 0 );
#endif
                // find the correct cluster for this edge
                auto pCluster = xPointerVec.get( )[ xPointerVec.to_physical_coord(
                    pEndJump->from_start_same_strand( ) + pEndJump->from_size( ) - 1, pEndJump->to_start( ) ) ];
                assert( pCluster != nullptr );
                assert( std::find( pCluster->vSupportingJumpIds.begin( ), pCluster->vSupportingJumpIds.end( ),
                                   pEndJump->iId ) != pCluster->vSupportingJumpIds.end( ) );
                pCluster->uiOpenEdges--;
                // check if we want to save the cluster
                if( pCluster->uiOpenEdges == 0 )
                {
                    for( size_t uiI = pCluster->uiToStart; uiI <= pCluster->uiToStart + pCluster->uiToSize; uiI++ )
                        xPointerVec.get( )[ uiI ] = nullptr;

#if DEBUG_LEVEL > 0
                    xActiveClusters.erase( std::remove_if( xActiveClusters.begin( ), xActiveClusters.end( ),
                                                           [&]( auto pX ) { return pX == pCluster; } ),
                                           xActiveClusters.end( ) );
#endif
                    if( pCluster->uiFromStart < uiForwStrandEnd && pCluster->uiFromStart >= uiForwStrandStart )
                    {
                        double estimatedCoverage =
                            std::min( vEstimatedCoverageList[ pPack->uiSequenceIdForPosition( pCluster->uiFromStart ) ],
                                      vEstimatedCoverageList[ pPack->uiSequenceIdForPosition(
                                          pCluster->vSupportingJumps[ 0 ]->to_start( ) ) ] );
                        if( pCluster->estimateCoverage( ) >=
                            std::max( dEstimateCoverageFactor * estimatedCoverage, 2.0 ) )
                            pRet->vContent.push_back( pCluster );
                    } // if
                } // if
            } // else
        } // while

#if DEBUG_LEVEL > 0
        // make sure that there is no open cluster left
        for( auto pCluster : xPointerVec.get( ) )
            assert( pCluster == nullptr );
#endif

        return pRet;
    } // method
}; // class

/**
 * @brief
 * @details
 */
class ExactCompleteBipartiteSubgraphSweep
    : public Module<CompleteBipartiteSubgraphClusterVector, false, CompleteBipartiteSubgraphClusterVector>
{
  public:
    std::shared_ptr<Pack> pPack;
    double dEstimateCoverageFactor;
    std::vector<double> vEstimatedCoverageList;
    int64_t iMaxInsertRatioDiff = 150;
    /**
     * @brief
     * @details
     */
    ExactCompleteBipartiteSubgraphSweep( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pSvDb,
                                         std::shared_ptr<Pack> pPack, int64_t iSvCallerRunId )
        : pPack( pPack ),
          dEstimateCoverageFactor( 0.1 ),
          vEstimatedCoverageList( pSvDb->pContigCovTable->getEstimatedCoverageList( iSvCallerRunId, pPack ) )
    {} // constructor

    inline void exact_sweep( std::vector<std::shared_ptr<SvJump>>& rvEdges, size_t uiStart, size_t uiEnd,
                             std::shared_ptr<CompleteBipartiteSubgraphClusterVector> pRet )
    {
        // squash sv_jump indices
        std::map<int64_t, size_t> xSquashedY;
        std::vector<std::shared_ptr<SvJump>> vEdgesStart;
        std::vector<std::shared_ptr<SvJump>> vEdgesEnd;
        while( uiStart < uiEnd )
        {
            std::shared_ptr<SvJump> pJmp = rvEdges[ uiStart ];
            xSquashedY[ pJmp->to_start( ) ] = 0;
            xSquashedY[ pJmp->sweep_end( ) ] = 0;
            vEdgesStart.push_back( pJmp );
            vEdgesEnd.push_back( pJmp );
            uiStart++;
        } // for

        // set the indices in the squashed map correctly
        {
            size_t uiI = 0;
            for( auto xIter = xSquashedY.begin( ); xIter != xSquashedY.end( ); ++xIter )
            {
                xIter->second = uiI;
                uiI++;
            } // for
        } // scope

        // create start list
        std::sort( vEdgesStart.begin( ), vEdgesStart.end( ),
                   []( std::shared_ptr<SvJump> pA, std::shared_ptr<SvJump> pB ) {
                       return pA->from_start( ) < pB->from_start( );
                   } ); // sort function call

        // create end list
        std::sort( vEdgesEnd.begin( ), vEdgesEnd.end( ), []( std::shared_ptr<SvJump> pA, std::shared_ptr<SvJump> pB ) {
            return pA->from_end( ) < pB->from_end( );
        } ); // sort function call

        // initialize the y-axis sweep pointer and counter vector
        std::vector<std::pair<std::shared_ptr<SvCall>, size_t>> vSweepVec( xSquashedY.size( ) );

        // do the actual sweep
        size_t uiI = 0;
        size_t uiJ = 0;
        while( uiJ < vEdgesEnd.size( ) )
        {
            if( uiI < vEdgesStart.size( ) && vEdgesStart[ uiI ]->from_start( ) <= vEdgesEnd[ uiJ ]->from_end( ) )
            {
                // create a cluster containing merely the current call
                auto pNewCluster = std::make_shared<SvCall>( vEdgesStart[ uiI ] );

                // turn tail edge lines into squares
                if( !vEdgesStart[ uiI ]->switch_strand_known( ) )
                    pNewCluster->uiToSize = vEdgesStart[ uiI ]->from_size( ) + 1;

                // join with all overlapping clusters
                size_t uiStart = xSquashedY[ vEdgesStart[ uiI ]->to_start( ) ];
                size_t uiIdx = uiStart;
                size_t uiEnd = xSquashedY[ vEdgesStart[ uiI ]->sweep_end( ) ];
                std::set<std::shared_ptr<SvCall>> xJoinedClusters;
                while( uiIdx <= uiEnd )
                {
                    if( vSweepVec[ uiIdx ].second > 0 &&
                        xJoinedClusters.find( vSweepVec[ uiIdx ].first ) == xJoinedClusters.end( ) )
                    {
                        pNewCluster->join( *vSweepVec[ uiIdx ].first );
                        xJoinedClusters.insert( vSweepVec[ uiIdx ].first );
                    } // if
                    uiIdx++;
                } // while

                // insert the newly computed cluster into the pointer vector and counter vector
                uiIdx = xSquashedY[ pNewCluster->uiToStart ];
                size_t uiEnd2 = xSquashedY[ pNewCluster->uiToStart + pNewCluster->uiToSize - 1 ];
                while( uiIdx <= uiEnd2 )
                {
                    if( uiStart <= uiIdx && uiIdx <= uiEnd )
                    {
                        vSweepVec[ uiIdx ].second++;
                        vSweepVec[ uiIdx ].first = pNewCluster;
                    } // if
                    else if( vSweepVec[ uiIdx ].second > 0 &&
                             xJoinedClusters.find( vSweepVec[ uiIdx ].first ) != xJoinedClusters.end( ) )
                        vSweepVec[ uiIdx ].first = pNewCluster;
                    uiIdx++;
                } // while

                uiI++;
            } // if
            else
            {
                size_t uiStart = xSquashedY[ vEdgesEnd[ uiJ ]->to_start( ) ];
                size_t uiEnd = xSquashedY[ vEdgesEnd[ uiJ ]->sweep_end( ) ];
                auto pCurrCluster = vSweepVec[ uiStart ].first;
                assert( pCurrCluster != nullptr );
                pCurrCluster->uiOpenEdges -= 1;
                // check if that closes the cluster
                if( pCurrCluster->uiOpenEdges == 0 )
                {
                    // remove jumps with equal read id's
                    std::map<int64_t, std::shared_ptr<SvJump>> xNewJumps;
                    for( auto pJump : pCurrCluster->vSupportingJumps )
                    {
                        if( xNewJumps.count( pJump->iReadId ) == 0 ||
                            pJump->query_distance( ) < xNewJumps[ pJump->iReadId ]->query_distance( ) )
                            xNewJumps[ pJump->iReadId ] = pJump;
                    } // for
                    pCurrCluster->vSupportingJumps.clear( );
                    for( auto& xTup : xNewJumps )
                        pCurrCluster->vSupportingJumps.push_back( xTup.second );

                    // check if the cluster still fulfills the required criteria
                    double estimatedCoverage =
                        std::min( vEstimatedCoverageList[ pPack->uiSequenceIdForPosition( pCurrCluster->uiFromStart ) ],
                                  vEstimatedCoverageList[ pPack->uiSequenceIdForPosition( pCurrCluster->uiToStart ) ] );
                    if( pCurrCluster->estimateCoverage( ) >=
                        std::max( dEstimateCoverageFactor * estimatedCoverage, 2.0 ) )
                    {
                        pCurrCluster->reEstimateClusterSize( );
                        // we messed up this counter by removing jumps, fix that
                        pCurrCluster->uiNumSuppNt = 0;
                        for( auto pJump : pCurrCluster->vSupportingJumps )
                            pCurrCluster->uiNumSuppNt += pJump->numSupportingNt( );
                        // save the cluster
                        pRet->vContent.push_back( pCurrCluster );
                    } // if
                } // if
                // decrement the counter vector
                while( uiStart <= uiEnd )
                    vSweepVec[ uiStart++ ].second--;
                uiJ++;
            } // else
        } // while


    } // method

    /// complete linkage clustering for jump distances
    /// we call sweep_sv_jumps for all insert_ratio clusters with a max dist of max_insert_ratio_diff
    /// this clustering is necessary because there might be an edge in the graph that has several different inserted
    /// sequences. We need to consider these sequences individually -> cluster by sequence length
    /// if the sequences are different by nucleotides, we need to figure that out later in the multialignment step...
    inline void lineSweep( std::shared_ptr<SvCall> pCluster, //
                           std::shared_ptr<CompleteBipartiteSubgraphClusterVector>
                               pRet )
    {
        std::vector<std::shared_ptr<SvJump>>& rvEdges = pCluster->vSupportingJumps;
        std::sort( rvEdges.begin( ), rvEdges.end( ), []( std::shared_ptr<SvJump> pA, std::shared_ptr<SvJump> pB ) {
            if( pA->insert_ratio( ) != pB->insert_ratio( ) )
                return pA->insert_ratio( ) < pB->insert_ratio( );
            return pA->query_distance( ) < pB->query_distance( );
        } ); // sort function call

        // i & j are the start and end positions of the clusters, respectiveley
        // this works by setting i to the start of the cluster and then gradually increasing j while it belongs to the
        // complete linkage cluster (we work on sorted linear data)
        size_t uiI = 0;
        size_t uiJ = 0;

        while( uiI < rvEdges.size( ) )
        {
            // increase j if the insert_ratio between the 'i' and 'j' object is closer than 'max_insert_ratio_diff'
            // if we reach a tail edge (those edges are sorted to the end since their insert
            // size is 'inf') we check if the current insert ratio is larger than the tail of the read that created the
            // edge if so we join the tail edge into the cluster.
            if( uiJ < rvEdges.size( ) &&
                rvEdges[ uiI ]->insert_ratio( ) >= ( rvEdges[ uiJ ]->switch_strand_known( )
                                                         ? rvEdges[ uiJ ]->insert_ratio( ) - iMaxInsertRatioDiff
                                                         : (int64_t)rvEdges[ uiJ ]->query_distance( ) ) )
                uiJ++;
            else
            {
                this->exact_sweep( rvEdges, uiI, uiJ, pRet );
                uiI = uiJ;
            } // else
        } // while
    } // method

    virtual std::shared_ptr<CompleteBipartiteSubgraphClusterVector>
        EXPORTED execute( std::shared_ptr<CompleteBipartiteSubgraphClusterVector> pClusters )
    {
        auto pRet = std::make_shared<CompleteBipartiteSubgraphClusterVector>( );

        for( auto pCluster : pClusters->vContent )
            this->lineSweep( pCluster, pRet );

        return pRet;
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
void exportSweepSvJump( py::module& rxPyModuleId );
#endif