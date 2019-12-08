/**
 * @file svJumpsFromSeeds.cpp
 * @author Markus Schmidt
 */
#include "module/svJumpsFromSeeds.h"
#include "module/module.h"
#include <cmath>
#include <csignal>

using namespace libMA;

// @todo this whole function can be simplified a lot
std::pair<libMA::Rectangle<nucSeqIndex>, libMA::Rectangle<nucSeqIndex>>
SvJumpsFromSeeds::getPositionsForSeeds( Seed& rLast, Seed& rNext, nucSeqIndex uiQStart, nucSeqIndex uiQEnd,
                                        std::shared_ptr<Pack> pRefSeq )
{
    // seeds are overlapping on query -> rectange size zero
    if( &rLast != &xDummySeed && &rNext != &xDummySeed && rNext.start( ) < rLast.end( ) )
        return std::make_pair( libMA::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ),
                               libMA::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );
    if( &rLast != &xDummySeed && rLast.end( ) >= uiQEnd )
        return std::make_pair( libMA::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ),
                               libMA::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );
    if( &rNext != &xDummySeed && rNext.start( ) <= uiQStart )
        return std::make_pair( libMA::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ),
                               libMA::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );

    int64_t iLastRef; // inclusive
    int64_t iNextRef; // exclusive

    if( &rLast == &xDummySeed )
    {
        if( rNext.bOnForwStrand )
        {
            // make sure we do not create a contig bridging rectangle
            int64_t iStartOfContig =
                pRefSeq->startOfSequenceWithId( pRefSeq->uiSequenceIdForPosition( rNext.start_ref( ) ) ); // inclusive
            // estimate size fo rectangle based on query position of existing seed
            int64_t iQueryDistance = ( int64_t )( ( rNext.start( ) - uiQStart ) * dExtraSeedingAreaFactor );
            if( iQueryDistance > (int64_t)pSelectedSetting->xMaxSizeReseed->get( ) / 2 )
                iQueryDistance = (int64_t)pSelectedSetting->xMaxSizeReseed->get( ) / 2;
            iLastRef = std::max( iStartOfContig, (int64_t)rNext.start_ref( ) - iQueryDistance );
        } // if
        else
        {
            // make sure we do not create a contig bridging rectangle
            int64_t iEndOfContig =
                pRefSeq->endOfSequenceWithId( pRefSeq->uiSequenceIdForPosition( rNext.start_ref( ) + 1 ) ); // exclusive
            // estimate size fo rectangle based on query position of existing seed
            int64_t iQueryDistance = ( int64_t )( ( rNext.start( ) - uiQStart ) * dExtraSeedingAreaFactor );
            if( iQueryDistance > (int64_t)pSelectedSetting->xMaxSizeReseed->get( ) / 2 )
                iQueryDistance = (int64_t)pSelectedSetting->xMaxSizeReseed->get( ) / 2;
            iLastRef = std::min( iEndOfContig, (int64_t)rNext.start_ref( ) + 1 + iQueryDistance );
        } // else
    } // if
    else if( rLast.bOnForwStrand )
        iLastRef = (int64_t)rLast.end_ref( );
    else
        iLastRef = (int64_t)rLast.start_ref( ) - rLast.size( ) + 1;

    if( &rNext == &xDummySeed )
    {
        if( rLast.bOnForwStrand )
        {
            // make sure we do not create a contig bridging rectangle
            int64_t iEndOfContig =
                pRefSeq->endOfSequenceWithId( pRefSeq->uiSequenceIdForPosition( rLast.end_ref( ) ) ); // exclusive
            // estimate size fo rectangle based on query position of existing seed
            int64_t iQueryDistance = ( int64_t )( ( uiQEnd - rLast.end( ) ) * dExtraSeedingAreaFactor );
            if( iQueryDistance > (int64_t)pSelectedSetting->xMaxSizeReseed->get( ) / 2 )
                iQueryDistance = (int64_t)pSelectedSetting->xMaxSizeReseed->get( ) / 2;
            iNextRef = std::min( iEndOfContig, (int64_t)rLast.end_ref( ) + iQueryDistance );
        } // if
        else
        {
            // make sure we do not create a contig bridging rectangle
            int64_t iStartOfContig = pRefSeq->startOfSequenceWithId(
                pRefSeq->uiSequenceIdForPosition( rLast.start_ref( ) + 1 - rLast.size( ) ) ); // inclusive
            // estimate size fo rectangle based on query position of existing seed
            int64_t iQueryDistance = ( int64_t )( ( uiQEnd - rLast.end( ) ) * dExtraSeedingAreaFactor );
            if( iQueryDistance > (int64_t)pSelectedSetting->xMaxSizeReseed->get( ) / 2 )
                iQueryDistance = (int64_t)pSelectedSetting->xMaxSizeReseed->get( ) / 2;
            iNextRef = std::max( iStartOfContig,
                                 (int64_t)rLast.start_ref( ) + 1 - ( int64_t )( rLast.size( ) + iQueryDistance ) );
        } // else
    } // if
    else if( rNext.bOnForwStrand )
        iNextRef = (int64_t)rNext.start_ref( );
    else
        iNextRef = (int64_t)rNext.start_ref( ) + 1;

    if( iLastRef == iNextRef )
        return std::make_pair( libMA::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ),
                               libMA::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );
    if( &rLast != &xDummySeed && &rNext != &xDummySeed )
    {
        int64_t iRefSize;
        if( rLast.bOnForwStrand && rNext.bOnForwStrand )
            iRefSize = iNextRef - iLastRef;
        else if( !rLast.bOnForwStrand && !rNext.bOnForwStrand )
            iRefSize = iLastRef - iNextRef;
        else
            iRefSize = -1; // seeds on different strands always need separate rectangles.

        if( iRefSize > (int64_t)pSelectedSetting->xMaxSizeReseed->get( ) || // check for too large rectangle
            iRefSize < 0 || // check for seeds pointing the wrong direction
            // check for rectangle bridging multiple contigs
            pRefSeq->uiSequenceIdForPosition( iLastRef ) != pRefSeq->uiSequenceIdForPosition( iNextRef - 1 ) )
            // if the rectangle between the seeds is too large or spans over tow contigs, split it in two
            return std::make_pair(
                getPositionsForSeeds( rLast, xDummySeed, rLast.end( ), rNext.start( ), pRefSeq ).first,
                getPositionsForSeeds( xDummySeed, rNext, rLast.end( ), rNext.start( ), pRefSeq ).first );
    }
    int64_t iRefStart = std::min( iLastRef, iNextRef );
    int64_t iRefEnd = std::max( iLastRef, iNextRef );
    int64_t iRefSize = iRefEnd - iRefStart;

    int64_t iRectQStart = &rLast != &xDummySeed ? rLast.end( ) : uiQStart;
    int64_t iRectQEnd = &rNext != &xDummySeed ? rNext.start( ) : uiQEnd;
    return std::make_pair( libMA::Rectangle<nucSeqIndex>( (nucSeqIndex)iRefStart, (nucSeqIndex)iRectQStart,
                                                          (nucSeqIndex)iRefSize,
                                                          ( nucSeqIndex )( iRectQEnd - iRectQStart ) ),
                           libMA::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );
} // method


void SvJumpsFromSeeds::computeSeeds( libMA::Rectangle<nucSeqIndex>& xArea, std::shared_ptr<NucSeq> pQuery,
                                     std::shared_ptr<Pack> pRefSeq, std::shared_ptr<Seeds> rvRet,
                                     HelperRetVal* pOutExtra )
{
    if( xArea.xXAxis.size( ) == 0 || xArea.xYAxis.size( ) == 0 )
    {
        if( pOutExtra != nullptr )
            pOutExtra->vRectangleReferenceAmbiguity.push_back( 0 );
        return;
    } // if
    auto pRef = pRefSeq->vExtract( xArea.xXAxis.start( ), xArea.xXAxis.end( ) );
    auto pRefRevComp = pRefSeq->vExtract( xArea.xXAxis.start( ), xArea.xXAxis.end( ) );
    pRefRevComp->vReverseAll( );
    pRefRevComp->vSwitchAllBasePairsToComplement( );
    /*
     * In order to kill parlindromes:
     * sample the k-mer size of forward and reverse strand together
     */
    auto uiSampledKMerSize = sampleKMerSize( *pRef, *pRefRevComp, this->dProbabilityForRandomMatch );
    if( pOutExtra != nullptr )
        pOutExtra->vRectangleReferenceAmbiguity.push_back( uiSampledKMerSize );
    if( uiSampledKMerSize <= uiMaxAddSeedSize )
    {
        HashMapSeeding xHashMapSeeder;
        xHashMapSeeder.uiSeedSize = getKMerSizeForRectangle( xArea, this->dProbabilityForRandomMatch );
        if( xHashMapSeeder.uiSeedSize > xArea.xXAxis.size( ) || xHashMapSeeder.uiSeedSize > xArea.xYAxis.size( ) )
            return;

        // @todo this is inefficient:
        auto pQuerySegment = std::make_shared<NucSeq>( pQuery->fromTo( xArea.xYAxis.start( ), xArea.xYAxis.end( ) ) );
        auto pSeeds = xHashMapSeeder.execute( pQuerySegment, pRef );

        auto pSeedsRev = xHashMapSeeder.execute( pQuerySegment, pRefRevComp );
        // fix seed positions on forward strand
        for( Seed& rSeed : *pSeeds )
        {
            rSeed.uiPosOnReference += xArea.xXAxis.start( );
            rSeed.iStart += xArea.xYAxis.start( );
            assert( rSeed.end( ) <= pQuery->length( ) );
        } // for
        // fix seed positions on reverse strand
        for( Seed& rSeed : *pSeedsRev )
        {
            rSeed.bOnForwStrand = false;
            // undo reversion of reference
            assert( xArea.xXAxis.size( ) >= rSeed.uiPosOnReference + 1 );
            // rSeed.uiPosOnReference = xArea.xXAxis.size( ) - rSeed.uiPosOnReference - 1;
            // rSeed.uiPosOnReference += xArea.xXAxis.start( );
            assert( xArea.xXAxis.end( ) - rSeed.uiPosOnReference >= 1 );
            rSeed.uiPosOnReference = xArea.xXAxis.end( ) - rSeed.uiPosOnReference - 1;
            rSeed.iStart += xArea.xYAxis.start( );
            assert( rSeed.end( ) <= pQuery->length( ) );
        } // for

        pSeeds->confirmSeedPositions( pQuery, pRefSeq, false );
        pSeedsRev->confirmSeedPositions( pQuery, pRefSeq, false );

        rvRet->append( pSeeds );
        rvRet->append( pSeedsRev );
    }
    else
    {
        auto pFAlignment = std::make_shared<Alignment>( xArea.xXAxis.start( ), xArea.xYAxis.start( ) );
        AlignedMemoryManager xMemoryManager; // @todo this causes frequent (de-)&allocations; move this outwards
        xNW.ksw( pQuery, pRef, xArea.xYAxis.start( ), xArea.xYAxis.end( ) - 1, 0, pRef->length( ) - 1, pFAlignment,
                 xMemoryManager );
        auto pForwSeeds = pFAlignment->toSeeds( );

        // and now the reverse strand seeds
        auto pRAlignment = std::make_shared<Alignment>( );
        xNW.ksw( pQuery, pRefRevComp, xArea.xYAxis.start( ), xArea.xYAxis.end( ) - 1, 0, pRefRevComp->length( ) - 1,
                 pRAlignment, xMemoryManager );
        auto pRevSeeds = pRAlignment->toSeeds( );
        for( Seed& rSeed : *pRevSeeds )
        {
            rSeed.bOnForwStrand = false;
            // undo reversion of reference
            assert( xArea.xXAxis.size( ) >= rSeed.uiPosOnReference + 1 );
            // rSeed.uiPosOnReference = xArea.xXAxis.size( ) - rSeed.uiPosOnReference - 1;
            // rSeed.uiPosOnReference += xArea.xXAxis.start( );
            assert( xArea.xXAxis.end( ) - rSeed.uiPosOnReference >= 1 );
            rSeed.uiPosOnReference = xArea.xXAxis.end( ) - rSeed.uiPosOnReference - 1;
            rSeed.iStart += xArea.xYAxis.start( );
            assert( rSeed.end( ) <= pQuery->length( ) );
        } // for
        if( pFAlignment->score() >= pRAlignment->score() )
            rvRet->append( pForwSeeds );
        else
            rvRet->append( pRevSeeds );
    } // else
} // method

std::shared_ptr<Seeds>
SvJumpsFromSeeds::computeSeeds( std::pair<libMA::Rectangle<nucSeqIndex>, libMA::Rectangle<nucSeqIndex>>& xAreas,
                                std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRefSeq, HelperRetVal* pOutExtra )
{
    auto pSeeds = std::make_shared<Seeds>( );
    computeSeeds( xAreas.first, pQuery, pRefSeq, pSeeds, pOutExtra );
    computeSeeds( xAreas.second, pQuery, pRefSeq, pSeeds, pOutExtra );

    if( pSeeds->size( ) == 0 )
        return pSeeds;

    // turn k-mers into maximally extended seeds
    auto pvLumpedSeeds = xSeedLumper.execute( pSeeds );
    // seeds that touch end of rectangle might not be maximally extended
    return SeedExtender( ).execute( pvLumpedSeeds, pQuery, pRefSeq );
} // method


void SvJumpsFromSeeds::makeJumpsByReseedingRecursive( Seed& rLast, Seed& rNext, std::shared_ptr<NucSeq> pQuery,
                                                      std::shared_ptr<Pack> pRefSeq,
                                                      std::shared_ptr<ContainerVector<SvJump>>& pRet, size_t uiLayer,
                                                      HelperRetVal* pOutExtra )
{
    // returns a (reference pos, query pos, width [on reference], height [on query])
    auto xRectangles = getPositionsForSeeds( rLast, rNext, 0, pQuery->length( ), pRefSeq );
    if( pOutExtra != nullptr )
    {
        pOutExtra->vRectangles.push_back( xRectangles.first );
        pOutExtra->vRectangles.push_back( xRectangles.second );
        xParlindromeFilter.keepParlindromes( );
    } // if

    // reseed and recursiveley call this function again
    std::shared_ptr<Seeds> pvSeeds =
        xParlindromeFilter.execute( computeSeeds( xRectangles, pQuery, pRefSeq, pOutExtra ) );
    if( pOutExtra != nullptr )
    {
        pOutExtra->vRectangleFillPercentage.push_back( rectFillPercentage( pvSeeds, xRectangles ) );
        pOutExtra->vRectangleFillPercentage.push_back( rectFillPercentage( pvSeeds, xRectangles ) );
    } // if

    std::sort( pvSeeds->begin( ), pvSeeds->end( ),
               []( const Seed& rA, const Seed& rB ) { return rA.start( ) < rB.start( ); } );
    if( pOutExtra != nullptr )
    {
        for( auto& rSeed : *pvSeeds )
        {
            pOutExtra->pSeeds->push_back( rSeed );
            pOutExtra->vLayerOfSeeds.push_back( uiLayer );
            pOutExtra->vParlindromeSeed.push_back( false );
        } // for
        for( auto& rSeed : *xParlindromeFilter.pParlindromes )
        {
            pOutExtra->pSeeds->push_back( rSeed );
            pOutExtra->vLayerOfSeeds.push_back( uiLayer );
            pOutExtra->vParlindromeSeed.push_back( true );
        } // for
    } // if

    // this auto locks using the DB; seeds need to be sorted by query position
    xCoverageInserter.insert( *pvSeeds, pQuery->length( ) );

    if( pvSeeds->size( ) > 0 )
    {
        Seed* pCurr = &rLast;
        for( auto& rSeed : *pvSeeds )
        {
            makeJumpsByReseedingRecursive( *pCurr, rSeed, pQuery, pRefSeq, pRet, uiLayer + 1, pOutExtra );
            pCurr = &rSeed;
        } // for
        makeJumpsByReseedingRecursive( *pCurr, rNext, pQuery, pRefSeq, pRet, uiLayer + 1, pOutExtra );

        // we found at least one seed so we do not need to create a jump between rLast and rNext
        return;
    } // if

    // compute a SvJump and terminate the recursion
    if( ( &rLast == &xDummySeed || &rNext == &xDummySeed ) && bDoDummyJumps )
    {
        // we have to insert a dummy jump if the seed is far enough from the end/start of the query
        if( &rNext != &xDummySeed && rNext.start( ) > uiMinDistDummy )
            pRet->emplace_back( pSelectedSetting, rNext, pQuery->length( ), false, pQuery->iId, uiMaxDistDummy );
        if( &rLast != &xDummySeed && rLast.end( ) + uiMinDistDummy <= pQuery->length( ) )
            pRet->emplace_back( pSelectedSetting, rLast, pQuery->length( ), true, pQuery->iId, uiMaxDistDummy );
    } // if
    else
    {
        // we have to insert a jump between two seeds
        if( SvJump::validJump( rNext, rLast, true ) )
            pRet->emplace_back( pSelectedSetting, rNext, rLast, true, pQuery->iId );
        if( SvJump::validJump( rLast, rNext, false ) )
            pRet->emplace_back( pSelectedSetting, rLast, rNext, false, pQuery->iId );
    } // else
} // function

std::shared_ptr<ContainerVector<SvJump>> SvJumpsFromSeeds::execute_helper( std::shared_ptr<SegmentVector> pSegments,
                                                                           std::shared_ptr<Pack>
                                                                               pRefSeq,
                                                                           std::shared_ptr<FMIndex>
                                                                               pFM_index,
                                                                           std::shared_ptr<NucSeq>
                                                                               pQuery,
                                                                           HelperRetVal* pOutExtra )
{
    AlignedMemoryManager xMemoryManager;
    auto pRet = std::make_shared<ContainerVector<SvJump>>( );

    // sort seeds by query position:
    std::sort( pSegments->begin( ), pSegments->end( ),
               []( const Segment& rA, const Segment& rB ) { return rA.start( ) < rB.start( ); } );
    auto pSeeds = std::make_shared<Seeds>( );
    // avoid multiple allocations (we can only guess the actual number of seeds here)
    pSeeds->reserve( pSegments->size( ) * 2 );

    // filter ambiguous segments -> 1; don't -> 0
#if 0
    /**
     * This filters segments that occur multiple times on the reference:
     * For each of those segments extract all seeds.
     * Then pick the one seed that is closest, via their delta pos, to:
     *  - either the last (on query) unique seed
     *  - or the next (on query) unique seed
     * This drastically reduces the number of seeds.
     * Further, it shall help dealing with ambiguous regions by eliminating all but the seeds most likely to fit
     * into a chain of seeds
     */
    std::vector<Segment*> vTemp;
    size_t uiNumSeedsTotal = 0;
    int64_t iLastUniqueDelta = -1;
    for( Segment& rSegment : *pSegments )
    {
        if( rSegment.size( ) < uiMinSeedSizeSV )
            continue;
        uiNumSeedsTotal += rSegment.saInterval( ).size( );

        if( rSegment.saInterval( ).size( ) == 1 )
            rSegment.forEachSeed( *pFM_index, [&]( Seed& rS ) {
                // deal with vTemp
                int64_t iCurrUniqueDelta = (int64_t)rS.start( ) - (int64_t)rS.start_ref( );
                for( Segment* pSeg : vTemp )
                {
                    Seed xBest;
                    int64_t iMinDist = -1;
                    pSeg->forEachSeed( *pFM_index, [&]( Seed& rS ) {
                        // @todo this should be useing the delta distances...
                        int64_t iDelta = (int64_t)rS.start( ) - (int64_t)rS.start_ref( );
                        int64_t iDist = std::abs( iDelta - iCurrUniqueDelta );
                        // if the segment is not the first segment this triggers
                        if( iLastUniqueDelta != -1 )
                            iDist = std::min( iDist, std::abs( iDelta - iLastUniqueDelta ) );
                        if( iMinDist == -1 || iDist < iMinDist )
                        {
                            xBest = rS;
                            iMinDist = iDist;
                        } // if
                        return true;
                    } ); // forEachSeed function call
                    assert( iMinDist != -1 );
                    pSeeds->push_back( xBest );
                } // for
                vTemp.clear( );

                iLastUniqueDelta = (int64_t)rS.start( ) - (int64_t)rS.start_ref( );
                pSeeds->push_back( rS );
                return true;
            } ); // forEachSeed function call \ if
        else
            vTemp.push_back( &rSegment );
    } // for
    // seeds at the end
    for( Segment* pSeg : vTemp )
    {
        Seed xBest;
        int64_t iMinDist = -1;
        pSeg->forEachSeed( *pFM_index, [&]( Seed& rS ) {
            int64_t iDelta = (int64_t)rS.start( ) - (int64_t)rS.start_ref( );
            int64_t iDist = std::abs( iDelta - iLastUniqueDelta );
            if( iMinDist == -1 || iDist < iMinDist )
            {
                xBest = rS;
                iMinDist = iDist;
            } // if
            return true;
        } ); // forEachSeed function call
        assert( iMinDist != -1 );
        pSeeds->push_back( xBest );
    } // for

    {
        std::lock_guard<std::mutex> xGuard( xLock );
        uiNumSeedsEliminatedAmbiguityFilter += uiNumSeedsTotal - pSeeds->size( );
        uiNumSeedsKeptAmbiguityFilter += pSeeds->size( );
    } // scope xGuard

#else
    pSegments->emplaceAllEachSeeds( *pFM_index, pQuery->length( ), /*uiMaxAmbiguitySv*/ 1, uiMinSeedSizeSV, *pSeeds,
                                    [&]( ) { return true; } );
#endif

    if( pOutExtra != nullptr )
        xParlindromeFilter.keepParlindromes( );
    auto pFilteredSeeds = xParlindromeFilter.execute( pSeeds );
    std::sort( pFilteredSeeds->begin( ), pFilteredSeeds->end( ),
               []( const Seed& rA, const Seed& rB ) { return rA.start( ) < rB.start( ); } );

    if( pOutExtra != nullptr )
    {
        for( auto& rSeed : *pFilteredSeeds )
        {
            pOutExtra->pSeeds->push_back( rSeed );
            pOutExtra->vLayerOfSeeds.push_back( 0 );
            pOutExtra->vParlindromeSeed.push_back( false );
        } // for
        for( auto& rSeed : *xParlindromeFilter.pParlindromes )
        {
            pOutExtra->pSeeds->push_back( rSeed );
            pOutExtra->vLayerOfSeeds.push_back( 0 );
            pOutExtra->vParlindromeSeed.push_back( true );
        } // for
    } // if

    // insert coverage
    xCoverageInserter.insert( *pFilteredSeeds, pQuery->length( ) );

    // actually compute the jumps
    Seed* pCurr = &xDummySeed;
    for( auto& rSeed : *pFilteredSeeds )
    {
        makeJumpsByReseedingRecursive( *pCurr, rSeed, pQuery, pRefSeq, pRet, 1, pOutExtra );
        pCurr = &rSeed;
    } // for
    makeJumpsByReseedingRecursive( *pCurr, xDummySeed, pQuery, pRefSeq, pRet, 1, pOutExtra );

    return pRet;
} // method

std::shared_ptr<ContainerVector<SvJump>> SvJumpsFromSeeds::execute( std::shared_ptr<SegmentVector> pSegments,
                                                                    std::shared_ptr<Pack>
                                                                        pRefSeq,
                                                                    std::shared_ptr<FMIndex>
                                                                        pFM_index,
                                                                    std::shared_ptr<NucSeq>
                                                                        pQuery )
{
    return execute_helper( pSegments, pRefSeq, pFM_index, pQuery, nullptr );
} // method

#ifdef WITH_PYTHON
void exportSvJumpsFromSeeds( py::module& rxPyModuleId )
{
    py::class_<Interval<nucSeqIndex>>( rxPyModuleId, "nucSeqInterval" )
        .def_readwrite( "start", &Interval<nucSeqIndex>::iStart )
        .def_readwrite( "size", &Interval<nucSeqIndex>::iSize );
    py::class_<libMA::Rectangle<nucSeqIndex>>( rxPyModuleId, "nucSeqRectangle" )
        .def_readwrite( "x_axis", &libMA::Rectangle<nucSeqIndex>::xXAxis )
        .def_readwrite( "y_axis", &libMA::Rectangle<nucSeqIndex>::xYAxis );
    py::class_<libMA::SvJumpsFromSeeds::HelperRetVal>( rxPyModuleId, "SvJumpsFromSeedsHelperRetVal" )
        .def_readwrite( "layer_of_seeds", &libMA::SvJumpsFromSeeds::HelperRetVal::vLayerOfSeeds )
        .def_readwrite( "seeds", &libMA::SvJumpsFromSeeds::HelperRetVal::pSeeds )
        .def_readwrite( "rectangles", &libMA::SvJumpsFromSeeds::HelperRetVal::vRectangles )
        .def_readwrite( "parlindrome", &libMA::SvJumpsFromSeeds::HelperRetVal::vParlindromeSeed )
        .def_readwrite( "rectangles_fill", &libMA::SvJumpsFromSeeds::HelperRetVal::vRectangleFillPercentage )
        .def_readwrite( "rectangle_ambiguity", &libMA::SvJumpsFromSeeds::HelperRetVal::vRectangleReferenceAmbiguity );
    py::bind_vector<std::vector<libMA::Rectangle<nucSeqIndex>>>( rxPyModuleId, "RectangleVector", "" );
    py::bind_vector<std::vector<double>>( rxPyModuleId, "DoubleVector", "" );
    py::bind_vector<std::vector<bool>>( rxPyModuleId, "BoolVector", "" );
    exportModule<SvJumpsFromSeeds, int64_t, std::shared_ptr<SV_DB>, std::shared_ptr<Pack>>(
        rxPyModuleId, "SvJumpsFromSeeds", []( auto&& x ) {
            x.def( "commit", &SvJumpsFromSeeds::commit ) //
                .def( "execute_helper", &SvJumpsFromSeeds::execute_helper_py );
        } );
} // function
#endif