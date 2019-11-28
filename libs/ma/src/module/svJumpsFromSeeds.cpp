/**
 * @file svJumsFromSeeds.cpp
 * @author Markus Schmidt
 */
#include "module/svJumpsFromSeeds.h"
#include "pybind11/stl.h"
#include <cmath>
#include <csignal>

using namespace libMA;

Rectangle<nucSeqIndex> SvJumpsFromSeeds::getPositionsForSeeds( Seed& rLast, Seed& rNext, nucSeqIndex uiQSize,
                                                               nucSeqIndex uiRSize )
{
    // seeds are overlapping on query -> rectange size zero
    if( rNext.start( ) < rLast.end( ) )
        return Rectangle<nucSeqIndex>( 0, 0, 0, 0 );

    int64_t iLastRef; // inclusive
    int64_t iNextRef; // inclusive

    if( &rLast == &xDummySeed )
    {
        if( rNext.bOnForwStrand )
            iLastRef = std::max( (int64_t)0, (int64_t)rNext.start_ref( ) -
                                                 ( int64_t )( 1 + rNext.start( ) * dExtraSeedingAreaFactor ) );
        else
            iLastRef =
                std::min( (int64_t)uiRSize,
                          (int64_t)rNext.start_ref( ) + ( int64_t )( rNext.start( ) * dExtraSeedingAreaFactor ) );
    }
    else if( rLast.bOnForwStrand )
        iLastRef = (int64_t)rLast.end_ref( );
    else
        iLastRef = (int64_t)rLast.start_ref( ) - rLast.size( ) + 1;

    if( &rNext == &xDummySeed ) // @todo this is not working...
    {
        if( rLast.bOnForwStrand )
            iNextRef = std::min( (int64_t)uiRSize,
                                 (int64_t)rLast.end_ref( ) +
                                     ( int64_t )( ( uiQSize - rLast.end( ) ) * dExtraSeedingAreaFactor ) );
        else
            iNextRef =
                std::max( (int64_t)0, (int64_t)rLast.start_ref( ) -
                                          ( int64_t )( 1 + ( uiQSize - rLast.end( ) ) * dExtraSeedingAreaFactor ) );
    }
    else if( rNext.bOnForwStrand )
        iNextRef = (int64_t)rNext.start_ref( );
    else
        iNextRef = (int64_t)rNext.start_ref( ) + 1;

    int64_t iRefStart = std::min( iLastRef, iNextRef );
    int64_t iRefEnd = std::max( iLastRef, iNextRef );
    if( iRefStart == iRefEnd )
        return Rectangle<nucSeqIndex>( 0, 0, 0, 0 );
    int64_t iRefSize = iRefEnd - iRefStart;
    return Rectangle<nucSeqIndex>( (nucSeqIndex)iRefStart, rLast.end( ), (nucSeqIndex)iRefSize,
                                   rNext.start( ) - rLast.end( ) );
} // method


void SvJumpsFromSeeds::computeSeeds( Rectangle<nucSeqIndex> xArea, std::shared_ptr<NucSeq> pQuery,
                                     std::shared_ptr<Pack> pRefSeq, std::shared_ptr<Seeds> rvRet )
{
    if( getKMerSizeForRectangle( xArea ) > xArea.xXAxis.size( ) ||
        getKMerSizeForRectangle( xArea ) > xArea.xYAxis.size( ) )
        return;

    HashMapSeeding xHashMapSeeder;
    xHashMapSeeder.uiSeedSize = getKMerSizeForRectangle( xArea );
    // @todo this is inefficient:
    auto pQuerySegment = std::make_shared<NucSeq>( pQuery->fromTo( xArea.xYAxis.start( ), xArea.xYAxis.end( ) ) );
    auto pRef = pRefSeq->vExtract( xArea.xXAxis.start( ), xArea.xXAxis.end( ) );
    auto pSeeds = xHashMapSeeder.execute( pQuerySegment, pRef );

    pRef->vReverseAll( );
    pRef->vSwitchAllBasePairsToComplement( );
    auto pSeedsRev = xHashMapSeeder.execute( pQuerySegment, pRef );
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
} // method

std::shared_ptr<Seeds> SvJumpsFromSeeds::computeSeeds( Rectangle<nucSeqIndex>& xArea, std::shared_ptr<NucSeq> pQuery,
                                                       std::shared_ptr<Pack> pRefSeq )
{
    auto pSeeds = std::make_shared<Seeds>( );
    if( xArea.xXAxis.size( ) < (nucSeqIndex)pSelectedSetting->xMaxSizeReseed->get( ) )
        computeSeeds( xArea, pQuery, pRefSeq, pSeeds );
    else
    {
        return pSeeds;
        // use slightly more than the query size as reference section length
        auto uiRefSecLen = (nucSeqIndex)std::min( xArea.xYAxis.start( ) * dExtraSeedingAreaFactor,
                                                  pSelectedSetting->xMaxSizeReseed->get( ) / 2.0 );
        computeSeeds( Rectangle<nucSeqIndex>(
                          xArea.xXAxis.start( ), xArea.xYAxis.start( ),
                          std::min( uiRefSecLen, pRefSeq->uiUnpackedSizeForwardPlusReverse( ) - xArea.xXAxis.start( ) ),
                          xArea.xYAxis.size( ) ),
                      pQuery, pRefSeq, pSeeds );
        computeSeeds( Rectangle<nucSeqIndex>( xArea.xXAxis.end( ) > uiRefSecLen ? xArea.xXAxis.end( ) - uiRefSecLen : 0,
                                              xArea.xYAxis.start( ),
                                              xArea.xXAxis.end( ) > uiRefSecLen ? uiRefSecLen : xArea.xXAxis.end( ),
                                              xArea.xYAxis.size( ) ),
                      pQuery, pRefSeq, pSeeds );
    } // else

    if( pSeeds->size( ) == 0 )
        return pSeeds;

    // turn k-mers into maximally extended seeds
    auto pvLumpedSeeds = xSeedLumper.execute( pSeeds );
    nucSeqIndex uiMaxSeedSize = 0;
    for( auto& rSeed : *pvLumpedSeeds )
        uiMaxSeedSize = std::max( uiMaxSeedSize, rSeed.size( ) );
    // remove seeds with size less than half than max seed size.
    pvLumpedSeeds->erase( std::remove_if( pvLumpedSeeds->begin( ), pvLumpedSeeds->end( ),
                                          [&]( auto& rSeed ) { return rSeed.size( ) < uiMaxSeedSize; } ),
                          pvLumpedSeeds->end( ) );
    // if( xArea.xXAxis.start( ) == 82431)
    //    std::raise(SIGINT);
    return pvLumpedSeeds;
} // method

/** @brief Determine the appropriate k-mers size for a "rectangle"
 * @detail The formula used over here is:
 * 1 - t <= (1 - 1/4^k)^( (w-k+1)*(h-k+1) )
 * where (1 - 1/4^k) is the probability that two k-sized nucleotide sequences do not match.
 *	     (w-k+1)*(h-k+1) ) is the number of possible K-mer combinations within the rectangle.
 */
nucSeqIndex SvJumpsFromSeeds::getKMerSizeForRectangle( Rectangle<nucSeqIndex>& rRect )
{
    auto w = rRect.xXAxis.size( ); // width of rectangle
    auto h = rRect.xYAxis.size( ); // height of rectangle
    auto& t = SvJumpsFromSeeds::dProbabilityForRandomMatch; // probabilty threshold

    int64_t iDenominator = 4 * 4 * 4;
    for( nucSeqIndex uiK = 3; uiK <= std::min( w, h ); uiK++ )
    {
        // p is the probability of NOT having random match in the rectange (wh, w).
        // So, if p is small, then the probability that we have a random match is high.
        // If p is large, the probabilty of having a random match is low.
        double p = std::pow( 1.0 - ( 1.0 / ( (double)iDenominator ) ), ( w - uiK + 1 ) * ( h - uiK + 1 ) );
        if( p >= 1 - t )
            // The probabilty for a match of length uiK is now below 1 - t.
            return uiK;
        iDenominator *= 4;
    } // for

    // Give up (we cannot reach the required probabilty) and return an impossible k-mer size with respect to matching.
    return std::min( w, h ) + 1;
} // method


void SvJumpsFromSeeds::makeJumpsByReseedingRecursive( Seed& rLast, Seed& rNext, std::shared_ptr<NucSeq> pQuery,
                                                      std::shared_ptr<Pack> pRefSeq,
                                                      std::shared_ptr<ContainerVector<SvJump>>& pRet, size_t uiLayer,
                                                      std::shared_ptr<Seeds> pSeeds,
                                                      std::vector<size_t>* pvLayerOfSeeds,
                                                      std::vector<Rectangle<nucSeqIndex>>* pvRectanglesOut )
{
    // returns a (reference pos, query pos, width [on reference], height [on query])
    auto xRectangle =
        getPositionsForSeeds( rLast, rNext, pQuery->length( ), pRefSeq->uiUnpackedSizeForwardPlusReverse( ) );
    if( pvRectanglesOut != nullptr )
        pvRectanglesOut->push_back( xRectangle );

    // check if there is enough space to reseed
    if( getKMerSizeForRectangle( xRectangle ) <= xRectangle.xXAxis.size( ) &&
        getKMerSizeForRectangle( xRectangle ) <= xRectangle.xYAxis.size( ) )
    {
#if 0
        if( !rLast.bOnForwStrand && !rNext.bOnForwStrand )
            std::cout << "rect: " << xRectangle.xXAxis.start( ) << ", " << xRectangle.xXAxis.end( ) << ", "
                      << xRectangle.xYAxis.start( ) << ", " << xRectangle.xYAxis.end( ) << " rLast: " << rLast.start( )
                      << ", " << rLast.start_ref( ) << ", " << rLast.size( ) << " rNext: " << rNext.start( )
                      << ", " << rNext.start_ref( ) << ", " << rNext.size( ) << std::endl;
#endif
        // reseed and recursiveley call this function again
        std::shared_ptr<Seeds> pvSeeds = computeSeeds( xRectangle, pQuery, pRefSeq );
        // if( xRectangle.xXAxis.start( ) == 82431)
        //    std::raise(SIGINT);
        // if( pvSeeds->size( ) > 0 )
        //{
        //    if( !rLast.bOnForwStrand && !rNext.bOnForwStrand )
        //        std::cout << "#";
        //    else
        //        std::cout << ".";
        //} // if
        std::sort( pvSeeds->begin( ), pvSeeds->end( ),
                   []( const Seed& rA, const Seed& rB ) { return rA.start( ) < rB.start( ); } );
        if( pSeeds != nullptr )
            for( auto& rSeed : *pvSeeds )
            {
                pSeeds->push_back( rSeed );
                pvLayerOfSeeds->push_back( uiLayer );
            } // for

        // this auto locks using the DB; seeds need to be sorted by query position
        xCoverageInserter.insert( *pvSeeds, pQuery->length( ) );

        if( pvSeeds->size( ) > 0 )
        {
            Seed* pCurr = &rLast;
            for( auto& rSeed : *pvSeeds )
            {
                makeJumpsByReseedingRecursive( *pCurr, rSeed, pQuery, pRefSeq, pRet, uiLayer + 1, pSeeds,
                                               pvLayerOfSeeds, pvRectanglesOut );
                pCurr = &rSeed;
            } // for
            makeJumpsByReseedingRecursive( *pCurr, rNext, pQuery, pRefSeq, pRet, uiLayer + 1, pSeeds, pvLayerOfSeeds,
                                           pvRectanglesOut );

            // we found at least one seed so we do not need to create a jump between rLast and rNext
            return;
        } // if
        // else: we did not find any futher seeds in between rLast and rNext -> create the jump between those seeds
    } // if

    // compute a SvJump and terminate the recursion
    if( &rLast == &xDummySeed || &rNext == &xDummySeed )
    {
        // we have to insert a dummy jump if the seed is far enough from the end/start of the query
        if( &rNext != &xDummySeed && rNext.start( ) > uiMinDistDummy )
            pRet->emplace_back( pSelectedSetting, rNext, pQuery->length( ), false, pQuery->iId );
        if( &rLast != &xDummySeed && rLast.end( ) + uiMinDistDummy <= pQuery->length( ) )
            pRet->emplace_back( pSelectedSetting, rLast, pQuery->length( ), true, pQuery->iId );
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

std::shared_ptr<ContainerVector<SvJump>>
SvJumpsFromSeeds::execute_helper( std::shared_ptr<SegmentVector> pSegments,
                                  std::shared_ptr<Pack>
                                      pRefSeq,
                                  std::shared_ptr<FMIndex>
                                      pFM_index,
                                  std::shared_ptr<NucSeq>
                                      pQuery,
                                  std::shared_ptr<Seeds>
                                      pSeeds,
                                  std::vector<size_t>* pvLayerOfSeeds,
                                  std::vector<Rectangle<nucSeqIndex>>* pvRectanglesOut )
{
    AlignedMemoryManager xMemoryManager;
    auto pRet = std::make_shared<ContainerVector<SvJump>>( );

    // sort seeds by query position:
    std::sort( pSegments->begin( ), pSegments->end( ),
               []( const Segment& rA, const Segment& rB ) { return rA.start( ) < rB.start( ); } );
    Seeds vSeeds;
    // avoid multiple allocations (we can only guess the actual number of seeds here)
    vSeeds.reserve( pSegments->size( ) * 2 );

    // filter ambiguous segments -> 1; don't -> 0
#if 1
    /**
     * This filters segments that occur multiple times on the reference:
     * For each of those segments extract all seeds.
     * Then pick the one seed that is closest, on the reference, to:
     *  - either the last (on query) unique seed
     *  - or the next (on query) unique seed
     * This drastically reduces the number of seeds.
     * Further, it shall help dealing with ambiguous regions by eliminating all but the seeds most likely to fit
     * into a chain of seeds
     */
    std::vector<Segment*> vTemp;
    size_t uiNumSeedsTotal = 0;
    int64_t iLastUniqueRefPos = -1;
    for( Segment& rSegment : *pSegments )
    {
        if( rSegment.size( ) < uiMinSeedSizeSV )
            continue;
        uiNumSeedsTotal += rSegment.saInterval( ).size( );

        if( rSegment.saInterval( ).size( ) == 1 )
            rSegment.forEachSeed( *pFM_index, [&]( Seed& rS ) {
                // deal with vTemp
                int64_t iCurrUniqueRefPos = (int64_t)rS.start_ref( );
                for( Segment* pSeg : vTemp )
                {
                    Seed xBest;
                    int64_t iMinDist = -1;
                    pSeg->forEachSeed( *pFM_index, [&]( Seed& rS ) {
                        int64_t iStart = (int64_t)rS.start_ref( );
                        int64_t iDist = std::abs( iStart - iCurrUniqueRefPos );
                        // if the segment is not the first segment this triggers
                        if( iLastUniqueRefPos != -1 )
                            iDist = std::min( iDist, std::abs( iStart - iLastUniqueRefPos ) );
                        if( iMinDist == -1 || iDist < iMinDist )
                        {
                            xBest = rS;
                            iMinDist = iDist;
                        } // if
                        return true;
                    } ); // forEachSeed function call
                    assert( iMinDist != -1 );
                    vSeeds.push_back( xBest );
                } // for
                vTemp.clear( );

                iLastUniqueRefPos = (int64_t)rS.end_ref( );
                vSeeds.push_back( rS );
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
            int64_t iStart = (int64_t)rS.start_ref( );
            int64_t iDist = std::abs( iStart - iLastUniqueRefPos );
            if( iMinDist == -1 || iDist < iMinDist )
            {
                xBest = rS;
                iMinDist = iDist;
            } // if
            return true;
        } ); // forEachSeed function call
        assert( iMinDist != -1 );
        vSeeds.push_back( xBest );
    } // for

    {
        std::lock_guard<std::mutex> xGuard( xLock );
        uiNumSeedsEliminatedAmbiguityFilter += uiNumSeedsTotal - vSeeds.size( );
        uiNumSeedsKeptAmbiguityFilter += vSeeds.size( );
    } // scope xGuard

#else
    pSegments->emplaceAllEachSeeds( *pFM_index, pQuery->length( ), uiMaxAmbiguitySv, uiMinSeedSizeSV, vSeeds,
                                    [&]( ) { return true; } );
#endif

    if( pSeeds != nullptr )
        for( auto& rSeed : vSeeds )
        {
            pSeeds->push_back( rSeed );
            pvLayerOfSeeds->push_back( 0 );
        } // for

    // insert coverage
    xCoverageInserter.insert( vSeeds, pQuery->length( ) );

    // actually compute the jumps
    Seed* pCurr = &xDummySeed;
    for( auto& rSeed : vSeeds )
    {
        makeJumpsByReseedingRecursive( *pCurr, rSeed, pQuery, pRefSeq, pRet, 1, pSeeds, pvLayerOfSeeds,
                                       pvRectanglesOut );
        pCurr = &rSeed;
    } // for
    makeJumpsByReseedingRecursive( *pCurr, xDummySeed, pQuery, pRefSeq, pRet, 1, pSeeds, pvLayerOfSeeds,
                                   pvRectanglesOut );

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
    return execute_helper( pSegments, pRefSeq, pFM_index, pQuery, nullptr, nullptr, nullptr );
} // method

#ifdef WITH_PYTHON
void exportSvJumpsFromSeeds( py::module& rxPyModuleId )
{
    py::class_<Interval<nucSeqIndex>>( rxPyModuleId, "nucSeqInterval" )
        .def_readwrite( "start", &Interval<nucSeqIndex>::iStart )
        .def_readwrite( "size", &Interval<nucSeqIndex>::iSize );
    py::class_<Rectangle<nucSeqIndex>>( rxPyModuleId, "nucSeqRectangle" )
        .def_readwrite( "x_axis", &Rectangle<nucSeqIndex>::xXAxis )
        .def_readwrite( "y_axis", &Rectangle<nucSeqIndex>::xYAxis );
    exportModule<SvJumpsFromSeeds, int64_t, std::shared_ptr<SV_DB>, std::shared_ptr<Pack>>(
        rxPyModuleId, "SvJumpsFromSeeds", []( auto&& x ) {
            x.def( "commit", &SvJumpsFromSeeds::commit ) //
                .def( "execute_helper", &SvJumpsFromSeeds::execute_helper_py );
        } );
} // function
#endif