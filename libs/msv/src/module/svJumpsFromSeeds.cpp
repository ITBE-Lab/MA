/**
 * @file svJumpsFromSeeds.cpp
 * @author Markus Schmidt
 */
#include "msv/module/svJumpsFromSeeds.h"
#include "ms/module/module.h"
#include "ms/util/pybind11.h"
#include "pybind11/stl.h"
#include <cmath>
#include <csignal>

using namespace libMSV;
using namespace libMS;

// @todo this whole function can be simplified a lot
std::pair<geom::Rectangle<nucSeqIndex>, geom::Rectangle<nucSeqIndex>>
SvJumpsFromSeeds::getPositionsForSeeds( Seed& rLast, Seed& rNext, nucSeqIndex uiQStart, nucSeqIndex uiQEnd,
                                        std::shared_ptr<Pack> pRefSeq )
{
    if( &rNext == &xDummySeed )
    {
        nucSeqIndex uiBottom = rLast.end( );

        if( uiBottom >= uiQEnd )
            return std::make_pair( geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ),
                                   geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );

        nucSeqIndex uiTop = uiQEnd;
        // limit height of rectangle
        if( (int64_t)uiTop - (int64_t)uiBottom > iMaxSizeReseed / 2 )
            uiTop = uiBottom + iMaxSizeReseed / 2;

        int64_t iLeft;
        int64_t iRight;
        if( rLast.bOnForwStrand )
        {
            iLeft = ( int64_t )( rLast.end_ref( ) );
            // exclusive
            int64_t iEndOfContig = pRefSeq->endOfSequenceWithId( pRefSeq->uiSequenceIdForPosition( iLeft ) );
            // make sure we do not create a contig bridging rectangle
            iRight = std::min( iLeft + iMaxSizeReseed / 2, iEndOfContig );

            // make square
            if( iRight - iLeft < (int64_t)uiTop - (int64_t)uiBottom )
                uiTop = uiBottom + ( iRight - iLeft );
            else
                iRight = iLeft + ( uiTop - uiBottom );
        } // else
        else
        {
            iRight = ( int64_t )( rLast.start_ref( ) - rLast.size( ) + 1 );
            // inclusive
            int64_t iStartOfContig = pRefSeq->startOfSequenceWithId( pRefSeq->uiSequenceIdForPosition( iRight ) );
            // make sure we do not create a contig bridging rectangle
            iLeft = std::max( iRight - iMaxSizeReseed / 2, iStartOfContig );

            // make square
            if( iRight - iLeft < (int64_t)uiTop - (int64_t)uiBottom )
                uiTop = uiBottom + ( iRight - iLeft );
            else
                iLeft = iRight - ( uiTop - uiBottom );
        } // else
        // can return as a single rectangle
        return std::make_pair( geom::Rectangle<nucSeqIndex>( (nucSeqIndex)iLeft, uiBottom,
                                                             ( nucSeqIndex )( iRight - iLeft ), uiTop - uiBottom ),
                               geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );
    } // if
    else if( &rLast == &xDummySeed )
    {
        nucSeqIndex uiTop = rNext.start( );

        if( uiTop <= uiQStart )
            return std::make_pair( geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ),
                                   geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );

        nucSeqIndex uiBottom = uiQStart;
        // limit height of rectangle
        if( (int64_t)uiTop - (int64_t)uiBottom > iMaxSizeReseed / 2 )
        {
            if( (int64_t)uiTop < iMaxSizeReseed / 2 )
                uiBottom = 0;
            else
                uiBottom = uiTop - iMaxSizeReseed / 2;
        } // if

        int64_t iLeft;
        int64_t iRight;
        if( rNext.bOnForwStrand )
        {
            iRight = ( int64_t )( rNext.start_ref( ) );
            // inclusive
            int64_t iStartOfContig = pRefSeq->startOfSequenceWithId( pRefSeq->uiSequenceIdForPosition( iRight ) );
            // make sure we do not create a contig bridging rectangle
            iLeft = std::max( iRight - iMaxSizeReseed / 2, iStartOfContig );

            // make square
            if( iRight - iLeft < (int64_t)uiTop - (int64_t)uiBottom )
                uiBottom = uiTop - ( iRight - iLeft );
            else
                iLeft = iRight - ( uiTop - uiBottom );
        } // else
        else
        {
            iLeft = ( int64_t )( rNext.start_ref( ) + 1 );
            // exclusive
            int64_t iEndOfContig = pRefSeq->endOfSequenceWithId( pRefSeq->uiSequenceIdForPosition( iLeft ) );
            // make sure we do not create a contig bridging rectangle
            iRight = std::min( iLeft + iMaxSizeReseed / 2, iEndOfContig );

            // make square
            if( iRight - iLeft < (int64_t)uiTop - (int64_t)uiBottom )
                uiBottom = uiTop - ( iRight - iLeft );
            else
                iRight = iLeft + ( uiTop - uiBottom );
        } // else
        // can return as a single rectangle
        return std::make_pair( geom::Rectangle<nucSeqIndex>( (nucSeqIndex)iLeft, uiBottom,
                                                             ( nucSeqIndex )( iRight - iLeft ), uiTop - uiBottom ),
                               geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );
    } // else if


    // if seeds are on same strand we have a chance to connect the rectangles
    else if( rLast.bOnForwStrand == rNext.bOnForwStrand )
    {
        nucSeqIndex uiBottom = rLast.end( );
        nucSeqIndex uiTop = rNext.start( );

        // if( uiTop <= uiBottom )
        //    return std::make_pair( geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ),
        //                           geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );

        nucSeqIndex uiLeft;
        if( rLast.bOnForwStrand )
            uiLeft = rLast.end_ref( );
        else
            uiLeft = rNext.start_ref( ) + 1;


        nucSeqIndex uiRight;
        if( rLast.bOnForwStrand )
            uiRight = rNext.start_ref( );
        else
            uiRight = rLast.start_ref( ) - rLast.size( ) + 1;

        // check if rectangle exists on reference
        if( uiRight >= uiLeft )
            // check if rectangle is small enough
            if( (int64_t)std::max( uiTop > uiBottom ? uiTop - uiBottom : 0, uiRight - uiLeft ) <= iMaxSizeReseed / 2 )
            {
                // check if rectangle is not bridging
                auto uiIDLeft = pRefSeq->uiSequenceIdForPosition( uiLeft );
                auto uiIDRight = pRefSeq->uiSequenceIdForPosition( uiRight );
                if( uiIDLeft == uiIDRight )
                {
                    // check if rectangle is not proportional (too wide)
                    if( uiRight - uiLeft > 2 * ( uiTop > uiBottom ? uiTop - uiBottom : 0 ) && uiRight - uiLeft > 20 )
                    {
                        // return two squares instead
                        auto uiCenter = ( uiRight + uiLeft ) / 2;
                        auto uiSize = uiCenter - uiLeft;
                        return std::make_pair(
                            geom::Rectangle<nucSeqIndex>(
                                uiLeft, rLast.bOnForwStrand ? uiBottom : ( uiTop > uiSize ? uiTop - uiSize : 0 ),
                                uiSize, uiSize ),
                            geom::Rectangle<nucSeqIndex>(
                                uiCenter, rLast.bOnForwStrand ? ( uiTop > uiSize ? uiTop - uiSize : 0 ) : uiBottom,
                                uiSize, uiSize ) );
                    } // if
                    else if( uiTop <= uiBottom )
                    {
                        return std::make_pair( geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ),
                                               geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );
                    } // else if
                    else
                        // can return as a single rectangle
                        return std::make_pair(
                            geom::Rectangle<nucSeqIndex>( uiLeft, uiBottom, uiRight - uiLeft, uiTop - uiBottom ),
                            geom::Rectangle<nucSeqIndex>( 0, 0, 0, 0 ) );
                } // if
                // else triggers return at end of function
            } // if
        // else triggers return at end of function
    } // else if
    // else triggers return at end of function


    // if we reach this point, the rectangle between the seeds is too large or spans over tow contigs
    // so we create two separate rectangles
    return std::make_pair( getPositionsForSeeds( rLast, xDummySeed, rLast.end( ), rNext.start( ), pRefSeq ).first,
                           getPositionsForSeeds( xDummySeed, rNext, rLast.end( ), rNext.start( ), pRefSeq ).first );
} // method


void SvJumpsFromSeeds::computeSeeds( geom::Rectangle<nucSeqIndex>& xArea, std::shared_ptr<NucSeq> pQuery,
                                     std::shared_ptr<Pack> pRefSeq, std::shared_ptr<Seeds> rvRet,
                                     HelperRetVal* pOutExtra )
{
    if( xArea.xXAxis.size( ) == 0 || xArea.xYAxis.size( ) == 0 )
    {
        if( pOutExtra != nullptr )
        {
            pOutExtra->vRectangleReferenceAmbiguity.push_back( 0 );
            pOutExtra->vRectangleUsedDp.push_back( false );
            pOutExtra->vRectangleKMerSize.push_back( 0 );
        } // if
        return;
    } // if
    auto pRef = pRefSeq->vExtract( xArea.xXAxis.start( ), xArea.xXAxis.end( ) );
    auto pRefRevComp = pRefSeq->vExtract( xArea.xXAxis.start( ), xArea.xXAxis.end( ) );
    pRefRevComp->vReverseAll( );
    pRefRevComp->vSwitchAllBasePairsToComplement( );
    // @todo this is inefficient:
    auto pQuerySegment = std::make_shared<NucSeq>( pQuery->fromTo( xArea.xYAxis.start( ), xArea.xYAxis.end( ) ) );
    auto pRevCompQuerySegment = std::make_shared<NucSeq>( *pQuerySegment );
    pRevCompQuerySegment->vReverseAll( );
    pRevCompQuerySegment->vSwitchAllBasePairsToComplement( );
    /*
     * In order to kill random seeds:
     * sample the k-mer size of forward and reverse strand together
     * do query as well as reference
     */
    auto uiSampledAmbiguity =
        std::max( std::max( sampleSequenceAmbiguity( *pRef, *pRef, this->dProbabilityForRandomMatch ),
                            sampleSequenceAmbiguity( *pRef, *pRefRevComp, this->dProbabilityForRandomMatch ) ),
                  std::max( sampleSequenceAmbiguity( *pQuerySegment, *pQuerySegment, this->dProbabilityForRandomMatch ),
                            sampleSequenceAmbiguity( *pQuerySegment, *pRevCompQuerySegment,
                                                     this->dProbabilityForRandomMatch ) ) );
    if( pOutExtra != nullptr )
        pOutExtra->vRectangleReferenceAmbiguity.push_back( uiSampledAmbiguity );
    // ignore reference ambiguity here
    if( true || uiSampledAmbiguity <= xArea.xXAxis.size( ) * ( 1 + dMaxSequenceSimilarity ) )
    {
        HashMapSeeding xHashMapSeeder;
        // purely statistics; i.e. does not look at the sequence at all here
        xHashMapSeeder.uiSeedSize = getKMerSizeForRectangle( xArea, this->dProbabilityForRandomMatch );
        if( pOutExtra != nullptr )
        {
            pOutExtra->vRectangleUsedDp.push_back( false );
            pOutExtra->vRectangleKMerSize.push_back( xHashMapSeeder.uiSeedSize );
        } // if
        if( xHashMapSeeder.uiSeedSize > xArea.xXAxis.size( ) || xHashMapSeeder.uiSeedSize > xArea.xYAxis.size( ) )
            return;

        auto pSeeds = xHashMapSeeder.execute( pQuerySegment, pRef );

        auto pSeedsRev = xHashMapSeeder.execute( pQuerySegment, pRefRevComp );
        // fix seed positions on forward strand
        for( Seed& rSeed : *pSeeds )
        {
            rSeed.uiPosOnReference += xArea.xXAxis.start( );
            rSeed.iStart += xArea.xYAxis.start( );
            assert( rSeed.end( ) <= pQuery->length( ) );
            libMA::ExtractSeeds::setDeltaOfSeed( rSeed, pQuery->length( ), *pRefSeq, true );
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
            libMA::ExtractSeeds::setDeltaOfSeed( rSeed, pQuery->length( ), *pRefSeq, true );
        } // for

        pSeeds->confirmSeedPositions( pQuery, pRefSeq, false );
        pSeedsRev->confirmSeedPositions( pQuery, pRefSeq, false );

        rvRet->append( pSeeds );
        rvRet->append( pSeedsRev );
    }
    else
    {
        if( pOutExtra != nullptr )
        {
            pOutExtra->vRectangleUsedDp.push_back( true );
            pOutExtra->vRectangleKMerSize.push_back( 0 );
        } // if
        return;
        auto pFAlignment = std::make_shared<Alignment>( xArea.xXAxis.start( ), xArea.xYAxis.start( ) );
        AlignedMemoryManager xMemoryManager; // @todo this causes frequent (de-)&allocations; move this outwards
        xNW.ksw( pQuery, pRef, xArea.xYAxis.start( ), xArea.xYAxis.end( ) - 1, 0, pRef->length( ) - 1, pFAlignment,
                 xMemoryManager );
        auto pForwSeeds = pFAlignment->toSeeds( pRefSeq );
        for( Seed& rSeed : *pForwSeeds )
            libMA::ExtractSeeds::setDeltaOfSeed( rSeed, pQuery->length( ), *pRefSeq, true );

        // and now the reverse strand seeds
        auto pRAlignment = std::make_shared<Alignment>( );
        xNW.ksw( pQuery, pRefRevComp, xArea.xYAxis.start( ), xArea.xYAxis.end( ) - 1, 0, pRefRevComp->length( ) - 1,
                 pRAlignment, xMemoryManager );
        auto pRevSeeds = pRAlignment->toSeeds( pRefSeq );
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
            libMA::ExtractSeeds::setDeltaOfSeed( rSeed, pQuery->length( ), *pRefSeq, true );
        } // for
        if( pFAlignment->score( ) >= pRAlignment->score( ) )
            rvRet->append( pForwSeeds );
        else
            rvRet->append( pRevSeeds );
    } // else
} // method

std::shared_ptr<Seeds> SvJumpsFromSeeds::computeSeeds( geom::Rectangle<nucSeqIndex> xArea,
                                                       std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRefSeq,
                                                       HelperRetVal* pOutExtra )
{
    auto pSeeds = std::make_shared<Seeds>( );
    computeSeeds( xArea, pQuery, pRefSeq, pSeeds, pOutExtra );

    if( pSeeds->size( ) == 0 )
        return pSeeds;

    // turn k-mers into maximally extended seeds
    return xSeedLumper.execute( pSeeds, pQuery, pRefSeq );
} // method

std::shared_ptr<Seeds>
SvJumpsFromSeeds::computeSeeds( std::pair<geom::Rectangle<nucSeqIndex>, geom::Rectangle<nucSeqIndex>>& xAreas,
                                std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRefSeq, HelperRetVal* pOutExtra )
{
    auto pSeeds = std::make_shared<Seeds>( );
    computeSeeds( xAreas.first, pQuery, pRefSeq, pSeeds, pOutExtra );
    computeSeeds( xAreas.second, pQuery, pRefSeq, pSeeds, pOutExtra );

    if( pSeeds->size( ) == 0 )
        return pSeeds;

    // turn k-mers into maximally extended seeds
    return xSeedLumper.execute( pSeeds, pQuery, pRefSeq );
} // method


#ifdef WITH_PYTHON
std::shared_ptr<Seeds> computeSeeds( const ParameterSetManager& rParameters, uint32_t uiX, uint32_t uiY, uint32_t uiW,
                                     uint32_t uiH, std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRefSeq )
{
    return SvJumpsFromSeeds( rParameters, pRefSeq )
        .computeSeeds( geom::Rectangle<nucSeqIndex>( uiX, uiY, uiW, uiH ), pQuery, pRefSeq, nullptr );
} // method

void exportSvJumpsFromSeeds( libMS::SubmoduleOrganizer& xOrganizer )
{
    py::class_<geom::Interval<nucSeqIndex>>( xOrganizer.util( ), "nucSeqInterval" )
        .def_readwrite( "start", &geom::Interval<nucSeqIndex>::iStart )
        .def_readwrite( "size", &geom::Interval<nucSeqIndex>::iSize );
    py::class_<geom::Rectangle<nucSeqIndex>>( xOrganizer.util( ), "nucSeqRectangle" )
        .def_readwrite( "x_axis", &geom::Rectangle<nucSeqIndex>::xXAxis )
        .def_readwrite( "y_axis", &geom::Rectangle<nucSeqIndex>::xYAxis );
    py::class_<libMA::HelperRetVal>( xOrganizer._util( ), "SvJumpsFromSeedsHelperRetVal" )
        .def_readwrite( "layer_of_seeds", &libMA::HelperRetVal::vLayerOfSeeds )
        .def_readwrite( "seeds", &libMA::HelperRetVal::pSeeds )
        .def_readwrite( "seed_removed", &libMA::HelperRetVal::pRemovedSeeds )
        .def_readwrite( "rectangles", &libMA::HelperRetVal::vRectangles )
        .def_readwrite( "rectangle_layers", &libMA::HelperRetVal::vRectangleLayers )
        .def_readwrite( "palindrome", &libMA::HelperRetVal::vParlindromeSeed )
        .def_readwrite( "overlapping", &libMA::HelperRetVal::vOverlappingSeed )
        .def_readwrite( "rectangles_fill", &libMA::HelperRetVal::vRectangleFillPercentage )
        .def_readwrite( "rectangle_ambiguity", &libMA::HelperRetVal::vRectangleReferenceAmbiguity )
        .def_readwrite( "rectangle_k_mer_size", &libMA::HelperRetVal::vRectangleKMerSize )
        .def_readwrite( "rectangle_used_dp", &libMA::HelperRetVal::vRectangleUsedDp )
        .def_readwrite( "jump_seeds", &libMA::HelperRetVal::vJumpSeeds );
    py::bind_vector<std::vector<geom::Rectangle<nucSeqIndex>>>( xOrganizer.util( ), "RectangleVector", "" );
    py::bind_vector<std::vector<double>>( xOrganizer._util( ), "DoubleVector", "" );
    py::bind_vector<std::vector<bool>>( xOrganizer._util( ), "BoolVector", "" );

    py::bind_vector_ext<ContainerVector<SvJump>, libMS::Container, std::shared_ptr<ContainerVector<SvJump>>>(
        xOrganizer._util( ), "JumpVector", "" );

    exportModule<SvJumpsFromSeeds, std::shared_ptr<Pack>>( xOrganizer, "SvJumpsFromSeeds", []( auto&& x ) {
        x.def( "execute_helper", &SvJumpsFromSeeds::execute_helper_py );
        x.def( "execute_helper", &SvJumpsFromSeeds::execute_helper_py2 );
        x.def( "execute_helper_no_reseed", &SvJumpsFromSeeds::execute_helper_py3 );
        x.def( "compute_jumps", &SvJumpsFromSeeds::computeJumpsPy );
    } );
    exportModule<RecursiveReseeding, std::shared_ptr<Pack>>( xOrganizer, "RecursiveReseeding" );
    exportModule<RecursiveReseedingSoCs, std::shared_ptr<Pack>>( xOrganizer, "RecursiveReseedingSoCs", []( auto&& x ) {
        x.def( "execute_helper", &RecursiveReseedingSoCs::execute_helper_py );
    } );
    exportModule<JumpsFilterContigBorder>( xOrganizer, "JumpsFilterContigBorder", []( auto&& x ) {
        x.def( "by_contig_border", &JumpsFilterContigBorder::jumpByContigBorder );
    } );
    exportModule<SvJumpsFromExtractedSeeds, std::shared_ptr<Pack>>( xOrganizer, "SvJumpsFromExtractedSeeds" );
    exportModule<ExtractSeedsFilter, std::shared_ptr<Pack>, nucSeqIndex, nucSeqIndex>( xOrganizer,
                                                                                       "ExtractSeedsFilter" );
    exportModule<RecursiveReseedingSegments, std::shared_ptr<Pack>>( xOrganizer, "RecursiveReseedingSegments" );
    exportModule<FilterJumpsByRegion, int64_t, int64_t>( xOrganizer, "FilterJumpsByRegion" );
    exportModule<FilterJumpsByRegionSquare, int64_t, int64_t>( xOrganizer, "FilterJumpsByRegionSquare" );
    exportModule<FilterJumpsByRefAmbiguity>( xOrganizer, "FilterJumpsByRefAmbiguity" );
    exportModule<CollectSeedCoverage>( xOrganizer, "CollectSeedCoverage" );

    xOrganizer.util( ).def( "compute_seeds_area", &computeSeeds );
} // function
#endif