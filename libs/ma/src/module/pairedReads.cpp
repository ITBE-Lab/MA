/**
 * @file pairedReads.cpp
 * @author Markus Schmidt
 */
#ifdef _MSC_VER
#define NOMINMAX
#endif

#include "module/pairedReads.h"
#include <limits>

using namespace libMA;

using namespace libMA::defaults;
extern int libMA::defaults::iMatch;

std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
PairedReads::execute( std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2,
                      std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> pAlignments1,
                      std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> pAlignments2,
                      std::shared_ptr<Pack> pPack )
{

    // assert that we have at least one alignment for each mate.
    if( pAlignments1->size( ) == 0 )
        return pAlignments2;
    if( pAlignments2->size( ) == 0 )
        return pAlignments1;
        // std::cerr << "new read" << std::endl;
#if 1
    // the indices of the best alignment pair
    std::vector<std::tuple<int64_t, bool, size_t, size_t>> vScores;
    // try out all combinations
    // this assumes that not more than three or four alignments are reported
    // otherwise this here might be a serious bottleneck
    for( unsigned int i = 0; i < pAlignments1->size( ); i++ )
    {
        std::shared_ptr<Alignment> pAlignment1 = ( *pAlignments1 )[ i ]; // dc
        if( pAlignment1->length( ) == 0 )
            continue;
        for( unsigned int j = 0; j < pAlignments2->size( ); j++ )
        {
            std::shared_ptr<Alignment> pAlignment2 = ( *pAlignments2 )[ j ]; // dc
            if( pAlignment2->length( ) == 0 )
                continue;

            int64_t iScore = pAlignment1->score( ) + pAlignment2->score( );
            bool bIsPaired = false;
            // illumina reads must be on opposite strands
            if( pPack->bPositionIsOnReversStrand( pAlignment1->beginOnRef( ) ) !=
                pPack->bPositionIsOnReversStrand( pAlignment2->beginOnRef( ) ) )
            {
                nucSeqIndex uiP1 = pAlignment1->beginOnRef( );
                // for illumina the reads are always on opposite strands
                nucSeqIndex uiP2 = pPack->uiPositionToReverseStrand( pAlignment2->beginOnRef( ) );

                // get the distance of the alignments on the reference
                nucSeqIndex d = uiP1 < uiP2 ? uiP2 - uiP1 : uiP1 - uiP2;

                // bonus score if alignments are paired
                if( ( (double)d ) >= ( (double)mean ) - std * 3 && ( (double)d ) <= ( (double)mean ) + std * 3 )
                {
                    iScore *= dUnpaired;
                    bIsPaired = true;
                } // if
            }
            vScores.emplace_back( iScore, bIsPaired, i, j );
        } // for
    } // for

    std::sort( vScores.begin( ), vScores.end( ),
               []( const std::tuple<int64_t, bool, size_t, size_t>& rtA,
                   const std::tuple<int64_t, bool, size_t, size_t>& rtB ) {
                   if( std::get<0>( rtA ) == std::get<0>( rtB ) )
                       return std::get<1>( rtA ) && !std::get<1>( rtB );
                   return std::get<0>( rtA ) > std::get<0>( rtB );
               } );
    assert( vScores.size( ) <= 1 || std::get<0>( vScores[ 0 ] ) >= std::get<0>( vScores[ 1 ] ) );

    // set the paired property in the respective alignment stats
    // set which read was first..
    size_t uiI1 = std::get<2>( vScores[ 0 ] );
    size_t uiI2 = std::get<3>( vScores[ 0 ] );
    ( *pAlignments1 )[ uiI1 ]->xStats.bFirst = true;
    ( *pAlignments2 )[ uiI2 ]->xStats.bFirst = false;
    ( *pAlignments1 )[ uiI1 ]->bSecondary = false;
    ( *pAlignments2 )[ uiI2 ]->bSecondary = false;
    ( *pAlignments1 )[ uiI1 ]->xStats.pOther = std::weak_ptr<Alignment>( ( *pAlignments2 )[ uiI2 ] );
    ( *pAlignments2 )[ uiI2 ]->xStats.pOther = std::weak_ptr<Alignment>( ( *pAlignments1 )[ uiI1 ] );
    if( std::get<1>( vScores[ 0 ] ) && vScores.size( ) > 1 )
    {
        float fMapQ =
            ( (float)( std::get<0>( vScores[ 0 ] ) - std::get<0>( vScores[ 1 ] ) ) ) / std::get<0>( vScores[ 0 ] );

        if(( *pAlignments1 )[ uiI1 ]->getNumSeeds() <= 1)
            fMapQ /= 2;
        if(( *pAlignments2 )[ uiI2 ]->getNumSeeds() <= 1)
            fMapQ /= 2;

        if(( *pAlignments1 )[ uiI1 ]->score() >= iMatch * pQ1->length( ) * 0.8 && pAlignments1->size() >= 3)
            fMapQ *= 2;
        if(( *pAlignments2 )[ uiI2 ]->score() >= iMatch * pQ2->length( ) * 0.8 && pAlignments2->size() >= 3)
            fMapQ *= 2;
        if(fMapQ > 1)
            fMapQ = 1;
    
        ( *pAlignments1 )[ uiI1 ]->fMappingQuality = fMapQ;
        ( *pAlignments2 )[ uiI2 ]->fMappingQuality = fMapQ;
    } // if

    auto pRet = std::make_shared<ContainerVector<std::shared_ptr<Alignment>>>( );
    pRet->push_back( ( *pAlignments1 )[ uiI1 ] );
    pRet->push_back( ( *pAlignments2 )[ uiI2 ] );

    return pRet;

    // set the paired property in the respective alignment stats
    // if( bPaired )
#else
    {
        // set which read was first..

        ( *pAlignments1 )[ 0 ]->xStats.bFirst = true; // dc
        ( *pAlignments2 )[ 0 ]->xStats.bFirst = false; // dc
        ( *pAlignments1 )[ 0 ]->xStats.pOther = std::weak_ptr<Alignment>( ( *pAlignments2 )[ 0 ] );
        ( *pAlignments2 )[ 0 ]->xStats.pOther = std::weak_ptr<Alignment>( ( *pAlignments1 )[ 0 ] );

        nucSeqIndex uiP1 = ( *pAlignments1 )[ 0 ]->beginOnRef( );
        // for illumina the reads are always on opposite strands
        nucSeqIndex uiP2 = pPack->uiPositionToReverseStrand( ( *pAlignments2 )[ 0 ]->beginOnRef( ) );

        // get the distance of the alignments on the reference
        nucSeqIndex d = uiP1 < uiP2 ? uiP2 - uiP1 : uiP1 - uiP2;
        if( ( (double)d ) >= ( (double)mean ) - std * 3 && ( (double)d ) <= ( (double)mean ) + std * 3 )
        {

            ( *pAlignments1 )[ 0 ]->fMappingQuality *= 2;
            if( ( *pAlignments1 )[ 0 ]->fMappingQuality < 0.1 )
                ( *pAlignments1 )[ 0 ]->fMappingQuality = 0.1;
            if( ( *pAlignments1 )[ 0 ]->fMappingQuality > 1 )
                ( *pAlignments1 )[ 0 ]->fMappingQuality = 1;

            ( *pAlignments2 )[ 0 ]->fMappingQuality *= 2;
            if( ( *pAlignments2 )[ 0 ]->fMappingQuality < 0.1 ) // @todo this can push us over!!!
                ( *pAlignments2 )[ 0 ]->fMappingQuality = 0.1;
            if( ( *pAlignments2 )[ 0 ]->fMappingQuality > 1 )
                ( *pAlignments2 )[ 0 ]->fMappingQuality = 1;
        } // if

        auto pRet = std::make_shared<ContainerVector<std::shared_ptr<Alignment>>>( );
        pRet->push_back( ( *pAlignments1 )[ 0 ] );
        pRet->push_back( ( *pAlignments2 )[ 0 ] );
        return pRet;
    } // if
#endif
#if 0
    // return the alignments
    auto pRet = std::make_shared<ContainerVector<std::shared_ptr<Alignment>>>( );
    for( auto pX : *pAlignments1 )
        pRet->push_back( pX );
    for( auto pX : *pAlignments2 )
        pRet->push_back( pX );
    return pRet;
#endif
} // function

#ifdef WITH_PYTHON

#ifdef BOOST_PYTHON
void exportPairedReads( )
{
    // export the PairedReads class
    exportModule<PairedReads>( "PairedReads" );
} // function
#else
void exportPairedReads( py::module& rxPyModuleId )
{
    // export the PairedReads class
    exportModule<PairedReads>( rxPyModuleId, "PairedReads" );
} // function
#endif

#endif