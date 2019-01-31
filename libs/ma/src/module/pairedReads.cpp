/**
 * @file pairedReads.cpp
 * @author Markus Schmidt
 */
#ifdef _MSC_VER
#define NOMINMAX
#endif

#include <limits>
#include "module/pairedReads.h"

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
    //std::cerr << "new read" << std::endl;

    // the indices of the best alignment pair
    unsigned int uiI1 = 0, uiI2 = 0;
    double maxScore = 0;
    bool bPaired = false;
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

            // illumina reads must be on opposite strands
            if( pPack->bPositionIsOnReversStrand( pAlignment1->beginOnRef( ) ) ==
                pPack->bPositionIsOnReversStrand( pAlignment2->beginOnRef( ) ) )
                continue;

            nucSeqIndex uiP1 = pAlignment1->beginOnRef( );
            // for illumina the reads are always on opposite strands
            nucSeqIndex uiP2 = pPack->uiPositionToReverseStrand( pAlignment2->beginOnRef( ) );

            // get the distance of the alignments on the reference
            nucSeqIndex d = uiP1 < uiP2 ? uiP2 - uiP1 : uiP1 - uiP2;

            //      r->low  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
            //      if (r->low < 1) r->low = 1;
            //      r->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);

            // std::cerr << "A1: " << pAlignment1->xStats.sName
            //           << " A2: " << pAlignment2->xStats.sName << std::endl;
            // if( d < 1000 )
            //     std::cout << "Candidates found" << std::endl;

            // compute the score by the formula given in:
            // "Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM"
            // q = (int)((v.a[i].y>>32) + (v.a[k].y>>32) + .721 * log(2. * erfc(fabs(ns) *
            // M_SQRT1_2)) * opt->iMatchScore + .499);

            // (Heng Li)
            double dPairedScore =
                pAlignment1->score( ) + pAlignment2->score( ) +
                0.721 * std::log( 2. * std::erfc( std::abs( ( d - mean ) / std ) / std::sqrt( 2 ) ) ) * iMatch + .499;
            if( dPairedScore < 0 )
                dPairedScore = 0;


            // std::cerr << "score: " << score << std::endl;

            // check if we have a new best pair
            //std::cerr << dPairedScore << " to " << pAlignment1->score( ) + pAlignment2->score( ) - dUnpaired << " pref " << maxScore << " blub "<< 2. * std::erfc( std::abs( ( d - mean ) / std ) / std::sqrt( 2 ) ) << std::endl;
            if( dPairedScore > maxScore )
            {
                // if so update the score and indices
                maxScore = dPairedScore;
                uiI1 = i;
                uiI2 = j;
                // it's possible that the best pair is two individual alignments
                bPaired = dPairedScore >= pAlignment1->score( ) + pAlignment2->score( ) - dUnpaired;
                bPaired = true;
            } // if
            // if( bPaired)
        } // for
    } // for


    // set the paired property in the respective alignment stats
    if( bPaired )
    {
        // set which read was first..

        ( *pAlignments1 )[ uiI1 ]->xStats.bFirst = true; // dc
        ( *pAlignments2 )[ uiI2 ]->xStats.bFirst = false; // dc
        ( *pAlignments1 )[ uiI1 ]->bSecondary = false;
        ( *pAlignments2 )[ uiI2 ]->bSecondary = false;
        ( *pAlignments1 )[ uiI1 ]->xStats.pOther = std::weak_ptr<Alignment>( ( *pAlignments2 )[ uiI2 ] );
        ( *pAlignments2 )[ uiI2 ]->xStats.pOther = std::weak_ptr<Alignment>( ( *pAlignments1 )[ uiI1 ] );

        auto pRet = std::make_shared<ContainerVector<std::shared_ptr<Alignment>>>( );
        pRet->push_back( ( *pAlignments1 )[ uiI1 ] );
        pRet->push_back( ( *pAlignments2 )[ uiI2 ] );
        return pRet;
    } // if
    // if we do not have paired reads...

    // return the alignments
    auto pRet = std::make_shared<ContainerVector<std::shared_ptr<Alignment>>>( );
    for( auto pX : *pAlignments1 )
        pRet->push_back( pX );
    for( auto pX : *pAlignments2 )
        pRet->push_back( pX );
    return pRet;
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