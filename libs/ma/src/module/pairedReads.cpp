/**
 * @file pairedReads.cpp
 * @author Markus Schmidt
 */
#ifdef _MSC_VER
#define NOMINMAX
#endif

#include "module/pairedReads.h"

using namespace libMA;

using namespace libMA::defaults;
extern int libMA::defaults::iMatch;

#define SQ_12 std::sqrt( 12 )

double PairedReads::p( nucSeqIndex d ) const // TO DO inline me
{
    if( bNormalDist )
        return 1 - 0.499 * ( 1 + std::erf( ( d - mean ) / ( std * std::sqrt( 2 ) ) ) );
    if( bUniformDist )
        /*
         * uniform distribution has following CDF:
         * (d - a) / (b - a)
         * where a and b are min,max.
         * the mean is thus: (b+a)/2
         * the std is: (b-a)*sqrt(12)
         */
        return 1 - std::max( 0.0, std::min( 1.0, ( d - mean - std / ( SQ_12 * 2 ) ) / ( std / SQ_12 ) ) );

    return 0.0;
} // function

std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
PairedReads::execute( std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> pAlignments1,
                      std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
                          pAlignments2,
                      std::shared_ptr<Pack>
                          pPack )
{

    // assert that we have at least one alignment for each mate.
    if( pAlignments1->size( ) == 0 )
        return pAlignments2;
    if( pAlignments2->size( ) == 0 )
        return pAlignments1;

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

            // get the correct start position on the forward strand

            nucSeqIndex uiP1 = pAlignment1->beginOnRef( );
            if( pPack->bPositionIsOnReversStrand( pAlignment1->beginOnRef( ) ) )
                uiP1 = pPack->iAbsolutePosition( pAlignment1->endOnRef( ) );

            // get the correct start position on the forward strand
            nucSeqIndex uiP2 = pAlignment2->beginOnRef( );
            if( pPack->bPositionIsOnReversStrand( pAlignment2->beginOnRef( ) ) )
                uiP2 = pPack->iAbsolutePosition( pAlignment2->endOnRef( ) );


            // get the distance of the alignments on the reference
            nucSeqIndex d = uiP1 < uiP2 ? uiP2 - uiP1 : uiP1 - uiP2;
            // std::cerr << "A1: " << pAlignment1->xStats.sName
            //           << " A2: " << pAlignment2->xStats.sName << std::endl;
            // if( d < 1000 )
            //     std::cout << "Candidates found" << std::endl;

            // compute the score by the formula given in:
            // "Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM"
            // q = (int)((v.a[i].y>>32) + (v.a[k].y>>32) + .721 * log(2. * erfc(fabs(ns) *
            // M_SQRT1_2)) * opt->iMatchScore + .499);

            // (Heng Li)
            double score = pAlignment1->score( ) + pAlignment2->score( ) -
                           std::min( -iMatch * std::log( p( d ) ) / std::log( 4 ), (double)u );


            // std::cerr << "score: " << score << std::endl;

            // check if we have a new best pair
            if( score > maxScore )
            {
                // if so update the score and indices
                maxScore = score;
                uiI1 = i;
                uiI2 = j;
                // it's possible that the best pair is two individual alignments

                // std::cerr << "s2: " << -iMatch * std::log(p(d))/std::log(4) << std::endl;
                // if (-iMatch * std::log(p(d))/std::log(4) > 0)
                // {
                //     std::cerr << " > 0" << std::endl;
                //     exit(0);
                // }

                bPaired = -iMatch * std::log( p( d ) ) / std::log( 4 ) < u;
                // if( !bPaired)
                //     std::cerr << "Distance: " << d << "Paired: " << bPaired << std::endl;
            } // if
        } // for
    } // for


    // set the paired property in the respective alignment stats
    if( bPaired )
    {
        // set which read was first..

        ( *pAlignments1 )[ uiI1 ]->xStats.bFirst = true; // dc
        ( *pAlignments2 )[ uiI2 ]->xStats.bFirst = false; // dc

        // { // temp code
        //     std::shared_ptr<Alignment> pAlignment1 = std::dynamic_pointer_cast<Alignment>(//
        //     (*pAlignments1)[uiI1]); // dc std::shared_ptr<Alignment> pAlignment2 =
        //     std::dynamic_pointer_cast<Alignment>(// (*pAlignments2)[uiI2]); // dc nucSeqIndex
        //     uiP1 = pAlignment1->beginOnRef(); if
        //     (pPack->bPositionIsOnReversStrand(pAlignment1->beginOnRef()))
        //         uiP1 = pPack->iAbsolutePosition(pAlignment1->endOnRef());
        //
        //     // get the correct start position on the forward strand
        //     nucSeqIndex uiP2 = pAlignment2->beginOnRef();
        //     if (pPack->bPositionIsOnReversStrand(pAlignment2->beginOnRef()))
        //         uiP2 = pPack->iAbsolutePosition(pAlignment2->endOnRef());
        //
        //     // get the distance of the alignments on the reference
        //     nucSeqIndex d = uiP1 < uiP2 ? uiP2 - uiP1 : uiP1 - uiP2;
        //     if (d > 4000)
        //         std::cerr << "WARNING: " << d << std::endl;
        // }

        // std::cout << " uiI1 " <<
        // std::dynamic_pointer_cast<Alignment>((*pAlignments1)[uiI1])->xStats.bFirst << "\n"
        //          << " uiI2 " <<
        //          std::dynamic_pointer_cast<Alignment>((*pAlignments2)[uiI2])->xStats.bFirst <<
        //          std::endl;
        for( auto pX : *pAlignments1 )
            pX->bSecondary = true;
        ( *pAlignments1 )[ uiI1 ]->bSecondary = false;
        ( *pAlignments1 )[ uiI1 ]->xStats.pOther = std::weak_ptr<Alignment>( ( *pAlignments2 )[ uiI2 ] );

        for( auto pX : *pAlignments2 )
            pX->bSecondary = true;
        ( *pAlignments2 )[ uiI2 ]->xStats.pOther = std::weak_ptr<Alignment>( ( *pAlignments1 )[ uiI1 ] );
        ( *pAlignments2 )[ uiI2 ]->bSecondary = false;
        
    } // if

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