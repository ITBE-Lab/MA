/**
 * @file stripOfConsideration.cpp
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */
#include "ma/module/stripOfConsideration.h"
#include "util/system.h"
using namespace libMA;
using namespace libMS;


std::shared_ptr<SoCPriorityQueue> StripOfConsiderationSeeds::execute(
    std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuerySeq, std::shared_ptr<Pack> pRefSeq )
{
    // make sure that we return at least an SoC set
    if( pSeeds->empty( ) )
        return std::make_shared<SoCPriorityQueue>( pSeeds );

    const nucSeqIndex uiQLen = pQuerySeq->length( );

    double fMinLen = std::max( (double)fGiveUp * uiQLen, (double)this->uiCurrHarmScoreMin );
    if( uiMinGenomeSize >= pRefSeq->uiUnpackedSizeForwardPlusReverse( ) )
        fMinLen = 0;

    /*
     * This is the formula from the paper
     * computes the size required for the strip so that we collect all relevent seeds.
     */
    const nucSeqIndex uiStripSize = this->getStripSize( uiQLen, pGlobalParams->iMatch->get( ),
                                                        pGlobalParams->iExtend->get( ), pGlobalParams->iGap->get( ) );

    // sort the seeds according to their initial positions
    std::sort( pSeeds->begin( ), pSeeds->end( ),
               [ & ]( const Seed& a, const Seed& b ) { return a.uiDelta < b.uiDelta; } // lambda
    ); // sort function call

#if DEBUG_LEVEL > 0
    // check order of seeds to make sure optimization is valid
    if( !bRectangular )
    {
        bool bNotSeedRevStrandYet = true;
        for( size_t uiI = 0; uiI < pSeeds->size( ); uiI++ )
        {
            if( !bNotSeedRevStrandYet && ( *pSeeds )[ uiI ].bOnForwStrand )
            {
                std::cerr << "seeds from forward and reverse strand are intermangled!" << std::endl;
                assert( false );
            } // if
            if( bNotSeedRevStrandYet && !( *pSeeds )[ uiI ].bOnForwStrand )
                bNotSeedRevStrandYet = false;
        } // for
    } // if
#endif

    // positions to remember the maxima
    auto pSoCs = std::make_shared<SoCPriorityQueue>( pSeeds );
    DEBUG( pSoCs->uiQLen = uiQLen; ) // DEBUG

    // find the SOC maxima
    SoCOrder xCurrScore;
    std::vector<Seed>::iterator xStripStart = pSeeds->begin( );
    std::vector<Seed>::iterator xStripEnd = pSeeds->begin( );
    size_t uiContigIdStart = pRefSeq->uiSequenceIdForPosition( xStripStart->start_ref( ) );
    size_t uiContigIdEnd = pRefSeq->uiSequenceIdForPosition( xStripEnd->start_ref( ) );
    assert( pRefSeq->isForwPositionInSequenceWithId( uiContigIdStart, xStripStart->start_ref( ) ) );
    assert( pRefSeq->isForwPositionInSequenceWithId( uiContigIdEnd, xStripEnd->start_ref( ) ) );
    while( xStripEnd != pSeeds->end( ) && xStripStart != pSeeds->end( ) )
    {
        // xStripStart might have moved forward so adjust uiContigIdStart.
        if( bRectangular || xStripStart->bOnForwStrand )
            assert( (int64_t)uiContigIdStart <= pRefSeq->uiSequenceIdForPosition( xStripStart->start_ref( ) ) );
        else
            assert( (int64_t)uiContigIdStart >= pRefSeq->uiSequenceIdForPosition( xStripStart->start_ref( ) ) );
        // adjust the contig id
        // this while loop is a runtime optimization as it avoids the need of a binary search for the contig id
        // (since we know that the contig id between two seeds in a soc does not change often)
        while( !pRefSeq->isForwPositionInSequenceWithId( uiContigIdStart, xStripStart->start_ref( ) ) )
            uiContigIdStart += ( bRectangular || xStripStart->bOnForwStrand ) ? 1 : -1;

        // move xStripEnd forwards while it is closer to xStripStart than uiStripSize
        while( xStripEnd != pSeeds->end( ) && xStripStart->uiDelta + uiStripSize >= xStripEnd->uiDelta &&
               uiContigIdStart == uiContigIdEnd &&
               ( bRectangular || ( xStripStart->bOnForwStrand == xStripEnd->bOnForwStrand ) ) )
        {
            // remember the additional score
            xCurrScore += *xStripEnd;
            // move the iterator forward
            xStripEnd++;
            if( xStripEnd != pSeeds->end( ) )
            {
                /* if we have a non rectangular SoC where we just moved xStripEnd from the forward to the reverse strand
                 * then we have to set uiContigIdEnd to the largest contig id.
                 * Otherwise we have a bug with counting down contig ids if:
                 * There is no seed for the last contig of the forward strand but one seed for the last contig of the
                 * reverse strand. In this case the uiContigIdEnd would have to be counted upwards...
                 * We avoid this by forcing uiContigIdEnd to be the highest possible number.
                 */
                if( !bRectangular && !xStripEnd->bOnForwStrand && xStripStart->bOnForwStrand )
                    uiContigIdEnd = pRefSeq->uiNumContigs( ) - 1; // force to last contig so that we can count backwards
#if DEBUG_LEVEL > 0
                // code in here makes sure that our optimization below does not mess up
                int64_t uiExpectedContigIdEnd = pRefSeq->uiSequenceIdForPosition( xStripEnd->start_ref( ) );
                if( ( bRectangular || xStripEnd->bOnForwStrand ) && (int64_t)uiContigIdEnd > uiExpectedContigIdEnd )
                {
                    std::cerr << "got weird contig id (1): " << uiContigIdEnd << " " << uiExpectedContigIdEnd
                              << std::endl;
                    assert( false );
                } // if
                else if( !( bRectangular || xStripEnd->bOnForwStrand ) &&
                         (int64_t)uiContigIdEnd < uiExpectedContigIdEnd )
                {
                    std::cerr << "got weird contig id (2): " << uiContigIdEnd << " " << uiExpectedContigIdEnd
                              << std::endl;
                    assert( false );
                } // else if
#endif
                // adjust the contig id
                // this while loop is a runtime optimization as it avoids the need of a binary search for the contig id
                // (since we know that the contig id between two seeds in a soc does not change often)
                while( !pRefSeq->isForwPositionInSequenceWithId( uiContigIdEnd, xStripEnd->start_ref( ) ) )
                    uiContigIdEnd += ( bRectangular || xStripEnd->bOnForwStrand ) ? 1 : -1;
            } // if
        } // while


        DEBUG( pSoCs->vScores.push_back( std::make_pair( xCurrScore.uiAccumulativeLength,
                                                         xStripStart->start_ref( ) ) ); ) // DEBUG

        // FILTER
        /*
         * if the SoC quality is lower than fGiveUp * uiQLen we do not consider this SoC at
         * all fGiveUp = 0 disables this.
         */
        if( fGiveUp == 0 || xCurrScore.uiAccumulativeLength >= fMinLen )
            pSoCs->push_back_no_overlap( xCurrScore, xStripStart, xStripEnd, xStripStart->uiDelta,
                                         ( xStripEnd - 1 )->uiDelta );
        // move xStripStart one to the right (this will cause xStripEnd to be adjusted)
        bool bLastOnForw = xStripStart->bOnForwStrand;
        xCurrScore -= *( xStripStart++ );
        /* if we have a non rectangular SoC where we just moved xStripStart from the forward to the reverse strand
         * then we have to set uiContigIdStart to the largest contig id.
         * Otherwise we have a bug with counting down contig ids if:
         * There is no seed for the last contig of the forward strand but one seed for the last contig of the
         * reverse strand. In this case the uiContigIdStart would have to be counted upwards...
         * We avoid this by forcing uiContigIdStart to be the highest possible number.
         */
        if( !bRectangular && xStripStart != pSeeds->end( ) && xStripEnd != pSeeds->end( ) &&
            !xStripStart->bOnForwStrand && bLastOnForw )
            uiContigIdStart = pRefSeq->uiNumContigs( ) - 1; // force to last contig so that we can count backwards
    } // while

    // make a max heap from the SOC starting points according to the scores,
    // so that we can extract the best SOC first
    pSoCs->make_heap( );


    if( bRectangular )
        pSoCs->rectangularSoC( );

    // return the strip collection
    return pSoCs;

} // function
std::shared_ptr<SoCPriorityQueue> StripOfConsideration::execute( std::shared_ptr<SegmentVector> pSegments,
                                                                 std::shared_ptr<NucSeq>
                                                                     pQuerySeq,
                                                                 std::shared_ptr<Pack>
                                                                     pRefSeq,
                                                                 std::shared_ptr<FMIndex>
                                                                     pFM_index )
{
    auto pSeeds = xExtractHelper.execute( pSegments, pFM_index, pQuerySeq, pRefSeq );

    return xHelper.execute( pSeeds, pQuerySeq, pRefSeq );
} // function

#ifdef WITH_PYTHON

void exportStripOfConsideration( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the StripOfConsideration class
    exportModule<StripOfConsiderationSeeds>( xOrganizer, "StripOfConsiderationSeeds" );
    exportModule<StripOfConsideration>( xOrganizer, "StripOfConsideration" );
    exportModule<ExtractSeeds>( xOrganizer, "ExtractSeeds" );
    exportModule<GetAllFeasibleSoCs>( xOrganizer, "GetAllFeasibleSoCs" );
    exportModule<GetAllFeasibleSoCsAsSet>( xOrganizer, "GetAllFeasibleSoCsAsSet" );
} // function
#endif