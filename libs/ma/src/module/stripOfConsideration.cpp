/**
 * @file stripOfConsideration.cpp
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */
#include "ma/module/stripOfConsideration.h"
#include "util/system.h"
using namespace libMA;
using namespace libMS;


std::shared_ptr<SoCPriorityQueue> StripOfConsiderationSeeds::execute( std::shared_ptr<Seeds> pSeeds,
                                                                      std::shared_ptr<NucSeq>
                                                                          pQuerySeq,
                                                                      std::shared_ptr<Pack>
                                                                          pRefSeq,
                                                                      std::shared_ptr<FMIndex>
                                                                          pFM_index )
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
               [&]( const Seed& a, const Seed& b ) { return a.uiDelta < b.uiDelta; } // lambda
    ); // sort function call


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
        assert( (int64_t)uiContigIdStart <= pRefSeq->uiSequenceIdForPosition( xStripStart->start_ref( ) ) );
        // adjust the contig id
        while( !pRefSeq->isForwPositionInSequenceWithId( uiContigIdStart, xStripStart->start_ref( ) ) )
            uiContigIdStart++;

        // move xStripEnd forwards while it is closer to xStripStart than uiStripSize
        while( xStripEnd != pSeeds->end( ) &&
               getPositionForBucketing( uiQLen, *xStripStart ) + uiStripSize >=
                   getPositionForBucketing( uiQLen, *xStripEnd ) &&
               uiContigIdStart == uiContigIdEnd )
        {
            // remember the additional score
            xCurrScore += *xStripEnd;
            // move the iterator forward
            xStripEnd++;
            if( xStripEnd != pSeeds->end( ) )
            {
#if DEBUG_LEVEL > 0
                int64_t uiExpectedContigIdEnd = pRefSeq->uiSequenceIdForPosition( xStripEnd->start_ref( ) );
                if( (int64_t)uiContigIdEnd > uiExpectedContigIdEnd )
                {
                    std::cerr << "got wierd contig id: " << uiContigIdEnd << " " << uiExpectedContigIdEnd << std::endl;
                    assert( false );
                } // if
#endif
                // adjust the contig id
                while( !pRefSeq->isForwPositionInSequenceWithId( uiContigIdEnd, xStripEnd->start_ref( ) ) )
                    uiContigIdEnd++;
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
            pSoCs->push_back_no_overlap( xCurrScore, xStripStart, xStripEnd,
                                         getPositionForBucketing( uiQLen, *xStripStart ),
                                         getPositionForBucketing( uiQLen, *( xStripEnd - 1 ) ) );
        // move xStripStart one to the right (this will cause xStripEnd to be adjusted)
        xCurrScore -= *( xStripStart++ );
    } // while

    // make a max heap from the SOC starting points according to the scores,
    // so that we can extract the best SOC first
    pSoCs->make_heap( );


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

    return xHelper.execute( pSeeds, pQuerySeq, pRefSeq, pFM_index );
} // function

#ifdef WITH_PYTHON

void exportStripOfConsideration( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the StripOfConsideration class
    exportModule<StripOfConsiderationSeeds>( xOrganizer, "StripOfConsiderationSeeds" );
    exportModule<StripOfConsideration>( xOrganizer, "StripOfConsideration" );
    exportModule<ExtractSeeds>( xOrganizer, "ExtractSeeds" );
    exportModule<GetAllFeasibleSoCs>( xOrganizer, "GetAllFeasibleSoCs" );
} // function
#endif