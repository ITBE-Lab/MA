/**
 * @file stripOfConsideration.cpp
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */

#include "module/stripOfConsideration.h"
#include "util/radix_sort.h"
#include "util/system.h"
using namespace libMA;

using namespace libMA::defaults;
extern int libMA::defaults::iGap;
extern int libMA::defaults::iExtend;
extern int libMA::defaults::iMatch;
extern int libMA::defaults::iMissMatch;


std::shared_ptr<SoCPriorityQueue> StripOfConsideration::execute( std::shared_ptr<SegmentVector> pSegments,
                                                                 std::shared_ptr<NucSeq>
                                                                     pQuerySeq,
                                                                 std::shared_ptr<Pack>
                                                                     pRefSeq,
                                                                 std::shared_ptr<FMIndex>
                                                                     pFM_index )
{
    const nucSeqIndex uiQLen = pQuerySeq->length( );

    double fMinLen = std::max( (double)fGiveUp * uiQLen, (double)this->uiCurrHarmScoreMin );
    if( uiMinGenomeSize >= pRefSeq->uiUnpackedSizeForwardPlusReverse( ) )
        fMinLen = 0;

    /*
     * This is the formula from the paper
     * computes the size required for the strip so that we collect all relevent seeds.
     */
    const nucSeqIndex uiStripSize = this->getStripSize( uiQLen, iMatch, iExtend, iGap );

    // extract the seeds
    auto pSeeds = std::make_shared<Seeds>( );
    pSeeds->xStats.sName = pQuerySeq->sName;
    pSeeds->xStats.bSetMappingQualityToZero = pSegments->bSetMappingQualityToZero;
    // rough estimate of how many seeds we will have
    // (trying to avoid multiple allocations)
    pSeeds->reserve( pSegments->size( ) * 3 );

#if MEASURE_DURATIONS == ( 1 )
    pSegments->fExtraction += metaMeasureDuration( [&]( ) {
#endif
                                  emplaceAllNonBridgingSeed( *pSegments, // Segment vector (outcome of seeding)
                                                             *pFM_index, *pRefSeq, *pSeeds,
                                                             uiQLen ); // emplace all function call
#if MEASURE_DURATIONS == ( 1 )
                              } // lambda
                                                   )
                                  .count( ) *
                              1000; // metaMeasureDuration function call
#endif


    // make sure that we return at least an SoC set
    if( pSeeds->empty( ) )
        return std::make_shared<SoCPriorityQueue>( pSeeds );

#if MEASURE_DURATIONS == ( 1 )
    pSegments->fSorting += metaMeasureDuration( [&]( ) {
#endif

#if DELTA_CACHE == ( 1 )
    // delta values already set by emplaceAllNonBridgingSeed
    radix_sort<8>( &pSeeds->front( ), &pSeeds->front( ) + pSeeds->size( ),
                    []( const Seed& rSeed ) { return rSeed.uiDelta; } );
#else
    // sort the seeds according to their initial positions
    std::sort( pSeeds->begin( ), pSeeds->end( ),
               [&]( const Seed& a, const Seed& b ) {
                   return getPositionForBucketing( uiQLen, a ) < getPositionForBucketing( uiQLen, b );
               } // lambda
    ); // sort function call
#endif


#if MEASURE_DURATIONS == ( 1 )
                           } // lambda
                                                )
                               .count( ) *
                           1000; // metaMeasureDuration function call
#endif

    // positions to remember the maxima
    auto pSoCs = std::make_shared<SoCPriorityQueue>( pSeeds );

#if MEASURE_DURATIONS == ( 1 )
    pSegments->fLinesweep +=
        metaMeasureDuration( [&]( ) {
#endif
            // find the SOC maxima
            SoCOrder xCurrScore;
            std::vector<Seed>::iterator xStripStart = pSeeds->begin( );
            std::vector<Seed>::iterator xStripEnd = pSeeds->begin( );
            while( xStripEnd != pSeeds->end( ) && xStripStart != pSeeds->end( ) )
            {
                // move xStripEnd forwards while it is closer to xStripStart than uiStripSize
                while( xStripEnd != pSeeds->end( ) &&
                       getPositionForBucketing( uiQLen, *xStripStart ) + uiStripSize >=
                           getPositionForBucketing( uiQLen, *xStripEnd )

#if CONTIG_ID_CACHE == ( 1 )
                       && xStripStart->uiContigId == xStripEnd->uiContigId
#endif
                )
                {
                    // remember the additional score
                    xCurrScore += *xStripEnd;
                    // move the iterator forward
                    xStripEnd++;
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


#if MEASURE_DURATIONS == ( 1 )
        } // lambda
                             )
            .count( ) *
        1000; // metaMeasureDuration function call
#endif

    // return the strip collection
    return pSoCs;

} // function

std::shared_ptr<SoCPriorityQueue> MinimizerStripOfConsideration::execute(
    std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuerySeq, std::shared_ptr<Pack> pPack )
{
    // make sure that we return at least an SoC set
    if( pSeeds->empty( ) )
        return std::make_shared<SoCPriorityQueue>( pSeeds );

    auto uiQLen = pQuerySeq->length( );
    double fMinLen = std::max( (double)fGiveUp * uiQLen, (double)this->uiCurrHarmScoreMin );

    /*
     * This is the formula from the paper
     * computes the size required for the strip so that we collect all relevent seeds.
     */
    const nucSeqIndex uiStripSize = this->getStripSize( uiQLen, iMatch, iExtend, iGap );

#if 0 || DELTA_CACHE == ( 0 )
    // sort the seeds according to their initial positions
    std::sort( pSeeds->begin( ), pSeeds->end( ),
               [&]( const Seed& rA, const Seed& rB ) {
                   if( rA.start_ref( ) - (int64_t)rA.start( ) != rB.start_ref( ) - (int64_t)rB.start( ) )
                       return rA.start_ref( ) - (int64_t)rA.start( ) < rB.start_ref( ) - (int64_t)rB.start( );
                   return rA.start( ) < rB.start( );
               } // lambda
    ); // sort function call
#else
    for( auto& rS : *pSeeds )
        rS.uiDelta = getPositionForBucketing( uiQLen, rS );
    radix_sort<8>( &pSeeds->front( ), &pSeeds->front( ) + pSeeds->size( ),
                   []( const Seed& rSeed ) { return rSeed.uiDelta; } );
#endif

    // combine seeds
    auto pSeeds2 = xLumper.execute( pSeeds );
    auto pSeeds3 = xExtender.execute( pSeeds2, pQuerySeq, pPack );

    // positions to remember the maxima
    auto pSoCs = std::make_shared<SoCPriorityQueue>( pSeeds3 );

    // find the SOC maxima
    SoCOrder xCurrScore;
    std::vector<Seed>::iterator xStripStart = pSeeds3->begin( );
    std::vector<Seed>::iterator xStripEnd = pSeeds3->begin( );
    while( xStripEnd != pSeeds3->end( ) && xStripStart != pSeeds3->end( ) )
    {
        // move xStripEnd forwards while it is closer to xStripStart than uiStripSize
        while( xStripEnd != pSeeds3->end( ) && getPositionForBucketing( uiQLen, *xStripStart ) + uiStripSize >=
                                                   getPositionForBucketing( uiQLen, *xStripEnd ) )
        {
            // remember the additional score
            xCurrScore += *xStripEnd;
            // move the iterator forward
            xStripEnd++;
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

    // return the strip collection
    return pSoCs;
} // method

#ifdef WITH_PYTHON

void exportStripOfConsideration( py::module& rxPyModuleId )
{
    // export the StripOfConsideration class
    exportModule<StripOfConsideration>( rxPyModuleId, "StripOfConsideration" );
    exportModule<MinimizerStripOfConsideration>( rxPyModuleId, "MinimizerStripOfConsideration" );
} // function

#endif
