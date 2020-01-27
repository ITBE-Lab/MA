/**
 * @file binarySeeding.cpp
 * @author Markus Schmidt
 */
#include "module/binarySeeding.h"

using namespace libMA;


#include <atomic>
#include <chrono>
#include <memory>
#include <vector>


/* this function implements the segmentation of the query
 *
 * the process is synchronized using a thread pool
 *
 * to avoid unnecessary locking and unlocking we make sure that only one thread can process one
 * interval at any given time. ( there are always several intervals ready and waiting to be
 * processed. ) this way it is only necessary to lock once we touch the structure of the list which
 * stores the individual intervals.
 *
 * segmentation technique:
 *        start in the middle of the interval try to extend in both directions
 *        split the interval in 3 parts: prev:perfectMatch:post
 *        queue the prev and post intervals into the thread pool
 *        save the perfect match for later clustering
 */
void BinarySeeding::procesInterval( Interval<nucSeqIndex> xAreaToCover,
                                    std::shared_ptr<SegmentVector>
                                        pSegmentVector,
                                    std::shared_ptr<SuffixArrayInterface>
                                        pFM_index,
                                    std::shared_ptr<NucSeq>
                                        pQuerySeq,
                                    size_t uiCnt )
{
    while( true )
    {
        DEBUG_2( std::cout << "interval (" << xAreaToCover.start( ) << "," << xAreaToCover.end( ) << ")" << std::endl; )

        Interval<nucSeqIndex> xAreaCovered;
        // performs extension and records any found seeds
        // here we use bLrExtension to choose the extension scheme
        if( bLrExtension )
            xAreaCovered =
                maximallySpanningExtension( xAreaToCover.center( ), pFM_index, pQuerySeq, pSegmentVector, uiCnt == 0 );
        else
            xAreaCovered = smemExtension( xAreaToCover.center( ), pFM_index, pQuerySeq, pSegmentVector );

        DEBUG_2( std::cout << "splitting interval (" << xAreaToCover.start( ) << "," << xAreaToCover.end( ) << ") at ("
                           << xAreaCovered.start( ) << "," << xAreaCovered.end( ) << ")" << std::endl; )

        // if the extension did not fully cover until uiStart:
        if( xAreaCovered.start( ) != 0 && xAreaToCover.start( ) + 1 < xAreaCovered.start( ) )
        {
            // enqueue procesInterval() for a new interval that spans from uiStart to
            // where the extension stopped
            procesInterval(
                Interval<nucSeqIndex>( xAreaToCover.start( ), xAreaCovered.start( ) - xAreaToCover.start( ) ),
                pSegmentVector, pFM_index, pQuerySeq, uiCnt + 1 );
        } // if
        // if the extension did not fully cover until uiEnd:
        if( xAreaToCover.end( ) > xAreaCovered.end( ) + 1 )
        {
            xAreaToCover.set( xAreaCovered.end( ), xAreaToCover.end( ) - xAreaCovered.end( ) );
            // REPLACED by while loop
            // enqueue procesInterval() for a new interval that spans from where the extension
            // stopped to uiEnd
            // procesInterval(
            //    Interval<nucSeqIndex>(uiTo + 1, xAreaToCover.end() - uiTo - 1)),
            //    pSegmentVector,
            //    pFM_index,
            //    pQuerySeq
            //);
        } // if
        else
            break;
        uiCnt++;
    } // while
} // function

std::shared_ptr<SegmentVector> BinarySeeding::execute( std::shared_ptr<SuffixArrayInterface> pFM_index,
                                                       std::shared_ptr<NucSeq> pQuerySeq )
{
    std::shared_ptr<SegmentVector> pSegmentVector( new SegmentVector( ) );

    if( pQuerySeq == nullptr )
        return pSegmentVector;
    if( pQuerySeq->length() == 0 )
        return pSegmentVector;

    DEBUG_2( std::cout << pQuerySeq->fastaq( ) << std::endl; )

    procesInterval( Interval<nucSeqIndex>( 0, pQuerySeq->length( ) ), pSegmentVector, pFM_index, pQuerySeq, 0 );

    /*
     * If we have extremeley few seeds we may not be able to compute more than one alignment.
     * In that case we cannot do a mapping quality estimation. Therefore we will do a reseeding in this case.
     */
#if 0
    if(pSegmentVector->size() <= 10)
    {
#if 0
        size_t uiNumSubsegments = 4;
        size_t uiSegVecSize = pSegmentVector->size();
        uint8_t* puiActualQueryStartPointer = pQuerySeq->pxSequenceRef;
        size_t uiActualQuerySize = pQuerySeq->uiSize;
        for(size_t i = 0; i < uiNumSubsegments; i++)
        {
            size_t uiSubSegStart = i * uiActualQuerySize / uiNumSubsegments;
            size_t uiSubSegEnd = (i+1) * uiActualQuerySize / uiNumSubsegments;
            pQuerySeq->uiSize = uiSubSegEnd - uiSubSegStart;
            pQuerySeq->pxSequenceRef = puiActualQueryStartPointer + uiSubSegStart;
            procesInterval( Interval<nucSeqIndex>( 0, pQuerySeq->length( ) ), pSegmentVector, pFM_index, pQuerySeq );
            // move the query start positions of the computet segments forward
            while(uiSegVecSize < pSegmentVector->size())
            {
                (*pSegmentVector)[uiSegVecSize].iStart += uiSubSegStart;
                uiSegVecSize++;
            } // while
        } // for
#else
        size_t uiSizeSubsegments = 30;
        size_t uiStepSubsegments = 10;
        size_t uiSegVecSize = pSegmentVector->size();
        uint8_t* puiActualQueryStartPointer = pQuerySeq->pxSequenceRef;
        size_t uiActualQuerySize = pQuerySeq->uiSize;
        for(size_t i = 0; i < uiActualQuerySize - uiSizeSubsegments; i+=uiStepSubsegments)
        {
            size_t uiSubSegStart = i;
            size_t uiSubSegEnd = i + uiSizeSubsegments;
            pQuerySeq->uiSize = uiSubSegEnd - uiSubSegStart;
            pQuerySeq->pxSequenceRef = puiActualQueryStartPointer + uiSubSegStart;
            procesInterval( Interval<nucSeqIndex>( 0, pQuerySeq->length( ) ), pSegmentVector, pFM_index, pQuerySeq );
            // move the query start positions of the computet segments forward
            while(uiSegVecSize < pSegmentVector->size())
            {
                (*pSegmentVector)[uiSegVecSize].iStart += uiSubSegStart;
                uiSegVecSize++;
            } // while
        } // for
#endif
        pQuerySeq->pxSequenceRef = puiActualQueryStartPointer;
        pQuerySeq->uiSize = uiActualQuerySize;
    } // if
#endif

#if 0
    // for the illumina alignment
    size_t uiLongest = 0;
    for(size_t i = 1; i < pSegmentVector->size(); i++)
        if((*pSegmentVector)[i].size() > (*pSegmentVector)[uiLongest].size())
            uiLongest = i;
    if((*pSegmentVector)[uiLongest].size() >= 20 && (*pSegmentVector)[uiLongest].saInterval().size() > 50)
        pSegmentVector->bSetMappingQualityToZero = true;
#endif


    /*
     * Observation:
     * If the minimum seed size is below uiMinSeedSizeDrop we can abort here already,
     * without loosing accuracy
     */
    if( !bDisableHeuristics && uiMinSeedSizeDrop != 0 &&
        pSegmentVector->numSeedsLarger( uiMinSeedSizeDrop ) < fRelMinSeedSizeAmount * pQuerySeq->length( ) &&
        uiMinGenomeSize < pFM_index->getRefSeqLength( ) )
        pSegmentVector->clear( );

    return pSegmentVector;
} // function

#ifdef WITH_PYTHON

void exportBinarySeeding( py::module& rxPyModuleId )
{
    // export the BinarySeeding class
    exportModule<BinarySeeding>( rxPyModuleId, "BinarySeeding",
                                 []( auto&& x ) { x.def( "seed", &BinarySeeding::seed ); } );
} // function
#endif