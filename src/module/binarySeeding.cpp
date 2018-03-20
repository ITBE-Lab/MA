/** 
 * @file binarySeeding.cpp
 * @author Markus Schmidt
 */
#include "module/binarySeeding.h"

using namespace libMA;


#include <vector>
#include <memory>
#include <atomic>
#include <chrono>
 

/* this function implements the segmentation of the query
 *
 * the process is synchronized using a thread pool
 *
 * to avoid unnecessary locking and unlocking we make sure that only one thread can process one interval at any given time.
 * ( there are always several intervals ready and waiting to be processed. )
 * this way it is only necessary to lock once we touch the structure of the list which stores the individual intervals.
 *
 * segmentation technique:
 *        start in the middle of the interval try to extend in both directions
 *        split the interval in 3 parts: prev:perfectMatch:post
 *        queue the prev and post intervals into the thread pool
 *        save the perfect match for later clustering
*/
void BinarySeeding::procesInterval(
            Interval<nucSeqIndex> xAreaToCover,
            std::shared_ptr<SegmentVector> pSegmentVector,
            std::shared_ptr<FMIndex> pFM_index,
            std::shared_ptr<NucSeq> pQuerySeq
        )
{
    while(true)
    {
        DEBUG_2(
                std::cout << "interval (" << xAreaToCover.start() << "," << xAreaToCover.end() << ")" << std::endl;
            )

        Interval<nucSeqIndex> xAreaCovered;
        // performs extension and records any found seeds
        // here we use bLrExtension to choose the extension scheme
        if(bLrExtension)
            xAreaCovered = maximallySpanningExtension(
                xAreaToCover.center(), 
                pFM_index, 
                pQuerySeq, 
                pSegmentVector
            );
        else
            xAreaCovered = smemExtension(
                xAreaToCover.center(), 
                pFM_index,
                pQuerySeq, 
                pSegmentVector
            );

        DEBUG_2(
            std::cout << "splitting interval (" << xAreaToCover.start() << "," << xAreaToCover.end() << ") at (" << xAreaCovered.start() << "," << xAreaCovered.end() << ")" << std::endl;
        )

        // if the extension did not fully cover until uiStart:
        if (xAreaCovered.start() != 0 && xAreaToCover.start() + 1 < xAreaCovered.start())
        {
            // enqueue procesInterval() for a new interval that spans from uiStart to 
            // where the extension stopped
            procesInterval(
                Interval<nucSeqIndex>(
                    xAreaToCover.start(), xAreaCovered.start() - xAreaToCover.start() - 1
                ),
                pSegmentVector,
                pFM_index,
                pQuerySeq
            );
        }//if
        // if the extension did not fully cover until uiEnd:
        if (xAreaToCover.end() > xAreaCovered.end() + 1)
        {
            xAreaToCover.set(xAreaCovered.end() + 1, xAreaToCover.end() - xAreaCovered.end() - 1);
            // REPLACED by while loop
            // enqueue procesInterval() for a new interval that spans from where the extension stopped
            // to uiEnd
            //procesInterval(
            //    Interval<nucSeqIndex>(uiTo + 1, xAreaToCover.end() - uiTo - 1)),
            //    pSegmentVector,
            //    pFM_index,
            //    pQuerySeq
            //);
        }//if
        else
            break;
    }//while
}//function


ContainerVector BinarySeeding::getInputType() const
{
    return ContainerVector{
            //the forward fm_index
            std::shared_ptr<Container>(new FMIndex()),
            //the query sequence
            std::shared_ptr<Container>(new NucSeq()),
        };
}
std::shared_ptr<Container> BinarySeeding::getOutputType() const
{
    return std::shared_ptr<Container>(new SegmentVector());
}


std::shared_ptr<Container> BinarySeeding::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<FMIndex> pFM_index = std::static_pointer_cast<FMIndex>((*vpInput)[0]);
    std::shared_ptr<NucSeq> pQuerySeq = 
        std::dynamic_pointer_cast<NucSeq>((*vpInput)[1]);

    std::shared_ptr<SegmentVector> pSegmentVector(new SegmentVector());
    if(pQuerySeq == nullptr)
        return pSegmentVector;

    DEBUG_2(
        std::cout << pQuerySeq->fastaq() << std::endl;
    )

    procesInterval(
        Interval<nucSeqIndex>(0, pQuerySeq->length()),
        pSegmentVector,
        pFM_index,
        pQuerySeq
    );

    return pSegmentVector;
}//function

void exportBinarySeeding()
{
    //export the BinarySeeding class
    boost::python::class_<
            BinarySeeding, 
            boost::python::bases<Module>,
            std::shared_ptr<BinarySeeding>
        >(
            "BinarySeeding",
            boost::python::init<bool>()
        )
        ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<BinarySeeding>,
        std::shared_ptr<Module> 
    >();

}//function