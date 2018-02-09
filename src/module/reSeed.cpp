
#include "module/reSeed.h"


#include <vector>
#include <memory>
#include <atomic>
#include <chrono>
using namespace libMA;
 
#define complement(x) (uint8_t)NucSeq::nucleotideComplement(x)

void ReSeed::extend(
        std::shared_ptr<SegmentVector> pxVector,
        nucSeqIndex min,
        nucSeqIndex max,
        std::shared_ptr<FMIndex> pFM_index,
        std::shared_ptr<NucSeq> pQuerySeq
    )
{
    assert(min >= 0);
    assert(max < pQuerySeq->length());
    assert(max > 0);
    // query sequence itself
    const uint8_t *q = pQuerySeq->pGetSequenceRef(); 
    
    /* Initialize ik on the foundation of the single base q[x].
     * In order to understand this initialization you should have a look 
     *to the corresponding PowerPoint slide.
     */
    // start I(q[x]) in T (start in BWT used for backward search) + 1, 
    // because very first string in SA-array starts with $
    // size in T and T' is equal due to symmetry
    SAInterval ik(
                        pFM_index->L2[(int)q[max]] + 1, 
                        pFM_index->L2[complement(q[max])] + 1, 
                        pFM_index->L2[(int)q[max] + 1] - pFM_index->L2[(int)q[max]]
                    );

    std::list<Segment> curr = std::list<Segment>();
    for(nucSeqIndex i = max-1; i >= min; i--)
    {
        DEBUG_2(
            std::cout << i-1 << " -> " << ik.start() << " " << ik.end() << std::endl;
            std::cout << i-1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
        )
        assert(ik.size() > 0);
        ik = pFM_index->extend_backward(ik, q[i]);

        DEBUG_2(
            std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
            std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
        )

        //unsigned -> so prevent underflows
        if(i == 0)
            break;
    }//for
    std::shared_ptr<Segment> pSeg(new Segment(min,max-min,ik));
    assert(pSeg->end() < pQuerySeq->length());
    pxVector->push_back(pSeg);
}//function


ContainerVector ReSeed::getInputType() const
{
    return ContainerVector{
            //the forward fm_index
            std::shared_ptr<Container>(new FMIndex()),
            //the forward fm_index
            std::shared_ptr<Container>(new SegmentVector()),
            //the query sequence
            std::shared_ptr<Container>(new NucSeq()),
        };
}

std::shared_ptr<Container> ReSeed::getOutputType() const
{
    return std::shared_ptr<Container>(new SegmentVector());
}


std::shared_ptr<Container> ReSeed::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<FMIndex> pFM_index = std::static_pointer_cast<FMIndex>((*vpInput)[0]);
    std::shared_ptr<SegmentVector> pSegments = std::static_pointer_cast<SegmentVector>((*vpInput)[1]);
    std::shared_ptr<NucSeq> pQuerySeq = 
        std::static_pointer_cast<NucSeq>((*vpInput)[2]);

    std::shared_ptr<SegmentVector> pSegmentVector(new SegmentVector());

    for(std::shared_ptr<Segment> pxNode : *pSegments)
    {
        pSegmentVector->push_back(pxNode);

        if(pxNode->size() > 2)
            extend(pSegmentVector, pxNode->start() + 1, pxNode->end() - 2, pFM_index, pQuerySeq);
    }//for

    return pSegmentVector;
}//function

void exportReSeed()
{
    //export the ReSeed class
    boost::python::class_<
            ReSeed,
            boost::python::bases<Module>,
            std::shared_ptr<ReSeed>
        >("ReSeed")
        .def_readwrite("min_split_len", &ReSeed::minSplitLen)
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<ReSeed>,
        std::shared_ptr<Module> 
    >();
}//function