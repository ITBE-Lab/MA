/** 
 * @file reSeed.cpp
 * @author Markus Schmidt
 */
#include "module/reSeed.h"


#include <vector>
#include <memory>
#include <atomic>
#include <chrono>
using namespace libMA;
 
#define complement(x) (uint8_t)NucSeq::nucleotideComplement(x)

void ReSeed::extend(
        std::shared_ptr<SegmentVector> pxVector,
        nucSeqIndex center,
        std::shared_ptr<FMIndex> pFM_index,
        std::shared_ptr<NucSeq> pQuerySeq
    )
{
    if(center == 0)
        return;
    if(center == pQuerySeq->length()-1)
        return;
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
        pFM_index->L2[(int)q[center]] + 1,
        pFM_index->L2[complement(q[center])] + 1,
        pFM_index->L2[(int)q[center] + 1] - pFM_index->L2[(int)q[center]]
    );

    std::list<Segment> curr = std::list<Segment>();
    bool bBackwards = true;
    nucSeqIndex lower = center - 1;
    nucSeqIndex higher = center + 1;
    /*
     * extend the segment on both ends until the ambiguity is small enough for us to use it.
     */
    while(ik.size() >= maxAmbiguity)
    {
        //unsigned so check for underflows of lower...
        if(lower > higher && bBackwards)
            return;
        //make sure we dont extend past the query length
        if(higher >= pQuerySeq->length() && !bBackwards)
            return;

        DEBUG_2(
            std::cout << i-1 << " -> " << ik.start() << " " << ik.end() << std::endl;
            std::cout << i-1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
        )

        assert(ik.size() > 0);
        /*
         * the extends forwads and backwards in turns,
         * while always reversing the resulting SA_Interval
         */
        ik = pFM_index->extend_backward(
            ik, 
            /*
             * if we extend backwards we just need to use the nucleotide at the correct position
             * otherwise we need the complement of the nucleotide at position "higher"
             */
            bBackwards ? q[lower--] : complement(q[higher++])
        ).revComp();

        bBackwards = !bBackwards;

        DEBUG_2(
            std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
            std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
        )
    }//for

    /*
    * make sure we always get the forward interval as the main one
    */
    if(!bBackwards)
        ik = ik.revComp();
    /*
        * save the correct segment
        * second parameter is the length but lower and higher always point to the elements
        * one past the last extended one therefore -2
        */
    std::shared_ptr<Segment> pSeg(new Segment(lower+1,higher-lower-2,ik));
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

        extend(pSegmentVector, pxNode->center(), pFM_index, pQuerySeq);
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
        >("ReSeed", boost::python::init<int>())
        .def_readwrite("max_ambiguity", &ReSeed::maxAmbiguity)
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<ReSeed>,
        std::shared_ptr<Module> 
    >();
}//function