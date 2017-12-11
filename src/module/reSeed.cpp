
#include "module/reSeed.h"


#include <vector>
#include <memory>
#include <atomic>
#include <chrono>
using namespace libLAuS;
 
#define complement(x) (uint8_t)NucSeq::nucleotideComplement(x)

/**
* Delivers 4 intervals for a single input interval.
* Here we use only two fields in the BWT_Interval.
*/
SAInterval ReSeed::extend_backward( 
        // current interval
        const SAInterval &ik,
        // the character to extend with
        const uint8_t c,
        std::shared_ptr<FMIndex> pFM_index
    )
{
    bwt64bitCounter cntk[4]; // Number of A, C, G, T in BWT until start of interval ik
    bwt64bitCounter cntl[4]; // Number of A, C, G, T in BWT until end of interval ik

    assert(ik.start() < ik.end());
    assert(ik.start() > 0);

    //here the intervals seem to be (a,b] while mine are [a,b)
    pFM_index->bwt_2occ4(
        // start of SA index interval
        ik.start() - 1,
        // end of SA index interval
        ik.end() - 1,
        cntk,                        // output: Number of A, C, G, T until start of interval
        cntl                        // output: Number of A, C, G, T until end of interval
    );

    for(unsigned int i = 0; i < 4; i++)
        assert(cntk[i] <= cntl[i]);

    bwt64bitCounter cnts[4]; // Number of A, C, G, T in BWT interval ik
    //the cnts calculated here might be off by one
    for(unsigned int i = 0; i < 4; i++)
        cnts[i] = cntl[i] - cntk[i];

    DEBUG_2(
        std::cout << cnts[0] << " + " << cnts[1] << " + " << cnts[2] << " + " << cnts[3] << " = " 
                  << (t_bwtIndex)(cnts[0] + cnts[1] + cnts[2] + cnts[3]) << " ?= " 
                  << ik.size() << "(-1)" << std::endl;
    )

    bwt64bitCounter cntk_2[4];
    cntk_2[0] = ik.startRevComp();
    /*
     * PROBLEM:
     * 
     * the representation of the $ in the count part of the FM_index is indirect
     *         done by storing the position of the $
     * if have two bwt indices k and l
     * the counts do not return the $ obviously...
     * 
     * The result may be off by one since sometimes we have a $ before the current pos 
     * sometimes we do not...
     *
     * lets adjust the sizes of the smaller intervals accordingly
     */
    if(
            ik.start() < pFM_index->primary && 
            ik.end() >= pFM_index->primary
        )
    {
        cntk_2[0]++;
        DEBUG_2(
            std::cout << "adjusted cntk_2[0] because of primary" << std::endl;
        )
        assert( (t_bwtIndex)(cnts[0] + cnts[1] + cnts[2] + cnts[3]) == ik.size() - 1 );
    }//if
    else{
        if( (t_bwtIndex)(cnts[0] + cnts[1] + cnts[2] + cnts[3]) != ik.size() )
        {
            std::cout << ik.start() << " " << ik.end() << " " << pFM_index->primary << std::endl;
            std::cout << cnts[0] << " + " << cnts[1] << " + " << cnts[2] << " + " <<
                cnts[3] << " = " << (t_bwtIndex)(cnts[0] + cnts[1] + cnts[2] + cnts[3]) << " ?= "
                << ik.size() << "(-1)" << std::endl;
        }//if
        assert( (t_bwtIndex)(cnts[0] + cnts[1] + cnts[2] + cnts[3]) == ik.size() );
    }//else
    //for all nucleotides
    for(unsigned int i = 1; i < 4; i++)
        cntk_2[i] = cntk_2[i-1] + cnts[complement(i-1)];



    //BWAs SA intervals seem to be (a,b] while mine are [a,b)
    //pFM_index->L2[c] start of nuc c in BWT
    //cntk[c] offset of new interval
    //cntl[c] end of new interval
    return SAInterval(pFM_index->L2[c] + cntk[c] + 1, cntk_2[complement(c)], cnts[c]);
} // method


void ReSeed::extend(
        std::shared_ptr<SegmentVector> pxVector,
        nucSeqIndex min,
        nucSeqIndex max,
        std::shared_ptr<FMIndex> pFM_index,
        std::shared_ptr<NucSeq> pQuerySeq
    )
{

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
                        pFM_index->L2[complement(q[min])] + 1, 
                        pFM_index->L2[(int)q[min]] + 1, 
                        pFM_index->L2[(int)q[min] + 1] - pFM_index->L2[(int)q[min]]
                    );

    std::list<Segment> curr = std::list<Segment>();
    for(nucSeqIndex i = min+1; i < max; i++)
    {
        DEBUG_2(
            std::cout << i-1 << " -> " << ik.start() << " " << ik.end() << std::endl;
            std::cout << i-1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
        )
        assert(ik.size() > 0);
        SAInterval ok = extend_backward(ik, complement(q[i]), pFM_index);

        DEBUG_2(
            std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
            std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
        )
        /*
        * In fact, if ok.getSize is zero, then there are no matches any more.
        */
        if (ok.size() == 0)
            max = i-1;
            break; // the SA-index interval size is too small to be extended further
        ik = ok;
    }//for
    std::shared_ptr<Segment> pSeg(new Segment(min,max-min,ik.revComp()));
    assert(min >= 0);
    assert(max < pQuerySeq->length());
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
        ContainerVector vpInput
    )
{
    std::shared_ptr<FMIndex> pFM_index = std::static_pointer_cast<FMIndex>(vpInput[0]);
    std::shared_ptr<SegmentVector> pSegments = std::static_pointer_cast<SegmentVector>(vpInput[1]);
    std::shared_ptr<NucSeq> pQuerySeq = 
        std::static_pointer_cast<NucSeq>(vpInput[2]);

    std::shared_ptr<SegmentVector> pSegmentVector(new SegmentVector());

    for(std::shared_ptr<Segment> pxNode : *pSegments)
    {
        pSegmentVector->push_back(pxNode);

        if(pxNode->size() > 2)
            extend(pSegmentVector, pxNode->start() + 1, pxNode->end() - 1, pFM_index, pQuerySeq);
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