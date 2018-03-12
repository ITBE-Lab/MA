/** 
 * @file binarySeeding.cpp
 * @author Markus Schmidt
 */
#include "module/otherSeeding.h"

using namespace libMA;

#define INCLUDE_SELF_TESTS (  1 )

//#include "BWTmem.h"
#include <vector>
//#include "analysis.h"
#include <memory>
#include <atomic>
#include <chrono>
//#include "assembly.h"
 
#define complement(x) (uint8_t)NucSeq::nucleotideComplement(x)


/// @todo  i should be my own module
void OtherSeeding::bowtieExtension(
        std::shared_ptr<FMIndex> pFM_index,
        std::shared_ptr<NucSeq> pQuerySeq,
        std::shared_ptr<SegmentVector> pSegmentVector
    )
{
    unsigned int iSize = 16;
    unsigned int iStep = 1;

    const uint8_t *q = pQuerySeq->pGetSequenceRef(); 
    for(nucSeqIndex i = 0; i < pQuerySeq->length()-iSize; i+= iStep)
    {

        SAInterval ik(
                            pFM_index->L2[complement(q[i])] + 1, 
                            pFM_index->L2[(int)q[i]] + 1, 
                            pFM_index->L2[(int)q[i] + 1] - pFM_index->L2[(int)q[i]]
                        );
        for(nucSeqIndex i2 = 1; i2 <= iSize; i2++)
        {
            // this is the extension
            // actually forward since we extend backwards on the reverse complement...
            ik = pFM_index->extend_backward(ik, complement(q[i2 + i]));
            if(ik.size() == 0)
                break;
        }//for
        if(ik.size() == 0)
            continue;
        // found a seed
        pSegmentVector->push_back(std::shared_ptr<Segment>(new Segment(i, iSize, ik.revComp())));
    }//for
}//function

/// @todo  i should be my own module
/**
 * This is what blasr does:
 * 
 * for each position in the query:
 *  extend maximally backwards 
 *  save a seed that is one shorter than the maximal extension
 *  if the seed is longer than K = 12 adn maxambiguity = ?
 */
void OtherSeeding::doBlasrExtension(
        std::shared_ptr<FMIndex> pFM_index,
        std::shared_ptr<NucSeq> pQuerySeq,
        std::shared_ptr<SegmentVector> pSegmentVector
    )
{
    const uint8_t *q = pQuerySeq->pGetSequenceRef(); 
    for(nucSeqIndex i = 0; i < pQuerySeq->length(); i+= 1)
    {

        SAInterval ik(
                            pFM_index->L2[(int)q[i]] + 1, 
                            pFM_index->L2[complement(q[i])] + 1, 
                            pFM_index->L2[(int)q[i] + 1] - pFM_index->L2[(int)q[i]]
                        );
        SAInterval lk;// will hold the maximal extension
        SAInterval llk;// will hold the interval one shorter than the maximal extension
        nucSeqIndex i2 = 0;
        while(i2 <= i)
        {
            llk = lk;
            lk = ik;
            ik = pFM_index->extend_backward(ik, q[i - i2]);
            if(ik.size() == 0)
                break;
            i2++;
        }//for
        if(i2 <= 12)
            continue;
        // found a seed
        pSegmentVector->push_back(std::shared_ptr<Segment>(new Segment(i - i2 + 1, i2-1, llk)));
    }//for
}//function


ContainerVector OtherSeeding::getInputType() const
{
    return ContainerVector{
            //the forward fm_index
            std::shared_ptr<Container>(new FMIndex()),
            //the query sequence
            std::shared_ptr<Container>(new NucSeq()),
        };
}
std::shared_ptr<Container> OtherSeeding::getOutputType() const
{
    return std::shared_ptr<Container>(new SegmentVector());
}


std::shared_ptr<Container> OtherSeeding::execute(
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

    if(bBowtie)
        bowtieExtension(pFM_index, pQuerySeq, pSegmentVector);
    else
        doBlasrExtension(pFM_index, pQuerySeq, pSegmentVector);

    return pSegmentVector;
}//function

void exportOtherSeeding()
{
    //export the OtherSeeding class
    boost::python::class_<
            OtherSeeding, 
            boost::python::bases<Module>,
            std::shared_ptr<OtherSeeding>
        >(
            "OtherSeeding",
            boost::python::init<bool>()
        )
        ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<OtherSeeding>,
        std::shared_ptr<Module> 
    >();

}//function