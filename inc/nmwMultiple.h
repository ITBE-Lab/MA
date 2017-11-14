//TODO: docu; make params acessible



#ifndef NMW_MULTIPLE_H
#define NMW_MULTIPLE_H

#include "needlemanWunsch.h"

/**
 * @brief Execute LineSweep for all given seeds.
 * @details 
 * The module calls LineSweep for all Seeds in the SeedsVector.
 * @ingroup module
 */
class NmwMultiple: public CppModule
{
    NeedlemanWunsch nmw;
public:
    unsigned int tryNmany = 10;

    //overload
    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput)
    {
        std::shared_ptr<SeedsVector> pSeedsVector = std::shared_ptr<SeedsVector>(
            std::static_pointer_cast<SeedsVector>(vpInput[0]));
        std::shared_ptr<NucleotideSequence> pQuery 
            = std::static_pointer_cast<NucleotideSequence>(vpInput[1]);
        std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefPack = 
            std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[2]);

        std::shared_ptr<Alignment> pRet = std::shared_ptr<Alignment>(new Alignment(0, 0));
        for(unsigned int i = 0; i < tryNmany && i < pSeedsVector->size(); i++)
        {
            std::vector<std::shared_ptr<Container>> vInput = {
                        (*pSeedsVector)[i],
                        pQuery,
                        pRefPack
                    };//vector
            std::shared_ptr<Alignment> pCurr = std::static_pointer_cast<Alignment>(nmw.execute(vInput));
            if(pCurr->score(10, -3, -4, -1) > pRet->score(10, -3, -4, -1))
                pRet = pCurr;
        }//for

        return pRet;
    }//function

    //overload
    std::vector<ContainerType> getInputType()
    {
        return std::vector<ContainerType>{
                ContainerType::seedsVector,
                //the query sequence
                ContainerType:nucSeq,
                //the reference sequence
                ContainerType:packedNucSeq,
            };
    }//function

    //overload
    ContainerType getOutputType()
    {
        return ContainerType::alignment;
    }//function
};//class

/**
 * @brief Exposes the SweepAllReturnBest @ref CppModule "module" to boost python.
 * @ingroup export
 */
void exportNmwMultiple();

#endif
