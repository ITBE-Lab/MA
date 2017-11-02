/**
 * @package LAuS.sweepAllReturnBest
 * @brief Implements @ref LAuS.sweepAllReturnBest.SweepAllReturnBest "SweepAllReturnBest".
 * @file sweepAllReturnBest.py
 * @brief Implements @ref LAuS.sweepAllReturnBest.SweepAllReturnBest "SweepAllReturnBest".
 * @author Markus Schmidt
*/

#ifndef SWEEP_ALL_RETURN_BEST
#define SWEEP_ALL_RETURN_BEST

#include "linesweep.h"
#include "threadPool.h"

/**
 * @brief Execute LineSweep for all given seeds.
 * @details 
 * The module calls LineSweep for all Seeds in the SeedsVector.
 * @ingroup module
 */
class SweepAllReturnBest: public CppModule
{
private:
    static LineSweep xLineSweep;

public:
    SweepAllReturnBest()
            :
        CppModule()
    {}//default constructor

    //overload
    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput)
    {
        std::shared_ptr<SeedsVector> pSeedsVector = std::shared_ptr<SeedsVector>(
            std::static_pointer_cast<SeedsVector>(vpInput[0]));
        std::shared_ptr<NucleotideSequence> pQuerySeq =
            std::static_pointer_cast<NucleotideSequence>(vpInput[1]);
        std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq =
            std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[2]);

        SeedsVector vTempResults = SeedsVector(pSeedsVector->size());
        {
            ThreadPool xPool(8);
            for(unsigned int i = 0; i < vTempResults.size(); i++)
            {
                vTempResults[i] = std::shared_ptr<Seeds>();
                xPool.enqueue(
                    []
                    (
                        size_t, 
                        std::shared_ptr<SeedsVector> pSeedsVector,
                        std::shared_ptr<NucleotideSequence> pQuerySeq,
                        std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq,
                        SeedsVector* vTempResults,
                        LineSweep& rLineSweep,
                        unsigned int i
                    )
                    {
                        std::vector<std::shared_ptr<Container>> vInput = {
                                pQuerySeq, 
                                pRefSeq, 
                                (*pSeedsVector)[i]
                            };
                        try
                        {
                            (*vTempResults)[i] = std::static_pointer_cast<Seeds>(
                                    rLineSweep.execute(vInput));
                        }
                        catch(NullPointerException e) 
                        {
                            std::cerr << e.what() << std::endl;
                        }
                        catch(std::exception e) 
                        {
                            std::cerr << e.what() << std::endl;
                        }
                        catch(...) 
                        {
                            std::cerr << "unknown exception when executing" << std::endl;
                        }
                        if((*vTempResults)[i] == nullptr)
                            throw NullPointerException("linesweep deleviered nullpointer as result");
                    },//lambda
                    pSeedsVector, pQuerySeq, pRefSeq, &vTempResults, xLineSweep, i
                );
            }//for
        }//scope xPool

        if(vTempResults.empty())
            return std::shared_ptr<Seeds>(new Seeds());

        std::shared_ptr<Seeds> pRet = vTempResults[0];
        assert(pRet != nullptr);
        for(std::shared_ptr<Seeds> pSeeds : vTempResults)
        {
            assert(pSeeds != nullptr);
            if(pSeeds->getScore() > pRet->getScore())
                pRet = pSeeds;
            }//for

        return pRet;
    }//function

    //overload
    std::vector<ContainerType> getInputType()
    {
        return std::vector<ContainerType>{
                ContainerType::seedsVector,
                ContainerType::nucSeq,
                ContainerType::packedNucSeq
            };
    }//function

    //overload
    ContainerType getOutputType()
    {
        return ContainerType::seeds;
    }//function
};//class

/**
 * @brief Exposes the SweepAllReturnBest @ref CppModule "module" to boost python.
 * @ingroup export
 */
void exportSweepAll();

#endif
