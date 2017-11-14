/**
 * TODO:
 * @file execOnVector.h
 */

#ifndef EXEC_ON_VECTOR_H
#define EXEC_ON_VECTOR_H

#include "cppModule.h"
#include "threadPool.h"

/**
 * @brief Execute LineSweep for all given seeds.
 * @details 
 * The module calls LineSweep for all Seeds in the SeedsVector.
 * @ingroup module
 */
class ExecOnVec: public CppModule
{
private:
    std::shared_ptr<CppModule> pModule;

public:
    ExecOnVec(std::shared_ptr<CppModule> pModule)
            :
        CppModule(),
        pModule(pModule)
    {}//default constructor

    //overload
    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput)
    {
        std::shared_ptr<ContainerVector> pSeedsVector = std::shared_ptr<SeedsVector>(
            std::static_pointer_cast<SeedsVector>(vpInput[0]));

        std::shared_ptr<SeedsVector> vTempResults = std::shared_ptr<SeedsVector>(
                new SeedsVector(pSeedsVector->size())
            );
        {
            ThreadPool xPool( NUM_THREADS_ALIGNER );
            for(unsigned int i = 0; i < vTempResults->size(); i++)
            {
                (*vTempResults)[i] = std::shared_ptr<Seeds>();
                xPool.enqueue(
                    [&vTempResults]
                    (
                        size_t, 
                        std::shared_ptr<SeedsVector> pSeedsVector,
                        std::shared_ptr<CppModule> pModule,
                        unsigned int i
                    )
                    {
                        std::vector<std::shared_ptr<Container>> vInput = {
                                (*pSeedsVector)[i]
                            };
                        try
                        {
                            (*vTempResults)[i] = std::static_pointer_cast<Seeds>(
                                    pModule->execute(vInput));
                            (*vTempResults)[i]->mem_score = (*vTempResults)[i]->getScore();
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
                            std::cerr << "linesweep deleviered nullpointer as result" << std::endl;
                        if((*vTempResults)[i] == nullptr)
                            throw NullPointerException("linesweep deleviered nullpointer as result");
                    },//lambda
                    pSeedsVector, pModule, i
                );
            }//for
        }//scope xPool

        if(vTempResults->empty())
            return std::shared_ptr<Seeds>(new Seeds());

        std::sort(
            vTempResults->begin(), vTempResults->end(),
            []
            (const std::shared_ptr<Seeds> a, const std::shared_ptr<Seeds> b)
            {
                return a->mem_score > b->mem_score;
            }//lambda
        );//sort function call

        return vTempResults;
    }//function
    
	std::vector<std::shared_ptr<Container>> getInputType();
    std::shared_ptr<Container> getOutputType();

};//class

/**
 * @brief Exposes the SweepAllReturnBest @ref CppModule "module" to boost python.
 * @ingroup export
 */
void exportExecOnVector();

#endif
