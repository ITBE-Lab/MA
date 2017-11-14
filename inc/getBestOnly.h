//TODO: docu


/*
 * FIXME: i'm super inefficient
 * would be better to make all following modules be able to deal with seedVectors instrad of seeds
 */

#ifndef GET_BEST_ONLY_H
#define GET_BEST_ONLY_H

#include "linesweep.h"

/**
 * @brief Execute LineSweep for all given seeds.
 * @details 
 * The module calls LineSweep for all Seeds in the SeedsVector.
 * @ingroup module
 */
class GetBestOnly: public CppModule
{
public:

    //overload
    std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput)
    {
        std::shared_ptr<SeedsVector> pSeedsVector = std::shared_ptr<SeedsVector>(
            std::static_pointer_cast<SeedsVector>(vpInput[0]));

        std::shared_ptr<Seeds> pRet = std::shared_ptr<Seeds>(new Seeds());
        if(!pSeedsVector->empty())
            pRet->append((*pSeedsVector)[0]);

        return pRet;
    }//function

    //overload
    std::vector<ContainerType> getInputType()
    {
        return std::vector<ContainerType>{
                ContainerType::seedsVector
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
void exportGetBestOnly();

#endif
