//TODO: docu

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

        return std::shared_ptr<Seeds>(new Seeds((*pSeedsVector)[0]));
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
        return ContainerType::seedsVector;
    }//function
};//class

/**
 * @brief Exposes the SweepAllReturnBest @ref CppModule "module" to boost python.
 * @ingroup export
 */
void exportGetBestOnly();

#endif
