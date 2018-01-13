#ifndef PIPE_H
#define PIPE_H

#include "module/module.h"

namespace libMABS
{

    /**
     * @brief Used to apply a computational graph to multile inputs
     * @ingroup module
     * @details
     * Requires n+1 pledges from the computation graph:
     * sets the first n pledges using its container vector input
     * Then gets the n+1th pledge.
     * once the result is obtained the first n pledges are reset 
     * with the next elements from the vector.
     * This way one output for each element in the vector is created.
     * @todo at the moment the lock for pipes is public this could be much finer...
     */
    class Pipe: public Module
    {
    public:
        unsigned int uiN;
        std::vector<std::shared_ptr<Pledge>> aIn;
        std::shared_ptr<Pledge> pOut;
        static std::mutex xLock;

        Pipe( 
                std::vector<std::shared_ptr<Pledge>> aIn,
                std::shared_ptr<Pledge> pOut
            )
                :
            aIn(aIn),
            pOut(pOut)
        {
            uiN = aIn.size();
            assert(uiN != 0);
        }//constructor

        std::shared_ptr<Container> execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - ContainerVector(x)
         * - ...
         * - ContainerVector(x)
         * n times
         * n is a parameter that can be chosen when the pledge is created
         */
        ContainerVector getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(x)
         */
        std::shared_ptr<Container> getOutputType() const;

        std::string getName() const
        {
            return "Pipe";
        }
    };//class
}//namspace libMABS

/**
 * @brief export the Pipe @ref Module "module" to python.
 * @ingroup export
 */
void exportPipe();


#endif
