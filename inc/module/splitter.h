#ifndef SPLITTER_H
#define SPLITTER_H

#include "module/module.h"
#include "container/nucSeq.h"
#include <fstream>

namespace libMA
{

    class Splitter: public Module
    {
    public:
        std::shared_ptr<ContainerVector> pVec;

        Splitter(std::shared_ptr<ContainerVector> pVec)
                :
            pVec(pVec)
        {}//constructor

        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Nil
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(NucSeq)
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "Splitter";
        }//function

        bool outputsVolatile() const
        {
            return true;
        }//function

    };//class

    class Lock: public Module
    {
    public:
        std::shared_ptr<Container> pType;

        Lock(std::shared_ptr<Container> pType)
                :
            pType(pType)
        {}//constructor

        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Nil
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(NucSeq)
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "Lock";
        }//function

    };//class

    class UnLock: public Module
    {
    public:
        std::shared_ptr<Pledge> pLockPledge;

        UnLock(std::shared_ptr<Pledge> pLockPledge)
                :
            pLockPledge(pLockPledge)
        {}//constructor

        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Nil
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(NucSeq)
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "UnLock";
        }//function

    };//class

}//namespace

void exportSplitter();

#endif