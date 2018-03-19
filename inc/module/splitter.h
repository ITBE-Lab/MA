#ifndef SPLITTER_H
#define SPLITTER_H

#include "module/module.h"
#include "container/nucSeq.h"
#include <fstream>
#include <mutex>

namespace libMA
{

    class Splitter: public Module
    {
    public:
        std::shared_ptr<Pledge> pVec;

        Splitter(std::shared_ptr<Pledge> pVec)
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

        std::string getFullDesc() const
        {
            return std::string("Splitter");
        }//function

    };//class

    class Collector: public Module
    {
    public:
        std::shared_ptr<ContainerVector> pVec;
        std::shared_ptr<std::mutex> pLock;

        Collector(std::shared_ptr<Container> pType)
                :
            pVec(new ContainerVector(pType)),
            pLock(new std::mutex)
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
            return "Collector";
        }//function

        std::string getFullDesc() const
        {
            return std::string("Collector");
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

        std::string getFullDesc() const
        {
            return std::string("Lock");
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

        std::string getFullDesc() const
        {
            return "UnLock";
        }//function

        bool outputsVolatile() const
        {
            return true;
        }//function
    };//class

}//namespace

void exportSplitter();

#endif