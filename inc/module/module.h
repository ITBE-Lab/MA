/** 
 * @file module.h
 * @brief Implements a abstract Module class.
 * @author Markus Schmidt
 */
#ifndef MODULE_H
#define MODULE_H

#include "container/container.h"
#include "util/debug.h"
#include <memory>
#include <Python.h>
#include <iostream>
#include <boost/python/list.hpp>
#include "util/threadPool.h"
#include <ctime>
#include <chrono>

#define PYTHON_MODULES_IN_COMP_GRAPH ( false )

/**
 * @brief the C++ code is in this namespace.
 */
namespace libMA
{
    extern std::mutex xPython;
    /**
     * @defgroup module
     * @brief All classes implementing some algorithm.
     * @details
     * These classes should all inherit from Module.
     */

    class Pledge;


    /**
     * @brief Checks weather the given container has the expected type.
     * @details
     * This requires a function since we have types like any of none.
     */
    bool EXPORTED typeCheck(
            std::shared_ptr<Container> pData, 
            std::shared_ptr<Container> pExpected
        );

    /**
     * @brief Checks weather given containers have the expected types.
     * @details
     * This requires a function since we have types like any of none.
     */
    bool EXPORTED typeCheck(
            ContainerVector vData, 
            ContainerVector vExpected
        );

    /**
     * @brief Abstract class intended for the implementaiton of various algorithms.
     * @details
     * All computing on data should inherit from this class
     * @see the Python implementation of @ref MA.aligner.Module "module".
     * @ingroup module
     */
    class Module
    {
    public:

        /**
         * @brief Execute the implemented algorithm.
         * @details
         * Expects the given containers to have the correct types.
         */
        virtual std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> pInput)
        {
            return nullptr;
        }//function

        /**
         * @brief The expected input types.
         * @details
         * Used for type checking the inputs before calling execute.
         */
        virtual ContainerVector EXPORTED getInputType() const
        {
            return ContainerVector{std::shared_ptr<Container>(new Nil())};
        }//function

        /**
         * @brief The expected output type.
         * @details
         * Used for type checking weather the module returns expected data.
         */
        virtual std::shared_ptr<Container> EXPORTED getOutputType() const
        {
            return std::shared_ptr<Container>(new Nil());
        }//function

        /**
         * @brief Return the name of the Module.
         * @details
         * Can be used to print a representation of the computational graph.
         * Is also used when creating error messages.
         */
        virtual std::string EXPORTED getName() const
        {
            return "Module";
        }//function

        /**
         * @brief Can the module return different output for the same input.
         * @details
         * For example a file reader may create different output if called twice on the same file.
         * If this is set to true the pledge will recompute the result each time get is called.
         */
        virtual bool EXPORTED outputsVolatile() const
        {
            return false;
        }//function

        /**
         * @brief Execute the implemented algorithm.
         * @details
         * Internally calls execute after checking the input types.
         * Also checks the result returned by Execute.
         */
        std::shared_ptr<Container> saveExecute(std::shared_ptr<ContainerVector> vInput)
        {
            if(!typeCheck(*vInput, getInputType()))
            {
                std::cerr << "input of module " << getName() << " had the wrong type" << std::endl;
                throw new ModuleIO_Exception("Input type and expected input type did not match.");
            }//if
            std::shared_ptr<Container> pRet = execute(vInput);
            if(!typeCheck(pRet, getOutputType()))
            {
                std::cerr << "output of module " << getName() << " had the wrong type" << std::endl;
                throw new ModuleIO_Exception("Module produced output of wrong type.");
            }//if
            return pRet;
        }//function

        /**
         * @brief Execute the implemented algorithm.
         * @details
         * Internally calls execute after checking the input types.
         * Also checks the result returned by Execute.
         */
        std::shared_ptr<Container> pyExecute(std::shared_ptr<ContainerVector> vInput)
        {
#if 1
            try
            {
                std::shared_ptr<Container> pRet = saveExecute(vInput);
                return pRet;
            } catch(ModuleIO_Exception e) {
                std::cerr << e.what() << std::endl;
                std::cerr << "cant return to python - cpp module failed" << std::endl;
                exit(0);
            } catch (const std::exception& e) {
                std::cerr << e.what() << std::endl;
            } catch (const std::string& e) {
                std::cerr << e << std::endl;
            } catch (...) {
                std::cerr << "unknown exception" << std::endl;
                boost::python::handle_exception();
            }//catch
            return nullptr;
#else
            throw AlignerException("python modules are not allowed to call cpp modules currently - sorry");
#endif
        }//function

        /**
         * @brief Make this module promise to execute it's function on the provided data.
         * @details
         * see Pledge for more information about the computational graph that can be build using 
         * promiseMe and Pledge.
         * @note This is static since we need to save a reference to the shared_ptr of the Module 
         * promising.
         */
        static std::shared_ptr<Pledge> EXPORTED promiseMe(
                std::shared_ptr<Module> pThis,
                std::vector<std::shared_ptr<Pledge>> vInput
            );
    };

    /**
     * @page comp_graph_page Computational Graph Quick Start
     *
     * Here is some python code that sets up the three main @ref Module "Modules" 
     * reqired for alignment, using a computation graph: <br>
     * Note that this setup mimics the one given in the @ref quick_start_sec "quick start" section
     * @code{.py}
     * # A module that creates seeds.
     * seg = BinarySeeding()
     * # A module that removes inconsistent seeds.
     * sweep = SweepAllReturnBest()
     * # A module that creates local alignments in the gaps between the seeds.
     * nmw = NeedlemanWunsch()
     * # A module that prints the alignment to the console.
     * printer = AlignmentPrinter()
     * @endcode
     *
     * Here we set up the @ref Module "modules" we need for the alignment process.
     * The @ref Module "modules" themselves do not store data. They can be used multiple times.
     *
     * @code{.py}
     * # Setup a pledge for the suffix array
     * fm_index_pledge = Pledge(ContainerType.fM_index)
     *
     * # Setup a pledge for the packed reference
     * ref_pledge = Pledge(ContainerType.packedNucSeq)
     *
     * # Setup a pledge for the query sequence.
     * query_pledge = Pledge(ContainerType.nucSeq)
     * @endcode
     *
     * Here we setup the @ref Pledge "pledges", promising that we will deliver some data.
     *
     * @code{.py}
     * # Add a segmentation module to the computation graph.
     * segments_pledge = seg.promise_me(fm_index_pledge, query_pledge)
     * # Call the line sweep module.
     * seeds_pledge = sweep.promise_me(segments_pledge)
     * # Call the local alignment module.
     * alignment_pledge = nmw.promise_me(seeds_pledge, query_pledge, ref_pledge)
     *
     * # Print the alignment
     * print_pledge = printer.promise_me(alignment_pledge)
     * @endcode
     *
     * Here we make add out modules to our computation graph. <br>
     *
     * @code{.py}
     * # Setup a container for the suffix array
     * fm_index = FMIndex()
     * # Load the array from a file.
     * fm_index.load("filename")
     *
     * # Setup a container for the packed reference
     * ref = Pack()
     * # Load the pack from a file.
     * ref.load("filename")
     *
     * # Create a query string.
     * query_string = "ACCTAA"
     * # Setup a container for the query sequence.
     * query = NucSeq(query_string)
     * @endcode
     *
     * Here we load our actual data.
     *
     * @code{.py}
     * #fullfill the FM-Index pledge we gave.
     * fm_index_pledge.set(fm_index)
     * #fullfill the reference pledge we gave.
     * ref_pledge.set(ref)
     * #fullfill the query pledge we gave.
     * query_pledge.set(query)
     *
     * # Trigger the alignment process.
     * print_pledge.next()
     * @endcode
     *
     * Here we fullfill the pledges we made and then trigger the alignment process.
     * Note that the AlignmentPrinter and SweepAllReturnBest modules are implemented in Python,
     * while all other modules are implemented in C++.
     * The computational graph is able to jump between python and
     * C++ modules and containers as needed.
     * Note that this setup does not perform well, since we skipped the seed filtering step.
     *
     */
    /**
     * @brief Abstract class intended to hold promises to data objects used by Modules.
     * @details
     * Use this class to set up a computational graph.
     * The graphs nodes are Modules (in python or cpp) that perform computational tasks.
     * Each node takes a vector of containers as input.
     * This container vector was created using the vector of Pledges that was given to the Module.
     * Then the Module returns a Container that fits the pledge returned by the Modules 
     * promiseMe function.
     *
     * Check out the @ref comp_graph_page "Computational Graph Quick Start"
     *
     * @note The advantage of the computational graph is that there is no unnecessary 
     * jumping between python and cpp code.
     *
     * @ingroup container
     */
    class Pledge : public Container
    {
    private:
        std::shared_ptr<Module> pledger;
        boost::python::object py_pledger;
        std::shared_ptr<Container> content;
        std::shared_ptr<Container> type;
        std::vector<std::shared_ptr<Pledge>> vPredecessors;
        std::vector<std::weak_ptr<Pledge>> vSuccessors;
        std::vector<std::shared_ptr<Pledge>> aSync;

    public:
        std::shared_ptr<std::mutex> pMutex;
        double execTime;
    private:

        /**
         * @brief Create a new pledge with a cpp module responsible to fullfill it.
         * @details
         * This means that this Pledge can be automatically fullfilled by the given module.
         */
        Pledge(
                std::shared_ptr<Module> pledger,
                std::shared_ptr<Container> type,
                std::vector<std::shared_ptr<Pledge>> vPredecessors
            )
                :
            pledger(pledger),
            py_pledger(),
            content(),
            type(type),
            vPredecessors(vPredecessors),
            vSuccessors(),
            aSync(),
            pMutex(new std::mutex),
            execTime(0)
        {}//constructor

        /**
         * @brief Create a new pledge with a Python module responsible to fullfill it.
         * @details
         * This means that this Pledge can be automatically fullfilled by the given module.
         */
        Pledge(
                boost::python::object py_pledger,
                std::shared_ptr<Container> type,
                std::vector<std::shared_ptr<Pledge>> vPredecessors
            )
                :
            pledger(),
            py_pledger(py_pledger),
            content(),
            type(type),
            vPredecessors(vPredecessors),
            vSuccessors(),
            aSync(),
            pMutex(new std::mutex),
            execTime(0)
        {}//constructor

        void execForGet()
        {
            bool bVolatile = false;
            //@todo same for py_pledger here
            if(pledger != nullptr)
                bVolatile = pledger->outputsVolatile();

            //in this case there is no need to execute again
            if(!bVolatile && content != nullptr)
                return;

            if(pledger == nullptr && py_pledger.is_none())
                throw ModuleIO_Exception("No pledger known");
            if(pledger != nullptr)
            {
                bool bDry = false;
                std::shared_ptr<ContainerVector> vInput(new ContainerVector());
                for(std::shared_ptr<Pledge> pFuture : vPredecessors)
                {
                    //here we execute all previous modules in the comp graph
                    vInput->push_back(pFuture->get());
                    if(vInput->back()->bDry)
                        bDry = true;
                    //keep the mutex locked while the execute runs
                    pFuture->pMutex->lock();
                }//for
                try
                {
                    auto timeStamp = std::chrono::system_clock::now();
                    //actually execute the module
                    content = pledger->execute(vInput);
                    std::chrono::duration<double> duration = std::chrono::system_clock::now() - timeStamp;
                    execTime = duration.count();
                    assert(typeCheck(content, type));
                    //content->bDry may be true already we do not want to overwrite it with false
                    if(bDry)
                        content->bDry = true;
                } catch(...)
                {
                    std::cerr << "unknown exception during execution" << std::endl;
                }//catch
                //unlock all locked mutex
                for(std::shared_ptr<Pledge> pFuture : vPredecessors)
                    pFuture->pMutex->unlock();
            }//if
            else
            {
                boost::python::list vInput;
                for(std::shared_ptr<Pledge> pFuture : vPredecessors)
                    vInput.append(pFuture->get());
                /*
                * here we jump to python code to call a function and resume the cpp code 
                * once python is done...
                */
                auto timeStamp = std::chrono::system_clock::now();
                content = boost::python::extract<
                        std::shared_ptr<Container>
                    >(
                        py_pledger.attr("save_execute")(vInput)
                    ); 
                std::chrono::duration<double> duration = std::chrono::system_clock::now() - timeStamp;
                execTime = duration.count();
                assert(typeCheck(content, type));
            }//else
        }//function

    public:
        /**
         * @brief Create a new pledge without a module giving the pledge.
         * @details
         * This means that this Pledge has to be fullfilled by calling set manually.
         */
        Pledge(
                std::shared_ptr<Container> type
            )
                :
            pledger(nullptr),
            py_pledger(),
            content(),
            type(type),
            vPredecessors(),
            vSuccessors(),
            pMutex(new std::mutex),
            execTime(0)
        {}//constructor

        /**
         * @brief this is required due to the use of mutex 
         */
        Pledge(const Pledge&) = delete;//copy constructor
        
        //overload
        bool canCast(std::shared_ptr<Container> c) const
        {
            return std::dynamic_pointer_cast<Pledge>(c) != nullptr;
        }//function

        //overload
        std::string getTypeName() const
        {
            return "Pledge";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            return type->getType();
        }//function

        const std::shared_ptr<Module> getPledger() const
        {
            return pledger;
        }//function

        static EXPORTED std::shared_ptr<Pledge> makePledge(
                std::shared_ptr<Module> pledger,
                std::vector<std::shared_ptr<Pledge>> vPredecessors
            );

        static inline std::shared_ptr<Pledge> makePyPledge(
                boost::python::object py_pledger,
                std::shared_ptr<Container> type,
                std::vector<std::shared_ptr<Pledge>> vPredecessors
            )
        {
            std::shared_ptr<Pledge> pRet = std::shared_ptr<Pledge>(
                new Pledge(py_pledger, type, vPredecessors)
            );
            
            for(std::shared_ptr<Pledge> pPredecessor : vPredecessors)
                pPredecessor->vSuccessors.push_back(std::weak_ptr<Pledge>(pRet));

            return pRet;
        }//contructor function

        /**
         * @brief Manually fullfill this pledge.
         * @details
         * Invalidates the content of all following pledges.
         * @note It would be better to have a "constant" or "variable" class instead of this function.
         */
        void set(std::shared_ptr<Container> c)
        {
            content = c;
            for(std::weak_ptr<Pledge> pSuccessor : vSuccessors)
            {
                std::shared_ptr<Pledge> lock = pSuccessor.lock();
                if(lock != nullptr)
                    lock->set(nullptr);
            }//for
        }//function

        /**
         * @brief Warpper for boost python
         * @details
         * transfers the ownership of the given container from python to the c++ code
         */
        void set_boost(std::shared_ptr<Container> c)
        {
            std::shared_ptr<Container> swap;
            swap.swap(c);
            set(swap);
        }//function

        /**
         * @brief Get the promised container.
         * @details
         * Checks weather the pledge has already been fullfilled for this call.
         * If necessary, uses the python or cpp module to fullfill the pledge,
         * and returns the respective container.
         */
        std::shared_ptr<Container> get()
        {
            //multithreading is possible thus a guard is required here.
            //deadlock prevention is trivial, since the computational graphs are essentially trees.
            std::lock_guard<std::mutex> xGuard(*pMutex);

            //execute
            execForGet();
            //also execute all the synchronized locks
            for(std::shared_ptr<Pledge> pSync : aSync)
                pSync->execForGet();

            return content;
        }//function

        /**
         * @brief Synchronize two pledges.
         * @details
         * They will always be executed right one after another.
         * This can be used to for paired reads,
         * where we want to read one read from file a and one from file b.
         * We need to make sure that we always read the pairs together.
         */
        static inline void synchronize(std::shared_ptr<Pledge> pA, std::shared_ptr<Pledge> pB)
        {
            //add each other to the caller list
            pA->aSync.push_back(pB);
            pB->aSync.push_back(pA);
            //sync the mutex of pA and pB
            pA->pMutex = pB->pMutex;
        }//function

        void forAllSyncs(std::function<void(std::shared_ptr<Pledge>)> fDo)
        {
            for(std::shared_ptr<Pledge> pSync : aSync)
                fDo(pSync);
        }//function

        /**
         * @brief Gets the given pledges simultaneously.
         * @details
         * If bLoop is true the threads will keep going until all volatile modules are dry
         * if numThreads is not specified numThreads will be equal to the amount of pledges given
         */
        static inline void simultaneousGet(
                std::vector<std::shared_ptr<Pledge>> vPledges,
                bool bLoop = false,
                unsigned int numThreads = 0
            )
        {
            if(numThreads == 0)
                numThreads = vPledges.size();
            DEBUG(
                std::cout <<"will cause crashes if used on python modules (@todo fix that)"<< std::endl;
            )
            {
                ThreadPool xPool(numThreads);
                for(std::shared_ptr<Pledge> pPledge : vPledges)
                {
                    xPool.enqueue(
                        [&bLoop]
                        (size_t, std::shared_ptr<Pledge> pPledge)
                        {
                            assert(pPledge != nullptr);

                            //if bLoop is true this loop will keep going until all 
                            //volatile modules are dry.
                            //otherwise one iteration is performed only
                            while(!pPledge->get()->bDry && bLoop);
                        },//lambda
                        pPledge
                    );
                }//for
            }//scope xPool
        }//function
    };//class
}//namespace libMA


/**
 * @brief Exposes the Module class to boost python.
 * @ingroup export
 */
void exportModule();

#endif