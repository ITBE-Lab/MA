/** 
 * @file cppModule.h
 * @brief Implements a abstract Module class.
 * @author Markus Schmidt
 */
#ifndef MODULE_H
#define MODULE_H

#include "container.h"
#include "debug.h"
#include <memory>
#include <Python.h>
#include <iostream>
#include <boost/python/list.hpp>
#include "threadPool.h"


extern std::mutex xPython;
/**
 * @defgroup module
 * @brief All classes implementing some algorithm.
 * @details
 * These classes should all inherit from Module.
 */

class Pledge;

/**
 * @brief Checks weather the given types are equal.
 * @details
 * This requires a function since we have types like any of none.
 */
bool typeCheck(
        ContainerType xData, 
        ContainerType xExpected
    );

/**
 * @brief Checks weather the given container has the expected type.
 * @details
 * This requires a function since we have types like any of none.
 */
bool typeCheck(
        std::shared_ptr<Container> pData, 
        ContainerType xExpected
    );

/**
 * @brief Checks weather given containers have the expected types.
 * @details
 * This requires a function since we have types like any of none.
 */
bool typeCheck(
        std::vector<std::shared_ptr<Container>> vData, 
        std::vector<ContainerType> vExpected
    );

/**
 * @brief Abstract class intended for the implementaiton of various algorithms.
 * @details
 * All computing on data should inherit from this class
 * @see the Python implementation of @ref LAuS.aligner.Module "module".
 * @ingroup module
 */
class CppModule
{
public:

    ~CppModule()
    {
        std::cerr << "Module destroyed" << std::endl;
    }

    /**
     * @brief Execute the implemented algorithm.
     * @details
     * Expects the given containers to have the correct types.
     */
    virtual std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> pInput)
    {
        return nullptr;
    }

    /**
     * @brief The expected input types.
     * @details
     * Used for type checking the inputs before calling execute.
     */
    virtual std::vector<ContainerType> getInputType()
    {
        return std::vector<ContainerType>{ContainerType::nothing};
    }

    /**
     * @brief The expected output type.
     * @details
     * Used for type checking weather the module returns expected data.
     */
    virtual ContainerType getOutputType()
    {
        return ContainerType::nothing;
    }

    /**
     * @brief Execute the implemented algorithm.
     * @details
     * Internally calls execute after checking the input types.
     * Also checks the result returned by Execute.
     */
    std::shared_ptr<Container> saveExecute(std::vector<std::shared_ptr<Container>> vInput)
    {
        if(!typeCheck(vInput, getInputType()))
            throw new ModuleIO_Exception("Input type and expected input type did not match.");
        std::shared_ptr<Container> pRet = execute(vInput);
        if(!typeCheck(pRet, getOutputType()))
            throw new ModuleIO_Exception("Module produced output of wrong type.");
        return pRet;
    }//function

    /**
     * @brief Execute the implemented algorithm.
     * @details
     * Internally calls execute after checking the input types.
     * Also checks the result returned by Execute.
     */
    std::shared_ptr<Container> pyExecute(std::vector<std::shared_ptr<Container>> vInput)
    {
        std::cerr << "CALLED BY PYTHON" << std::endl;
        xPython.unlock();
        std::shared_ptr<Container> pRet = saveExecute(vInput);
        std::cerr << "TRYING TO RETURN TO PYTHON" << std::endl;
        xPython.lock();
        std::cerr << "RETURNING TO PYTHON" << std::endl;
        return pRet;
    }//function

    /**
     * @brief Make this module promise to execute it's function on the provided data.
     * @details
     * see Pledge for more information about the computational graph that can be build using 
     * promiseMe and Pledge.
     */
    static std::shared_ptr<Pledge> promiseMe(
            std::shared_ptr<CppModule> pThis,
            std::vector<std::shared_ptr<Pledge>> vInput
        );
};

/**
 * @brief Abstract class intended to hold promises to data objects used by Modules.
 * @details
 * Us this class to set up a computational graph.
 * The graphs nodes are Modules (in python or cpp) that perform computational tasks.
 * Each node takes a vector of containers as input.
 * This container vector was created using the vector of Pledges that was given to the Module.
 * Then the Module returns a Container that fits the pledge returned by the Modules 
 * promiseMe function.
 *
 * @note The advantage of the computational graph is that there is no unnecessary 
 * jumping between python and cpp code.
 * 
 * @section comp_graph_sec Computational Graph Quick Start
 * 
 * Here is some python code that sets up the three main @ref Module "Modules" 
 * reqired for alignment, using a computation graph: <br>
 * Note that this setup mimics the one given in the @ref quick_start_sec "quick start" section
 * @code{.py}
 * # A module that creates seeds.
 * seg = Segmentation()
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
 * # Add the segmentation module to the computation graph.
 * segments_pledge = seg.promise_me((fm_index_pledge, query_pledge)) #(*)
 * # Call the line sweep module.
 * seeds_pledge = sweep.promise_me((segments_pledge,)) #(*)
 * # Call the local alignment module.
 * alignment_pledge = nmw.promise_me((seeds_pledge, query_pledge, ref_pledge)) #(*)
 * 
 * # Print the alignment
 * print_pledge = printer.promise_me((alignment_pledge, )) #(*)
 * @endcode
 *
 * Here we make add out modules to our computation graph. <br>
 * (*) All @ref Module "modules" use a tuple or list of @ref Pledge "pledges" as input.
 * Therefore the @ref Module::promiseMe "promise_me" function calls use the syntax shown above.
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
 * # Trigger the aligment process.
 * print_pledge.next()
 * @endcode
 * 
 * Here we fullfill the pledges we made and then trigger the alignment process.
 * Note that the AlignmentPrinter and SweepAllReturnBest modules are implemented in Python,
 * while all other modules are implemented in C++.
 * The computational graph is able to jump between python and C++ modules and containers as needed.
 *
 * @ingroup container
 */
class Pledge : public Container
{
private:
    std::shared_ptr<CppModule> pledger;
    boost::python::object py_pledger;
    std::shared_ptr<Container> content;
    ContainerType type;
    std::vector<std::shared_ptr<Pledge>> vPredecessors;
    std::vector<std::weak_ptr<Pledge>> vSuccessors;
    std::mutex xMutex;

    /**
     * @brief Create a new pledge with a cpp module responsible to fullfill it.
     * @details
     * This means that this Pledge can be automatically fullfilled by the given module.
     */
    Pledge(
            std::shared_ptr<CppModule> pledger,
            ContainerType type,
            std::vector<std::shared_ptr<Pledge>> vPredecessors
        )
            :
        pledger(pledger),
        py_pledger(),
        content(),
        type(type),
        vPredecessors(vPredecessors),
        vSuccessors(),
        xMutex()
    {
        std::cerr << "Pledge created " << type << std::endl;
    }//constructor

    /**
     * @brief Create a new pledge with a Python module responsible to fullfill it.
     * @details
     * This means that this Pledge can be automatically fullfilled by the given module.
     */
    Pledge(
            boost::python::object py_pledger,
            ContainerType type,
            std::vector<std::shared_ptr<Pledge>> vPredecessors
        )
            :
        pledger(),
        py_pledger(py_pledger),
        content(),
        type(type),
        vPredecessors(vPredecessors),
        vSuccessors(),
        xMutex()
    {
        std::cerr << "Pledge created " << type << std::endl;
    }//constructor
public:
    /**
     * @brief Create a new pledge without a module giving the pledge.
     * @details
     * This means that this Pledge has to be fullfilled by calling set manually.
     */
    Pledge(
            ContainerType type
        )
            :
        pledger(nullptr),
        py_pledger(),
        content(),
        type(type),
        vPredecessors(),
        vSuccessors(),
        xMutex()
    {
        std::cerr << "Pledge created " << type << std::endl;
    }//constructor
    
    ~Pledge()
    {
        std::cerr << "Pledge destroyed " << type << std::endl;
    }//destructor

    /**
     * @brief this is required due to the use of mutex 
    */
    Pledge(const Pledge&) = delete;//copy constructor

    static inline std::shared_ptr<Pledge> makePledge(
            std::shared_ptr<CppModule> pledger,
            std::vector<std::shared_ptr<Pledge>> vPredecessors
        )
    {
        std::shared_ptr<Pledge> pRet = std::shared_ptr<Pledge>(
            new Pledge(pledger, pledger->getOutputType(), vPredecessors)
        );
        
        for(std::shared_ptr<Pledge> pPredecessor : vPredecessors)
            pPredecessor->vSuccessors.push_back(std::weak_ptr<Pledge>(pRet));

        return pRet;
    }//contructor function

    static inline std::shared_ptr<Pledge> makePyPledge(
            boost::python::object py_pledger,
            ContainerType type,
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
        std::lock_guard<std::mutex> xGuard(xMutex);
        if(content != nullptr)
            return content;
        if(pledger == nullptr && py_pledger.is_none())
            throw ModuleIO_Exception("No pledger known");
        if(pledger != nullptr)
        {
            std::vector<std::shared_ptr<Container>> vInput;
            for(std::shared_ptr<Pledge> pFuture : vPredecessors)
                vInput.push_back(pFuture->get());
            content = (std::shared_ptr<Container>)pledger->execute(vInput);
            assert(typeCheck(content, type));
        }//if
        else
        {
            boost::python::list vInput;
            for(std::shared_ptr<Pledge> pFuture : vPredecessors)
            {
                std::shared_ptr<Container> pInput = pFuture->get();
                assert(pInput != nullptr);
                vInput.append(pInput);
            }
            /*
            * here we jump to python code to call a function and resume the cpp code 
            * once python is done...
            */
            try
            {
                std::cerr << "WAITING ON PYTHON LOCK" << std::endl;
                xPython.lock();
                std::cerr << "CALLING PYTHON" << std::endl;
                Container* extract = boost::python::extract<
                        Container*
                    >(
                        py_pledger.attr("execute")(vInput)
                    ); 
                std::cerr << "PYTHON RETURNED" << std::endl;
                assert(extract != nullptr);
                content = std::shared_ptr<Container>(extract);
            } catch (const std::system_error& e) {
                std::cerr << e.what() << std::endl;
            } catch (const std::exception& e) {
                std::cerr << e.what() << std::endl;
                xPython.unlock();
            } catch (const std::string& e) {
                std::cerr << e << std::endl;
            } catch (...) {
                std::cerr << "unknown exception" << std::endl;
                boost::python::handle_exception();
                exit(0);
            }//catch
            xPython.unlock();
            std::cerr << "PYTHON RETURNED" << std::endl;
            assert(typeCheck(content, type));
        }//else
        return content;
    }//function

    static inline std::vector<std::shared_ptr<Container>> simultaneousGet(
            std::vector<std::shared_ptr<Pledge>> vPledges,
            unsigned int numThreads
        )
    {
        std::vector<std::shared_ptr<Container>> vRet = std::vector<std::shared_ptr<Container>>(
                vPledges.size()
            );
        {
            ThreadPool xPool(numThreads);
            for(unsigned int i = 0; i < vPledges.size(); i++)
            {
                xPool.enqueue(
                    []
                    (size_t, std::shared_ptr<Pledge> pPledge, std::shared_ptr<Container> pResult)
                    {
                        assert(pPledge != nullptr);
                        pResult = pPledge->get();
                    },//lambda
                    vPledges[i], vRet[i]
                );
            }//for
        }//scope xPool

        return vRet;
    }//function

    //overload
    ContainerType getType()
    {
        return type;
    }//function
};//class



/**
 * @brief Exposes the Module class to boost python.
 * @ingroup export
 */
void exportModule();

#endif