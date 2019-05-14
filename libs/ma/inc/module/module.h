/**
 * @file module.h
 * @brief Implements a abstract Module class.
 * @author Markus Schmidt
 */
#ifndef MODULE_H
#define MODULE_H

#include "container/container.h"
#include "util/default_parameters.h"
#include "util/parameter.h"
#include "util/threadPool.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#ifdef WITH_PYTHON
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif

#include <algorithm>
#include <bitset>
#include <chrono>
#include <ctime>
#include <iostream>
#include <thread>
#include <typeinfo>
/// @endcond

#define PYTHON_MODULES_IN_COMP_GRAPH ( false )

/**
 * @brief the C++ code is in this namespace.
 */
namespace libMA
{

/* +++++++++++++++++++++++++++++++++++++++++++++++++
 * 4. Unpack a tuple to call a matching function.
 * FIX ME: Better naming scheme, more documentation.
 * http://stackoverflow.com/questions/7858817/unpacking-a-tuple-to-call-a-matching-function-pointer
 */
template <int...> struct seq
{};

template <int N, int... S> struct gens : gens<N - 1, N - 1, S...>
{};

template <int... S> struct gens<0, S...>
{
    typedef seq<S...> type;
};


/**
 * @defgroup module
 * @brief All classes implementing some algorithm.
 * @details
 * These classes should all inherit from Module.
 */

/**
 * @brief Abstract class intended for the implementaiton of various algorithms.
 * @details
 * All computing on data should inherit from this class
 * @see the Python implementation of @ref MA.aligner.Module "module".
 * @ingroup module
 */
template <typename TP_RETURN_, bool IS_VOLATILE_, typename... TP_ARGUMENTS> class Module
{
  private:
    std::bitset<1> xFlags;
    static constexpr size_t uiFinishedFlag = 0;

  public:
    typedef TP_RETURN_ TP_RETURN;
    static constexpr bool IS_VOLATILE = IS_VOLATILE_;
    static constexpr size_t NUM_ARGUMENTS = sizeof...( TP_ARGUMENTS );
    typedef std::tuple<std::shared_ptr<TP_ARGUMENTS>...> TP_TUPLE_ARGS;

  private:
    template <size_t IDX> struct ArgNamesCollector
    {
        bool operator( )( std::string& sArgNames, TP_TUPLE_ARGS& rArgs )
        {
            sArgNames += type_name( std::get<IDX>( rArgs ) ) + ", ";
            return true;
        } // operator
    }; // struct

  public:
    /**
     * @brief Execute the implemented algorithm.
     * @details
     * Expects the given containers to have the correct types.
     */
    virtual std::shared_ptr<TP_RETURN> EXPORTED execute( std::shared_ptr<TP_ARGUMENTS>... args )
    {
        std::string sArgNames;
        auto rTup = std::make_tuple( args... );
        TemplateLoop<NUM_ARGUMENTS, ArgNamesCollector>::iterate( sArgNames, rTup );
        throw AnnotatedException( type_name( this ) + " did not implement execute with return type: " +
                                  type_name<TP_RETURN>( ) + " and parameters: " + sArgNames );
        return nullptr;
    } // method

    inline void setFinished( )
    {
        if( IS_VOLATILE )
            xFlags.set( uiFinishedFlag );
    } // method

  private:
    template <int... S> std::shared_ptr<TP_RETURN> callFunc( const TP_TUPLE_ARGS& tParams, seq<S...> )
    {
        return this->execute( std::get<S>( tParams )... );
    } // method

  public:
    virtual std::shared_ptr<TP_RETURN> EXPORTED executeTup( const TP_TUPLE_ARGS& tParams )
    {
        /*
         * gens is used to gereate a list of integers: 0, 1, 2, 3
         * that is as long as sizeof...( TP_ARGUMENTS )
         * Therefore we can use these ints in order to unpack the tuple tParams and pass it's content to execute
         */
        return this->callFunc( tParams, typename gens<sizeof...( TP_ARGUMENTS )>::type( ) );
    } // method

    virtual bool isFinished( ) const
    {
        if( !IS_VOLATILE )
            return false;
        return xFlags[ uiFinishedFlag ];
    } // method

    virtual bool requiresLock( ) const
    {
        return false;
    } // method
}; // class

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


class BasePledge
{
  public:
    /**
     * @brief Reset the pledge
     */
    virtual void reset( )
    {
        throw AnnotatedException( type_name( this ) + " did not implement reset" );
    } // method

    virtual std::shared_ptr<Container> getAsBaseType( )
    {
        throw AnnotatedException( type_name( this ) + " did not implement getAsBaseType" );
        return nullptr;
    } // method

    virtual void addSuccessor( BasePledge* pX )
    {
        throw AnnotatedException( type_name( this ) + " did not implement addSuccessor" );
    } // method

    virtual void removeSuccessor( BasePledge* pX )
    {
        throw AnnotatedException( type_name( this ) + " did not implement removeSuccessor" );
    } // method

    virtual bool hasVolatile( ) const
    {
        throw AnnotatedException( type_name( this ) + " did not implement hasVolatile" );
        return false;
    } // method

    virtual bool isFinished( ) const
    {
        throw AnnotatedException( type_name( this ) + " did not implement isFinished" );
        return false;
    } // method
    /**
     * @brief Gets the given pledges simultaneously.
     * @details
     * If bLoop is true the threads will keep going until all volatile modules are dry
     * if numThreads is not specified numThreads will be set to the amount of pledges given.
     * callback is called by one of the worker threads everytime the thread has finished one task.
     * If callback returns false the threadspool is destroyed.
     */
    static inline void simultaneousGet( std::vector<std::shared_ptr<BasePledge>> vPledges,
                                        std::function<bool( )> callback = []( ) { return true; },
                                        unsigned int numThreads = 0 )
    {
        if( numThreads == 0 )
            numThreads = (unsigned int)( vPledges.size( ) );

        /*
         * If there is a python module in the comp. graph we can only use one single thread.
         * This is due to a limitation in the Python interpreter. There can only ever be one
         * active Python thread.
         * So, here we check if there is such a module and set the number of threads to one if
         * necessary. @todo
         */
        // if( numThreads > 1 )
        //     for( std::shared_ptr<BasePledge> pPledge : vPledges )
        //         if( pPledge->hasPythonPledger( ) )
        //         {
        //             numThreads = 1;
        //             DEBUG( std::cout << "Detected python module. Cannot use more than one thread." << std::endl; )
        //             break;
        //         } // if

        std::mutex xExceptionMutex;
        std::string sExceptionMessageFromWorker;

        {
            /*
             * This variable is volatile so that every thread has to load if from memory each loop.
             * This way, if callback returns false thread 0 can set bContinue to false and all threads will stop after
             * their next iteration.
             */
            std::atomic_bool bContinue( true );
            // set up a threadpool
            ThreadPool xPool( numThreads );

            // enqueue a task that executes the comp. graph for each thread in the pool.
            for( std::shared_ptr<BasePledge> pPledge : vPledges )
            {
                xPool.enqueue(
                    [&callback, &bContinue, &xExceptionMutex, &sExceptionMessageFromWorker, &xPool](
                        size_t uiTid, std::shared_ptr<BasePledge> pPledge ) {
                        assert( pPledge != nullptr );

                        /*
                         * Set bLoop = true if there is a volatile module in the graph.
                         * This enables looping.
                         * If bLoop is not set to true here we only compute the result once,
                         * as all further computations would yield the same result.
                         */
                        bool bLoop = pPledge->hasVolatile( );

                        // execute the loop if bLoop == true or do just one iteration otherwise.
                        do
                        {
                            try
                            {
                                /*
                                 * Execute one iteration of the comp. graph.
                                 * Then set bLoop to false if the execution returned pEoFContainer.
                                 * If bLoop is false already (due to there beeing no
                                 * volatile module) then leave it as false.
                                 */
                                bLoop &= pPledge->getAsBaseType( ) != nullptr;

                                // this callback function can be used to set a progress bar
                                // for example.
                                if( uiTid == 0 )
                                    bContinue = callback( );
                            } // try
                            catch( const std::exception& rxException )
                            {
                                std::lock_guard<std::mutex> xExceptionGuard( xExceptionMutex );
                                if( !sExceptionMessageFromWorker.empty( ) )
                                    std::cerr
                                        << "Drop exception (different thread threw already): " << rxException.what( )
                                        << std::endl;
                                else
                                {
                                    sExceptionMessageFromWorker = rxException.what( );
                                    bContinue = false;
                                } // else
                                return;
                            } // catch
                            catch( ... )
                            {
                                std::lock_guard<std::mutex> xExceptionGuard( xExceptionMutex );
                                if( !sExceptionMessageFromWorker.empty( ) )
                                    std::cerr << "Drop unknown exception (different thread threw already)."
                                              << std::endl;
                                else
                                {
                                    sExceptionMessageFromWorker = "Unknown exception";
                                    bContinue = false;
                                } // else
                                return;
                            } // catch
                        } while( bLoop && bContinue );
                        DEBUG( std::cout << "Thread " << uiTid << " finished." << std::endl; )
                    }, // lambda
                    pPledge );
            } // for
            // wait for the pool to finish it's work
        } // scope xPool

        if( !sExceptionMessageFromWorker.empty( ) )
        {
            std::cerr << "Throw exception" << std::endl;
            throw std::runtime_error( sExceptionMessageFromWorker );
        } // if
    } // function
}; // class


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
// TP_TYPE must be subtype of Container!
template <class TP_TYPE, bool IS_VOLATILE = false, typename... TP_DEPENDENCIES> class Pledge : public BasePledge
{
  public:
    double execTime = 0;
    typedef TP_TYPE TP_CONTENT;
    /*
     * typename TP_DEPENDENCIES::TP_CONTENT...
     * this does expand the parameter pack and refers to TP_CONTENT within each individual parameter
     */
    typedef Module<TP_CONTENT, IS_VOLATILE, typename TP_DEPENDENCIES::TP_CONTENT...> TP_PLEDGER;
    friend TP_PLEDGER;
    typedef std::tuple<std::shared_ptr<TP_DEPENDENCIES>...> TP_PREDECESSORS;
    typedef std::tuple<std::shared_ptr<typename TP_DEPENDENCIES::TP_CONTENT>...> TP_INPUT;

  protected:
    std::shared_ptr<TP_PLEDGER> pPledger;
    // all pledges that have this pledge as predecessor
    std::vector<BasePledge*> vSuccessors;

  private:
    // content of the pledge
    std::shared_ptr<TP_CONTENT> pContent = nullptr;
    // this puts the shared_ptr around all elements of the variadic template individually
    TP_PREDECESSORS tPredecessors;

  public:
    void addSuccessor( BasePledge* pX )
    {
        vSuccessors.push_back( pX );
    } // method

    void removeSuccessor( BasePledge* pX )
    {
        // moves the element(s) in question to the end of the vector
        auto itNewEnd =
            std::remove_if( vSuccessors.begin( ), vSuccessors.end( ), [&]( BasePledge* pY ) { return pX == pY; } );
        // removes the element(s) in question
        vSuccessors.resize( itNewEnd - vSuccessors.begin( ) );
    } // method

  private:
    std::shared_ptr<std::mutex> pMutex;

    /**
     * @brief calls the get function of the predecessor IDX.
     * @details
     * Can be used with template loops.
     */
    template <size_t IDX> struct GetCaller
    {
        bool operator( )( TP_PREDECESSORS& tPredecessors, TP_INPUT& tInput )
        {
            std::get<IDX>( tInput ) = std::get<IDX>( tPredecessors )->get( );
            // handle EoF case
            if( std::get<IDX>( tInput ) == nullptr )
                return false;
            return true;
        } // operator
    }; // struct

    /**
     * @brief calls the hasVolatile function of the predecessor IDX.
     * @details
     * Can be used with template loops.
     */
    template <size_t IDX> struct NonVolatileCaller
    {
        bool operator( )( const TP_PREDECESSORS& tPredecessors )
        {
            return !std::get<IDX>( tPredecessors )->hasVolatile( );
        } // operator
    }; // struct

    /**
     * @brief calls the isFinished function of the predecessor IDX.
     * @details
     * Can be used with template loops.
     */
    template <size_t IDX> struct NotFinishedCaller
    {
        bool operator( )( const TP_PREDECESSORS& tPredecessors )
        {
            return !std::get<IDX>( tPredecessors )->isFinished( );
        } // operator
    }; // struct

    /**
     * @brief calls the addSuccessor function of the predecessor IDX.
     * @details
     * Can be used with template loops.
     */
    template <size_t IDX> struct AddToSuccessorCaller
    {
        bool operator( )( Pledge* pThis, TP_PREDECESSORS& tPredecessors )
        {
            std::get<IDX>( tPredecessors )->addSuccessor( pThis );
            return true;
        } // operator
    }; // struct

    /**
     * @brief calls the removeSuccessor function of the predecessor IDX.
     * @details
     * Can be used with template loops.
     */
    template <size_t IDX> struct RemFromSuccessorCaller
    {
        bool operator( )( Pledge* pThis, TP_PREDECESSORS& tPredecessors )
        {
            std::get<IDX>( tPredecessors )->removeSuccessor( pThis );
            return true;
        } // operator
    }; // struct

    /**
     * @brief Used to synchronize the execution of pledges in the comp. graph.
     * @details
     * Locks a mutex if this pledge can be reached from multiple leaves in the graph;
     * Does not lock otherwise.
     * In either case fDo is called.
     */
    template <typename FUNCTOR> inline std::shared_ptr<TP_CONTENT> lockIfNecessary( FUNCTOR&& fDo ) const
    {
        // if(vSuccessors.size() > 1) @todo @fixme this should be here
        if( pPledger->requiresLock( ) )
        {
            // multithreading is possible thus a guard is required here.
            // deadlock prevention is trivial,
            // since computational graphs are essentially trees.
            std::lock_guard<std::mutex> xGuard( *pMutex );
            return fDo( );
        } // if
        else
            return fDo( );
    } // function

  public:
    /**
     * @brief Create a new pledge with a cpp module responsible to fullfill it.
     * @details
     * This means that this Pledge can be automatically fullfilled by the given module.
     */
    Pledge( std::shared_ptr<TP_PLEDGER> pPledger, std::shared_ptr<TP_DEPENDENCIES>... tPredecessors )
        : pPledger( pPledger ), tPredecessors( tPredecessors... ), pMutex( new std::mutex )
    {
        // create a variable so that we can take a reference of it.
        Pledge* pThis = this;
        // add to the successor lists of the predecessors
        TemplateLoop<sizeof...( TP_DEPENDENCIES ), AddToSuccessorCaller>::iterate( pThis, this->tPredecessors );
    } // constructor

    /**
     * @brief Create a new pledge without a module giving the pledge.
     * @details
     * This means that this Pledge has to be fullfilled by calling set manually.
     */
    Pledge( )
    {
        static_assert( sizeof...( TP_DEPENDENCIES ) == 0, "cannot default constuct pledge with dependencies" );
    } // constructor

    /**
     * @brief this is required due to the use of mutex
     */
    Pledge( const Pledge& ) = delete; // copy constructor

    ~Pledge( )
    {
        // create a variable so that we can take a reference of it.
        Pledge* pThis = this;
        // remove this pledge from the successor lists of all predecessors
        TemplateLoop<sizeof...( TP_DEPENDENCIES ), RemFromSuccessorCaller>::iterate( pThis, this->tPredecessors );
    } // deconstructor

    inline const std::shared_ptr<TP_PLEDGER> getPledger( ) const
    {
        return pPledger;
    } // function

    /**
     * @brief Reset the pledge
     * @note override
     */
    virtual void reset( )
    {
        if( pContent == nullptr )
            return;
        pContent = nullptr;
        for( BasePledge* pSuccessor : vSuccessors )
            if( pSuccessor != nullptr )
                pSuccessor->reset( );
    } // function

    /**
     * @brief Manually fullfill this pledge.
     * @details
     * Invalidates the content of all following pledges.
     * @note It would be better to have a "constant" or "variable" class instead of this function.
     */
    inline void set( std::shared_ptr<TP_CONTENT> pC )
    {
        // improves runtime (mostly when resetting module).
        if( pContent == pC )
            return;
        pContent = pC;
        // the content of this pledge changed, therefore all following ones are invalidated and need to be reset.
        for( BasePledge* pSuccessor : vSuccessors )
            if( pSuccessor != nullptr )
                pSuccessor->reset( );
    } // function

    /**
     * @brief checks wether there is a volatile module upstream in the comp. graph.
     */
    virtual bool hasVolatile( ) const
    {
        if( IS_VOLATILE )
            return true;
        // double negation necessary here, so that the loop iterates until the first volatile module is found
        // (in this case the NonVolatileCaller returns false, therefore breaks the loop, and !false is returned )
        // or returns true at the end which then in turn is negated to become false (i.e. no volatile module found).
        return !TemplateLoop<sizeof...( TP_DEPENDENCIES ), NonVolatileCaller>::iterate( tPredecessors );
    } // method

    virtual bool isFinished( ) const
    {
        if( pPledger != nullptr && pPledger->isFinished( ) )
            return true;
        return !TemplateLoop<sizeof...( TP_DEPENDENCIES ), NotFinishedCaller>::iterate( tPredecessors );
    } // method

    /**
     * @brief Get the promised container.
     * @details
     * Checks weather the pledge has already been fullfilled for this call.
     * If necessary, uses the python or cpp module to fullfill the pledge,
     * and returns the respective container.
     */
    virtual std::shared_ptr<TP_CONTENT> get( )
    {
        if( pPledger == nullptr && pContent == nullptr )
            throw AnnotatedException( "No pledger known for unfulfilled pledge of type " + type_name( this ) +
                                      "; With container type: " + type_name( pContent ) );
        // in this case there is no need to execute again
        if( pContent != nullptr && !IS_VOLATILE )
            return pContent;

        // locks a mutex if this pledge can be reached from multiple leaves in the graph
        // does not lock otherwise...
        return lockIfNecessary( [&]( ) {
            TP_INPUT tInput;
            // execute all dependencies
            if( !TemplateLoop<sizeof...( TP_DEPENDENCIES ), GetCaller>::iterate( tPredecessors, tInput ) )
                /*
                 * if one dependency returns EoF then stop executing
                 * We have a volatile module that's dry.
                 * In such cases we cannot compute the next element and therefore set
                 * the content of this module to EoF as well.
                 */
                pContent = nullptr;
            // handle the EoF case for this module
            else if( IS_VOLATILE && pPledger->isFinished( ) )
                pContent = nullptr;
            else
            {
                auto timeStamp = std::chrono::system_clock::now( );

                // actually execute the module
                pContent = pPledger->executeTup( tInput );

                std::chrono::duration<double> duration = std::chrono::system_clock::now( ) - timeStamp;
                // increase the total executing time for this pledge
                execTime += duration.count( );

                // safety check
                if( pContent == nullptr )
                    throw AnnotatedException( "A module is not allowed to return nullpointers in execute; throw an "
                                              "exception instead or return an empty container! Module type:" +
                                              type_name( pPledger ) );
            } // else

            // return must be here and returned throught the lambda to avoid data race...
            return pContent;
        } // lambda
        ); // function call
    } // function

    virtual std::shared_ptr<Container> getAsBaseType( )
    {
        return std::dynamic_pointer_cast<Container>( this->get( ) );
    } // method
}; // class

/**
 * @brief generate a pledge from a module pointer and the dependencies.
 * @details
 * Does nothing more than calling the correct constructor of Pledge.
 * Use this function so that the compiler can figure out the templates itself.
 */
template <class TP_MODULE, class... TP_PLEDGES>
std::shared_ptr<Pledge<typename TP_MODULE::TP_RETURN, TP_MODULE::IS_VOLATILE, TP_PLEDGES...>>
promiseMe( std::shared_ptr<TP_MODULE> pModule, std::shared_ptr<TP_PLEDGES>... pPledges )
{
    return std::make_shared<Pledge<typename TP_MODULE::TP_RETURN, TP_MODULE::IS_VOLATILE, TP_PLEDGES...>>(
        pModule, pPledges... );
} // function

template <class TP_MODULE, class... TP_PLEDGES>
std::shared_ptr<Pledge<typename TP_MODULE::TP_RETURN, TP_MODULE::IS_VOLATILE, TP_PLEDGES...>>
apply( std::shared_ptr<TP_PLEDGES>... pPledges )
{
    auto pModule = std::make_shared<TP_MODULE>( );
    return std::make_shared<Pledge<typename TP_MODULE::TP_RETURN, TP_MODULE::IS_VOLATILE, TP_PLEDGES...>>(
        pModule, pPledges... );
} // function

template <class TP_CONTAINER, class... TP_ARGS> std::shared_ptr<Pledge<TP_CONTAINER>> makePledge( TP_ARGS&&... args )
{
    auto pRet = std::make_shared<Pledge<TP_CONTAINER>>( );
    pRet->set( std::make_shared<TP_CONTAINER>( args... ) );
    return pRet;
} // function

/**
 * @brief executes a volatile module over and over.
 * @note cannot deal with modules that have an input at the moment.
 */
template <typename TP_MODULE, typename... TP_CONSTR_PARAMS>
class Cyclic : public Module<typename TP_MODULE::TP_RETURN, true>
{
    const ParameterSetManager& rParameters;
    std::tuple<TP_CONSTR_PARAMS...> tParams;
    std::unique_ptr<TP_MODULE> pModule;

    template <int... S> void resetHelper( seq<S...> )
    {
        pModule = std::make_unique<TP_MODULE>( rParameters, std::get<S>( tParams )... );
        assert( !pModule->isFinished( ) );
    } // method

    void reset( )
    {
        resetHelper( typename gens<sizeof...( TP_CONSTR_PARAMS )>::type( ) );
    } // method

  public:
    Cyclic( const ParameterSetManager& rParameters, TP_CONSTR_PARAMS... params )
        : rParameters( rParameters ), tParams( params... )
    {
        static_assert( TP_MODULE::IS_VOLATILE == true, "can only declare Cyclic on volatile modules!" );
        reset( );
    } // constructor

    std::shared_ptr<typename TP_MODULE::TP_RETURN> execute( )
    {
        auto pRet = pModule->execute( );
        if( pModule->isFinished( ) )
            reset( );
        return pRet;
    } // method

    // override
    bool requiresLock( ) const
    {
        return true;
    } // method
}; // class


/**
 * @brief The module class that is exported to Python.
 */
template <bool IS_VOLATILE> class PyModule : public Module<Container, IS_VOLATILE, PyContainerVector>
{}; // class

#ifdef WITH_PYTHON
/**
 * @brief The module class that is exported to Python and allows python to define it's own modules.
 * @details
 * pybind11 calls this a trampoline class, since it redirects calls back to python.
 *
 */
template <bool IS_VOLATILE> class ModuleWrapperPyToCpp : public PyModule<IS_VOLATILE>
{
  public:
    ModuleWrapperPyToCpp( )
    {} // default constructor
    /**
     * @brief Execute the implemented algorithm.
     * @details
     * This makes execute overridable by python
     */
    std::shared_ptr<Container> execute( std::shared_ptr<PyContainerVector> pArgs ) override
    {
        PYBIND11_OVERLOAD( PYBIND11_TYPE( std::shared_ptr<Container> ),
                           PYBIND11_TYPE( PyModule<IS_VOLATILE> ),
                           execute,
                           pArgs ); // PYBIND11_OVERLOAD
    } // method
}; // class
#endif

/**
 * @brief Input pledge for all python modules
 * @details
 * This is used to break the typesaveness of the cpp code.
 * We cannot have the types for python due to the dynamic types used there.
 * This vector allows us to combine all pledges adn them send them as a single pledge (i.e. an object of this class).
 * Each python module then outputs a container that can in turn be pushed into one of these Pledge vectors.
 */
class PyPledgeVector : public Pledge<PyContainerVector>
{
  private:
    std::vector<std::shared_ptr<BasePledge>> vPledges;

  public:
    using Pledge<PyContainerVector>::Pledge;

    ~PyPledgeVector( )
    {
        for( std::shared_ptr<BasePledge> pPledge : vPledges )
            pPledge->removeSuccessor( this );
    } // deconstructor

    void append( std::shared_ptr<BasePledge> pX )
    {
        if( pX == nullptr )
            throw AnnotatedException( "Cannot use a nullpointer as pledge" );
        vPledges.emplace_back( pX );
        pX->addSuccessor( this );
    } // method

    // overrides the base function
    virtual std::shared_ptr<PyContainerVector> get( )
    {
        auto pRet = std::make_shared<PyContainerVector>( );
        for( std::shared_ptr<BasePledge> pPledge : vPledges )
        {
            pRet->push_back( pPledge->getAsBaseType( ) );
            // handle EoF case
            if( pRet->back( ) == nullptr )
                return nullptr;
        } // for
        return pRet;
    } // method

    // overrides the base function
    virtual void reset( )
    {
        // force reset on successors
        for( BasePledge* pSuccessor : vSuccessors )
            if( pSuccessor != nullptr )
                pSuccessor->reset( );
    } // function

    /**
     * @brief checks wether there is a volatile module upstream in the comp. graph.
     * @details
     * overrides the base function
     */
    virtual bool hasVolatile( ) const
    {
        for( std::shared_ptr<BasePledge> pPledge : vPledges )
            if( pPledge->hasVolatile( ) )
                return true;
        return false;
    } // method

    virtual bool isFinished( ) const
    {
        // if( Pledge<PyContainerVector>::isFinished( ) )
        //    return true;
        for( std::shared_ptr<BasePledge> pPledge : vPledges )
            if( pPledge->isFinished( ) )
                return true;
        return false;
    } // method

    inline void simultaneousGetPy( unsigned int numThreads = 0 )
    {
        BasePledge::simultaneousGet( vPledges, []( ) { return true; }, numThreads );
    } // method
}; // class

template <class TP_MODULE, typename... TP_CONSTR_PARAMS>
class ModuleWrapperCppToPy : public PyModule<TP_MODULE::IS_VOLATILE>
{
  public:
    TP_MODULE xModule;

  private:
    template <size_t IDX> struct TupleFiller
    {
        bool operator( )( typename TP_MODULE::TP_TUPLE_ARGS& tTup, PyContainerVector& vIn )
        {
            if( vIn[ IDX ] == nullptr )
                throw AnnotatedException( "Wrong type for module (" + type_name<TP_MODULE>( ) +
                                          ") input (parameter index: " + std::to_string( IDX ) +
                                          "). Expected: " + type_name( std::get<IDX>( tTup ) ) +
                                          " but got: nullptr of type " + type_name( vIn[ IDX ] ) );
            // convert the content in the vector to the element types of the tuple
            auto pCasted = std::dynamic_pointer_cast<
                typename std::tuple_element<IDX, typename TP_MODULE::TP_TUPLE_ARGS>::type::element_type>( vIn[ IDX ] );
            // if the cast is not possible we get a nullptr here
            if( pCasted == nullptr )
                throw AnnotatedException( "Wrong type for module (" + type_name<TP_MODULE>( ) +
                                          ") input (parameter index: " + std::to_string( IDX ) + "). Expected: " +
                                          type_name( std::get<IDX>( tTup ) ) + " but got: " + type_name( vIn[ IDX ] ) );
            // if the cast is successfull we assign the tuple element
            std::get<IDX>( tTup ) = pCasted;
            return true;
        } // operator
    }; // struct

  public:
    ModuleWrapperCppToPy( TP_CONSTR_PARAMS... params ) : xModule( params... )
    {} // constructor

    virtual std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<PyContainerVector> pIn )
    {
        if( pIn->size( ) != TP_MODULE::NUM_ARGUMENTS )
            throw AnnotatedException(
                type_name<TP_MODULE>( ) + "'s pyExecute was called with the wrong amount of parameters. Expected: " +
                std::to_string( TP_MODULE::NUM_ARGUMENTS ) + " but got: " + std::to_string( pIn->size( ) ) );

        typename TP_MODULE::TP_TUPLE_ARGS tTyped;
        TemplateLoop<TP_MODULE::NUM_ARGUMENTS, TupleFiller>::iterate( tTyped, *pIn );
        return xModule.executeTup( tTyped );
    } // method

    virtual bool isFinished( ) const
    {
        return xModule.isFinished( );
    } // method

    virtual bool requiresLock( ) const
    {
        return xModule.requiresLock( );
    } // method
}; // class


#ifdef WITH_PYTHON

template <class TP_MODULE>
void exportModule( pybind11::module& xPyModuleId, // pybind module variable
                   const std::string&& sName, // module name
                   std::function<void( py::class_<TP_MODULE>&& )> fExportMembers =
                       []( py::class_<TP_MODULE>&& ) {} // default lambda
)
{ // Export Class members
    fExportMembers( py::class_<TP_MODULE>( xPyModuleId,
                                           ( std::string( "__" ) + sName ).c_str( ) ) // Python class name
    );

    typedef ModuleWrapperCppToPy<TP_MODULE, const ParameterSetManager&> TP_TO_EXPORT;
    // Export MA-Module Class
    py::class_<TP_TO_EXPORT, // derived class
               PyModule<TP_MODULE::IS_VOLATILE>, // base class
               std::shared_ptr<TP_TO_EXPORT>> // reference holder
        ( xPyModuleId, sName.c_str( ) )
            .def( py::init<const ParameterSetManager&>( ) )
            .def_readonly( "cpp_module", &TP_TO_EXPORT::xModule );

    py::implicitly_convertible<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>>( );
} // function

template <class TP_MODULE, typename TP_CONSTR_PARAM_FIRST, typename... TP_CONSTR_PARAMS>
void exportModule( pybind11::module& xPyModuleId, // pybind module variable
                   const std::string&& sName, // module name
                   std::function<void( py::class_<TP_MODULE>&& )> fExportMembers =
                       []( py::class_<TP_MODULE>&& ) {} // default lambda
)
{
    typedef ModuleWrapperCppToPy<TP_MODULE, const ParameterSetManager&, TP_CONSTR_PARAM_FIRST, TP_CONSTR_PARAMS...>
        TP_TO_EXPORT;
    fExportMembers( py::class_<TP_MODULE>( xPyModuleId, ( std::string( "__" ) + sName ).c_str( ) ) );

    py::class_<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>, std::shared_ptr<TP_TO_EXPORT>>( xPyModuleId,
                                                                                               sName.c_str( ) )
        .def( py::init<const ParameterSetManager&, TP_CONSTR_PARAM_FIRST, TP_CONSTR_PARAMS...>( ) )
        .def_readonly( "cpp_module", &TP_TO_EXPORT::xModule );
    py::implicitly_convertible<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>>( );
} // function

template <class TP_MODULE, typename TP_CONSTR_PARAM_FIRST, typename... TP_CONSTR_PARAMS>
void exportModuleAlternateConstructor( pybind11::module& xPyModuleId, // pybind module variable
                                       const std::string&& sName // module name
)
{
    typedef ModuleWrapperCppToPy<TP_MODULE, const ParameterSetManager&, TP_CONSTR_PARAM_FIRST, TP_CONSTR_PARAMS...>
        TP_TO_EXPORT;

    py::class_<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>, std::shared_ptr<TP_TO_EXPORT>>( xPyModuleId,
                                                                                               sName.c_str( ) )
        .def( py::init<const ParameterSetManager&, TP_CONSTR_PARAM_FIRST, TP_CONSTR_PARAMS...>( ) )
        .def_readonly( "cpp_module", &TP_TO_EXPORT::xModule );
    py::implicitly_convertible<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>>( );
} // function

#endif

/*
 * @brief for testing purposes
 */
class Test : public Module<Container, true>
{
    size_t uiNumIt = 5;
    int iTest;
    std::mutex xMutex;

  public:
    Test( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<Container> execute( )
    {
        int iTestCpy = std::rand( );
        {
            std::lock_guard<std::mutex> xGuard( xMutex );
            iTest = iTestCpy;
        } // scope of xGuard
        std::this_thread::sleep_for( std::chrono::seconds( 1 ) );
        {
            std::lock_guard<std::mutex> xGuard( xMutex );
            if( iTest != iTestCpy )
                exit( EXIT_FAILURE );
        } // scope of xGuard
        uiNumIt--;
        if( uiNumIt == 0 )
            setFinished( );
        return std::make_shared<Container>( );
    } // method

    // override
    bool requiresLock( ) const
    {
        return true;
    } // method
}; // class

} // namespace libMA


#ifdef WITH_PYTHON
void exportModuleClass( py::module& rxPyModuleId );
#endif

#endif