/**
 * @file get_inserter_container_module.h
 * @brief templates for database inserters
 * @author Markus Schmidt
 */
#pragma once

#include "container/sv_db/pool_container.h"
#include "module/module.h"

namespace libMA
{

/**
 * @brief container that holds the table that is inserted to
 * @details
 * this is part of a 3 class structure that should be used together:
 * - InserterContainer or BulkInserterContainer (container that holds the table that is inserted to)
 * - GetInserterContainerModule (module that generates an InserterContainer/BulkInserterContainer)
 * - InserterModule (module that calls the insert function of the InserterContainer/BulkInserterContainer)
 * The purpose of a class that extends this template is to implement an insert function that takes arguments of
 * InsertTypes... type and inserts them into a table of type TableType.
 * This class holds the table pTable for this purpose.
 * This inserter is intended for e.g. inserting reads, where each inserted row has a foreign key if belongs to.
 * E.g. the sequencer id for the read. This foreign key is saved in iId and available to the extending class.
 */
template <typename DBCon, typename TableType, typename... InsertTypes>
class AbstractInserterContainer : public Container
{
  public:
    using InsertTypesForw = TypePack<InsertTypes...>;
    using DBConForwarded = DBCon;

    const int iConnectionId;

    const int64_t iId;

    std::shared_ptr<TableType> pInserter;

  private:
    // needs to be below pTable, so that the transactions destructor is called first
    typename DBCon::sharedGuardedTrxnType pTransaction;
    // needs to be below pTransaction @todo discuss arne: sharedGuardedTrxnType does not keep DB connection alive
    std::shared_ptr<DBCon> pConnection;

  public:
    AbstractInserterContainer(
        std::tuple<typename DBCon::sharedGuardedTrxnType, int, std::shared_ptr<TableType>, std::shared_ptr<DBCon>> xFromRun, int64_t iId )
        : iConnectionId( std::get<1>( xFromRun ) ),
          iId( iId ),
          pInserter( std::get<2>( xFromRun ) ),
          pTransaction( std::get<0>( xFromRun ) ),
          pConnection( std::get<3>( xFromRun ) )
    {}

    /**
     * @brief insert an element into the table TableType
     * @details
     * Whenever you inherit from this class your job is to implement this function.
     */
    virtual void EXPORTED insert( std::shared_ptr<InsertTypes>... pArgs )
    {
        throw std::runtime_error( "insert function of AbstractInserterContainer was not defined." );
    } // method

  private:
    static void insert_helper( std::shared_ptr<DBCon> pConnection, AbstractInserterContainer* pThis,
                               std::shared_ptr<InsertTypes>... pArgs )
    {
        pThis->insert( pArgs... );
    } // method

  public:
    void pool_save_insert( std::shared_ptr<PoolContainer<DBCon>> pPool, std::shared_ptr<InsertTypes>... pArgs )
    {
        pPool->xPool.run( iConnectionId, &AbstractInserterContainer::insert_helper, this, pArgs... );
    } // method

    /**
     * @brief destruct the table / close the bulk inserter
     * @details
     * Once this is called any following calls to insert result in undefined behaviour.
     * In a c++ application this should never be called, since the deconstructor of this class fullfills the same
     * purpose.
     * Use this from python in order to make sure all elements are inserted.
     */
    virtual void close( )
    {
        pInserter.reset( );
        pTransaction.reset( );
        pConnection.reset();
    } // method
}; // class


/**
 * @brief container that holds a regular inserter for a table
 * @details
 * @see AbstractInserterContainer for docu
 */
template <typename DBCon, template <typename T> typename TableType, typename... InsertTypes>
class InserterContainer : public AbstractInserterContainer<DBCon, TableType<DBCon>, InsertTypes...>
{
  public:
    InserterContainer( std::shared_ptr<PoolContainer<DBCon>> pPool, int64_t iId )
        : AbstractInserterContainer<DBCon, TableType<DBCon>, InsertTypes...>(
              pPool->xPool.run( pPool->xPool.getDedicatedConId( ),
                                [this]( auto pConnection ) //
                                {
                                    return std::make_tuple( nullptr /*pConnection->sharedGuardedTrxn( )*/,
                                                            (int)pConnection->getTaskId( ),
                                                            std::make_shared<TableType<DBCon>>( pConnection ),
                                                            pConnection );
                                } ),
              iId )
    {} // constructor

}; // class

/// @brief a short typename for the bulk inserter type
template <typename TableType> using BulkInserterType = typename TableType::template SQLBulkInserterType<500>;
/**
 * @brief container that holds a bulk inserter for a table
 * @details
 * @see AbstractInserterContainer
 * This class does the same thing as InserterContainer. But instead of supplying the table it supplies a bulk inserter
 * to the table via pInserter.
 */
template <typename DBCon, template <typename T> typename TableType, typename... InsertTypes>
class BulkInserterContainer
    : public AbstractInserterContainer<DBCon, BulkInserterType<TableType<DBCon>>, InsertTypes...>
{
  public:
    BulkInserterContainer( std::shared_ptr<PoolContainer<DBCon>> pPool, int64_t iId )
        // here we create the bulk inserter. This forces us to construct an object of the table that is inserter to
        // This construction makes sure that the table exists in the database.
        : AbstractInserterContainer<DBCon, BulkInserterType<TableType<DBCon>>, InsertTypes...>(
              pPool->xPool.run( pPool->xPool.getDedicatedConId( ),
                                [this]( auto pConnection ) //
                                {
                                    return std::make_tuple(
                                        nullptr /*pConnection->sharedGuardedTrxn( )*/,
                                        (int)pConnection->getTaskId( ),
                                        TableType<DBCon>( pConnection ).template getBulkInserter<500>( ),
                                        pConnection );
                                } ),
              iId )
    {} // constructor

}; // class


/**
 * @brief Module that generates an InserterContainer/BulkInserterContainer
 * @details
 * this is part of a 3 class structure that should be used together:
 * - InserterContainer or BulkInserterContainer (container that holds the table that is inserted to)
 * - GetInserterContainerModule (module that generates an InserterContainer/BulkInserterContainer)
 * - InserterModule (module that calls the insert function of the InserterContainer/BulkInserterContainer)
 * This module is used in the computational graph in order to create a InserterContainer/BulkInserterContainer from
 * a PoolContainer.
 * InserterContainer's documentation says: "each inserted row has a foreign key if belongs to. E.g. the sequencer id for
 * the read."
 * This module has two constructors: one where the value of the foreign key can be given explicitly.
 * The other inserts a new row and therefore creates a new value for the foreign key.
 */
template <template <typename T> typename InserterContainerType, typename DBCon, typename DBConInit,
          template <typename T> typename TableType, typename = typename TableType<DBConInit>::ColTypesForw>
class GetInserterContainerModule
{}; // class

template <template <typename T> typename InserterContainerType, typename DBCon, typename DBConInit,
          // we want to extract ColumnTypes from TableType::ColTypesForw
          // in order to achieve this we use TypePack<ColumnTypes...> in combination with the template above
          // the above template makes it so, that the compiler can infer ColumnTypes on it's own.
          template <typename T> typename TableType, typename... ColumnTypes>
class GetInserterContainerModule<InserterContainerType, DBCon, DBConInit, TableType, TypePack<ColumnTypes...>>
    : public Module<InserterContainerType<DBCon>, false, PoolContainer<DBCon>>
{
  public:
    const int64_t iId;

    // export the template types so we can access them from outside
    // there really should be a way to just get a template parameter from a type in standard c++.   -.-
    // @todo InserterContainerTypeForw -> ForwardedInserterContainerType
    using InserterContainerTypeForw = InserterContainerType<DBCon>;
    using DBConForw = DBCon;
    using DBConInitForw = DBConInit;
    using TableTypeForw = TableType<DBConInit>;

    ///@brief creates a new inserter from given arguments
    GetInserterContainerModule( const ParameterSetManager& rParameters, std::shared_ptr<DBConInit> pConnection,
                                ColumnTypes... xArguments )
        : iId( TableType<DBConInit>( pConnection ).insert( xArguments... ) )
    {} // constructor

    ///@brief creates an inserter for an existing db element
    GetInserterContainerModule( const ParameterSetManager& rParameters, int64_t iId ) : iId( iId )
    {} // constructor

    std::shared_ptr<InserterContainerType<DBCon>> execute( std::shared_ptr<PoolContainer<DBCon>> pPool )
    {
        return std::make_shared<InserterContainerType<DBCon>>( pPool, iId );
    } // method
}; // class

/**
 * @brief Module that calls the insert function of the InserterContainer/BulkInserterContainer
 * @details
 * this is part of a 3 class structure that should be used together:
 * - InserterContainer or BulkInserterContainer (container that holds the table that is inserted to)
 * - GetInserterContainerModule (module that generates an InserterContainer/BulkInserterContainer)
 * - InserterModule (module that calls the insert function of the InserterContainer/BulkInserterContainer)
 * The actual insertion into a table is done via InserterContainer/BulkInserterContainer. However, both these
 * classes are containers. In order to use them in a comp. graph, we need a module to call them.
 * This template provides a module that does nothing else than calling the insert function of the inserter containers.
 * In order to do so it extracts the parameters of the insert function from InserterContainerType::InsertTypesForw.
 * When using this in a computational graph, supply the inserter container pledge as first container, followed by all
 * containers, you specified in the InserterContainer::InsertTypes... template.
 */
template <typename InserterContainerType, typename = typename InserterContainerType::InsertTypesForw>
class InserterModule
{}; // class

/**
 * @brief Wraps a jump inserter, so that it can become part of a computational graph.
 */
template <typename InserterContainerType, typename... InsertTypes>
class InserterModule<InserterContainerType, TypePack<InsertTypes...>>
    // we want to extract InsertTypes from InserterContainerType::InsertTypesForw
    // in order to achieve this we use TypePack<InsertTypes...> in combination with the template above
    // the above template makes it so, that the compiler can inferr InsertTypes on it's own.
    : public Module<Container, false, InserterContainerType,
                    PoolContainer<typename InserterContainerType::DBConForwarded>, InsertTypes...>
{
  public:
    InserterModule( const ParameterSetManager& rParameters )
    {} // constructor

    InserterModule( )
    {} // default constructor

    std::shared_ptr<Container> execute( std::shared_ptr<InserterContainerType> pInserter,
                                        std::shared_ptr<PoolContainer<typename InserterContainerType::DBConForwarded>>
                                            pPool,
                                        std::shared_ptr<InsertTypes>... pArgs )
    {
        pInserter->pool_save_insert( pPool, pArgs... );

        return std::make_shared<Container>( );
    } // method
}; // class

#ifdef WITH_PYTHON

/**
 * @brief helper for exportInserterContainer
 * @see exportInserterContainer
 */
template <class TP_MODULE, typename... TP_CONSTR_PARAMS>
class ModuleWrapperCppToPy2 : public ModuleWrapperCppToPy<TP_MODULE, const ParameterSetManager&, int64_t>
{
  public:
    ModuleWrapperCppToPy2( const ParameterSetManager& xParams,
                           std::shared_ptr<typename TP_MODULE::DBConInitForw> pConnection, TP_CONSTR_PARAMS... args )
        : ModuleWrapperCppToPy<TP_MODULE, const ParameterSetManager&, int64_t>(
              xParams,
              // in order to reuse the code from ModuleWrapperCppToPy (<- no 2 at the end) we need to reduce both
              // constructors down to a single one. This can be done by reducing the row insert constructor to the one
              // that simply takes the foreign key.
              // However, in order to preserve functionality we need to insert the the row here.
              typename TP_MODULE::TableTypeForw( pConnection ).insert( args... ) )
    {} // constructor

    ModuleWrapperCppToPy2( const ParameterSetManager& xParams, int64_t iId )
        // simply redirect to the ModuleWrapperCppToPy (<- no 2 at the end) constructor
        : ModuleWrapperCppToPy<TP_MODULE, const ParameterSetManager&, int64_t>( xParams, iId )
    {} // constructor
}; // class

/**
 * @brief helper for exportInserterContainer
 * @see exportInserterContainer
 */
template <class TP_MODULE, typename... TP_CONSTR_PARAMS>
void exportModule2( pybind11::module& xPyModuleId, // pybind module variable
                    const std::string& sName, // module name
                    TypePack<TP_CONSTR_PARAMS...> // empty struct just to transport TP_CONSTR_PARAMS type
)
{
    // export the GetInserterContainerModule
    py::class_<TP_MODULE>( xPyModuleId, ( std::string( "__Get" ) + sName ).c_str( ) )
        .def_readonly( "id", &TP_MODULE::iId );

    // export GetInserterContainerModule python wrapper
    // ModuleWrapperCppToPy2 is needed in order to allow two constructors
    typedef ModuleWrapperCppToPy2<TP_MODULE, TP_CONSTR_PARAMS...> TP_TO_EXPORT;
    py::class_<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>, std::shared_ptr<TP_TO_EXPORT>>(
        xPyModuleId, ( std::string( "Get" ) + sName ).c_str( ) )
        // export the constructor that generates the foreign key by creating a new row
        .def( py::init<const ParameterSetManager&, std::shared_ptr<typename TP_MODULE::DBConInitForw>,
                       TP_CONSTR_PARAMS...>( ) )
        // export the constructor that simply takes the foreign key as argument
        .def( py::init<const ParameterSetManager&, int64_t>( ) )
        .def_readonly( "cpp_module", &TP_TO_EXPORT::xModule );

    py::implicitly_convertible<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>>( );
} // function

/**
 * @brief helper for exportInserterContainer
 * @see exportInserterContainer
 */
template <typename InserterContainerType, typename DBCon, typename... InsertTypes>
void exportHelper( pybind11::module& xPyModuleId, // pybind module variable
                   const std::string& sName, // module name
                   TypePack<InsertTypes...> // empty struct just to transport InsertTypes type
)
{
    py::class_<InserterContainerType, Container, std::shared_ptr<InserterContainerType>>( xPyModuleId, sName.c_str( ) )
        // inserter container should never be constructed directly in python
        // .def( py::init<std::shared_ptr<PoolContainer<DBCon>>, int64_t>( ) )
        .def( "insert",
              // insert might be overloaded in some child classes
              // in order to pick the correct insert we cast it.
              ( void ( InserterContainerType::* )( std::shared_ptr<InsertTypes>... ) ) & InserterContainerType::insert )
        .def( "close", &InserterContainerType::close );
} // function

/**
 * @brief exports a child of GetInserterContainerModule and the InserterContainerType from it's template
 * @details
 * The InserterContainer is experted under the name rName.
 * For the GetInserterContainerModule "Get" is used as prefix.
 */
template <typename GetInserterContainerModuleType>
inline void exportInserterContainer( py::module& rxPyModuleId, const std::string& rName )
{
    // export the templated class
    exportHelper<typename GetInserterContainerModuleType::InserterContainerTypeForw,
                 typename GetInserterContainerModuleType::DBConForw>(
        rxPyModuleId, rName, typename GetInserterContainerModuleType::InserterContainerTypeForw::InsertTypesForw( ) );

    // TP_MODULE::TableTypeForw::ColTypesForw( ) is passed as last argument in order to transprot the constructor
    // parameters ColTypesForw is a completely empty structure
    exportModule2<GetInserterContainerModuleType>(
        rxPyModuleId, rName, typename GetInserterContainerModuleType::TableTypeForw::ColTypesForw( ) );
} // function

#endif

}; // namespace libMA
