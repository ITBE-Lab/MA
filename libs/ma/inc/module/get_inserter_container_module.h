/**
 * @file get_inserter_container_module.h
 * @brief templates for database inserters
 * @author Markus Schmidt
 */
#pragma once

#include "container/sv_db/connection_container.h"
#include "module/module.h"

namespace libMA
{

template <typename DBCon, template <typename T> typename TableType, typename... InsertTypes>
class InserterContainer : public Container
{
  public:
    using insertTypes_ = pack<InsertTypes...>;

    const static size_t BUFFER_SIZE = 500;
    using InserterType = typename TableType<DBCon>::template SQLBulkInserterType<BUFFER_SIZE>;
    std::shared_ptr<InserterType> pInserter;

    const int64_t iId;

    InserterContainer( std::shared_ptr<ConnectionContainer<DBCon>> pConnection, int64_t iId )
        : pInserter( std::make_shared<InserterType>( TableType<DBCon>( pConnection ) ) ), iId( iId )
    {} // constructor

    virtual void EXPORTED insert( std::shared_ptr<InsertTypes>... pArgs )
    {
        throw std::runtime_error("insert function of InserterContainer was not defined.");
    } // method

    virtual void close( )
    {
        pInserter.reset( );
    } // method
}; // class


template <template <typename T> typename InserterContainerType, typename DBCon, typename DBConInit,
          template <typename T> typename TableType, typename = typename TableType<DBConInit>::columnTypes>
class GetInserterContainerModule
{};

/**
 * @brief creates a call inserter container
 */
template <template <typename T> typename InserterContainerType, typename DBCon, typename DBConInit,
          template <typename T> typename TableType, typename... ColumnTypes>
class GetInserterContainerModule<InserterContainerType, DBCon, DBConInit, TableType, pack<ColumnTypes...>>
    : public Module<InserterContainerType<DBCon>, false, ConnectionContainer<DBCon>>
{
  public:
    const int64_t iId;

    using InserterContainerType_ = InserterContainerType<DBCon>;
    using DBCon_ = DBCon;
    using DBConInit_ = DBConInit;
    using TableType_ = TableType<DBConInit>;

    ///@brief creates a new inserter from given arguments
    GetInserterContainerModule( const ParameterSetManager& rParameters, std::shared_ptr<DBConInit> pConnection,
                                ColumnTypes... xArguments )
        : iId( TableType<DBConInit>( pConnection ).insert( xArguments... ) )
    {} // constructor

    ///@brief creates an inserter for an existing db element
    GetInserterContainerModule( const ParameterSetManager& rParameters, int64_t iId ) : iId( iId )
    {} // constructor

    /// @brief create a new jump inserter container
    std::shared_ptr<InserterContainerType<DBCon>> execute( std::shared_ptr<ConnectionContainer<DBCon>> pConnection )
    {
        return std::make_shared<InserterContainerType<DBCon>>( pConnection, iId );
    } // method
}; // class


template <typename InserterContainerType, typename = typename InserterContainerType::insertTypes_> class InserterModule
{};

/**
 * @brief Wraps a jump inserter, so that it can become part of a computational graph.
 */
template <typename InserterContainerType, typename... InsertTypes>
class InserterModule<InserterContainerType, pack<InsertTypes...>>
    : public Module<Container, false, InserterContainerType, InsertTypes...>
{
  public:
    InserterModule( const ParameterSetManager& rParameters )
    {} // constructor

    InserterModule( )
    {} // default constructor

    std::shared_ptr<Container> execute( std::shared_ptr<InserterContainerType> pInserter,
                                        std::shared_ptr<InsertTypes>... pArgs )
    {
        pInserter->insert( pArgs... );

        return std::make_shared<Container>( );
        // end of score for xGuard
    } // method
}; // class

#ifdef WITH_PYTHON

template <class TP_MODULE, typename... TP_CONSTR_PARAMS>
class ModuleWrapperCppToPy2 : public ModuleWrapperCppToPy<TP_MODULE, const ParameterSetManager&, int64_t>
{
  public:
    ModuleWrapperCppToPy2( const ParameterSetManager& xParams,
                           std::shared_ptr<typename TP_MODULE::DBConInit_> pConnection, TP_CONSTR_PARAMS... args )
        : ModuleWrapperCppToPy<TP_MODULE, const ParameterSetManager&, int64_t>(
              xParams, typename TP_MODULE::TableType_( pConnection ).insert( args... ) )
    {} // constructor

    ModuleWrapperCppToPy2( const ParameterSetManager& xParams, int64_t iId )
        : ModuleWrapperCppToPy<TP_MODULE, const ParameterSetManager&, int64_t>( xParams, iId )
    {} // constructor
}; // class

template <class TP_MODULE, typename... TP_CONSTR_PARAMS>
void exportModule2( pybind11::module& xPyModuleId, // pybind module variable
                    const std::string& sName, // module name
                    pack<TP_CONSTR_PARAMS...> // empty struct just to transport TP_CONSTR_PARAMS type
)
{
    // export the GetInserterContainerModule
    py::class_<TP_MODULE>( xPyModuleId, ( std::string( "__Get" ) + sName ).c_str( ) )
        .def_readonly( "id", &TP_MODULE::iId );

    // export GetInserterContainerModule python wrapper
    typedef ModuleWrapperCppToPy2<TP_MODULE, TP_CONSTR_PARAMS...> TP_TO_EXPORT;
    py::class_<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>, std::shared_ptr<TP_TO_EXPORT>>(
        xPyModuleId, ( std::string( "Get" ) + sName ).c_str( ) )
        .def( py::init<const ParameterSetManager&, std::shared_ptr<typename TP_MODULE::DBConInit_>,
                       TP_CONSTR_PARAMS...>( ) )
        .def( py::init<const ParameterSetManager&, int64_t>( ) )
        .def_readonly( "cpp_module", &TP_TO_EXPORT::xModule );

    py::implicitly_convertible<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>>( );
} // function

template <class TP_MODULE>
void exportModule2( pybind11::module& xPyModuleId, // pybind module variable
                    const std::string& sName // module name
)
{
    exportModule2<TP_MODULE>( xPyModuleId, sName, typename TP_MODULE::TableType_::columnTypes( ) );
};

template <typename GetInserterContainerModuleType>
inline void exportInserterContainer( py::module& rxPyModuleId, const std::string& rName )
{
    // export the templated class
    using InserterContainerType = typename GetInserterContainerModuleType::InserterContainerType_;
    py::class_<InserterContainerType, Container, std::shared_ptr<InserterContainerType>>( rxPyModuleId, rName.c_str( ) )
        .def( py::init<std::shared_ptr<ConnectionContainer<typename GetInserterContainerModuleType::DBCon_>>,
                       int64_t>( ) )
        .def( "insert", &InserterContainerType::insert )
        .def( "close", &InserterContainerType::close );

    exportModule2<GetInserterContainerModuleType>( rxPyModuleId, rName );
} // function

#endif

}; // namespace libMA
