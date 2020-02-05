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
        : pInserter( TableType<DBCon>( pConnection ).getBulkInserter<BUFFER_SIZE>( ) ), iId( iId )
    {} // constructor

    /// @brief cannot be copied
    InserterContainer( const InserterContainer& rOther ) = delete; // delete copy constructor

    virtual void insert( std::shared_ptr<InsertTypes>... pArgs );

    virtual void close( )
    {
        pInserter.reset( );
    } // method
}; // class


template <template <typename T> typename InserterContainerType, typename DBCon, typename DBConInit,
          template <typename T> typename TableType, typename = typename TableType<DBConInit>::columnTypes>
class GetInserterContainerModule : public Module<InserterContainerType<DBCon>, false, ConnectionContainer<DBCon>>
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


template <typename InserterContainerType, typename = typename InserterContainerType::insertTypes_,
          typename... InsertTypes>
class InserterModule : public Module<Container, false, InserterContainerType, InsertTypes...>
{};

/**
 * @brief Wraps a jump inserter, so that it can become part of a computational graph.
 */
template <typename InserterContainerType, typename... InsertTypes>
class InserterModule<InserterContainerType, pack<InsertTypes...>, InsertTypes...>
{
  public:
    InserterModule( const ParameterSetManager& rParameters )
    {} // constructor

    std::shared_ptr<Container> execute( std::shared_ptr<InserterContainerType> pInserter,
                                        std::shared_ptr<InsertTypes>... pArgs )
    {
        pInserter->insert( pArgs... );

        return std::make_shared<Container>( );
        // end of score for xGuard
    } // method
}; // class

#ifdef WITH_PYTHON

template <class TP_MODULE, typename TP_CONSTR_PARAM_FIRST, typename... TP_CONSTR_PARAMS>
void exportModule2( pybind11::module& xPyModuleId, // pybind module variable
                    const std::string& sName, // module name
                    pack<TP_CONSTR_PARAMS...> // empty struct just to transport TP_CONSTR_PARAMS type
)
{
    typedef ModuleWrapperCppToPy<TP_MODULE, const ParameterSetManager&, TP_CONSTR_PARAM_FIRST, TP_CONSTR_PARAMS...>
        TP_TO_EXPORT;
    py::class_<TP_MODULE>( xPyModuleId, ( std::string( "__Get" ) + sName ).c_str( ) ).def( "id", &TP_MODULE::iId );

    py::class_<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>, std::shared_ptr<TP_TO_EXPORT>>(
        xPyModuleId, ( std::string( "Get" ) + sName ).c_str( ) )
        .def( py::init<const ParameterSetManager&, TP_CONSTR_PARAM_FIRST, TP_CONSTR_PARAMS...>( ) )
        .def( py::init<const ParameterSetManager&, int64_t>( ) )
        .def_readonly( "cpp_module", &TP_TO_EXPORT::xModule );

    py::implicitly_convertible<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>>( );
} // function

template <class TP_MODULE>
void exportModule2( pybind11::module& xPyModuleId, // pybind module variable
                    const std::string& sName // module name
)
{
    exportModule2<TP_MODULE, std::shared_ptr<typename TP_MODULE::DBConInit_>>(
        xPyModuleId, sName, typename TP_MODULE::TableType_::columnTypes( ) );
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
