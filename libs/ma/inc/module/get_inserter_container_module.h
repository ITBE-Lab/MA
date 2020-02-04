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
    const static size_t BUFFER_SIZE = 500;
    using InserterType = typename TableType<DBCon>::template SQLBulkInserterType<BUFFER_SIZE>;
    std::shared_ptr<InserterType> pInserter;

  public:
    const int64_t iId;

    InserterContainer( std::shared_ptr<ConnectionContainer<DBCon>> pConnection, int64_t iId )
        : pInserter( TableType<DBCon>( pConnection ).template getBulkInserter<BUFFER_SIZE>( ) ), iId( iId )
    {} // constructor

    /// @brief cannot be copied
    InserterContainer( const InserterContainer& rOther ) = delete; // delete copy constructor

    virtual void insert( std::shared_ptr<InsertTypes>... pArgs );

    virtual void close( )
    {
        pInserter->reset( );
    } // method
}; // class

/**
 * @brief creates a call inserter container
 * @details
 * @todo ArgTypes can be removed here using std::tuple and the gens structure defined in module.h
 */
template <template <typename T> typename InserterContainerType, typename DBCon, typename DBConInit,
          template <typename T> typename TableType, typename... ArgTypes>
class GetInserterContainerModule : public Module<InserterContainerType<DBCon>, false, ConnectionContainer<DBCon>>
{
    const int64_t iId;

  public:
    using InserterContainerType_ = InserterContainerType<DBCon>;
    using DBCon_ = DBCon;
    using DBConInit_ = DBConInit;

    ///@brief creates a new inserter from given arguments
    GetInserterContainerModule( const ParameterSetManager& rParameters, std::shared_ptr<DBConInit> pConnection,
                                ArgTypes... xArguments )
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

/**
 * @brief Wraps a jump inserter, so that it can become part of a computational graph.
 */
template <typename DBCon, typename InserterContainerType, typename... ArgTypes>
class InserterModule : public Module<Container, false, InserterContainerType, ArgTypes...>
{
  public:
    InserterModule( const ParameterSetManager& rParameters )
    {} // constructor

    /// @brief insert
    std::shared_ptr<Container> execute( std::shared_ptr<InserterContainerType> pInserter,
                                        std::shared_ptr<ArgTypes>... pArgs )
    {
        pInserter->insert( pArgs... );

        return std::make_shared<Container>( );
        // end of score for xGuard
    } // method
}; // class

#ifdef WITH_PYTHON

template <class TP_MODULE, typename TP_CONSTR_PARAM_FIRST, typename... TP_CONSTR_PARAMS>
void exportModule2( pybind11::module& xPyModuleId, // pybind module variable
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
        .def( py::init<const ParameterSetManager&, int64_t>( ) )
        .def_readonly( "cpp_module", &TP_TO_EXPORT::xModule );
    py::implicitly_convertible<TP_TO_EXPORT, PyModule<TP_MODULE::IS_VOLATILE>>( );
} // function

template <typename GetInserterContainerModuleType, typename... ArgTypes>
inline void exportInserterContainer( py::module& rxPyModuleId, const std::string& rName )
{
    // export the templated class
    py::class_<typename GetInserterContainerModuleType::InserterContainerType_, Container,
               std::shared_ptr<typename GetInserterContainerModuleType::InserterContainerType_>>( rxPyModuleId, rName )
        .def( py::init<std::shared_ptr<typename GetInserterContainerModuleType::DBConType_>, int64_t>( ) )
        .def( "insert", GetInserterContainerModuleType::InserterContainerType::insert )
        .def( "close", GetInserterContainerModuleType::InserterContainerType::close );

    exportModule2<GetInserterContainerModuleType, std::shared_ptr<typename GetInserterContainerModuleType::DBConInit_>,
                  ArgTypes...>( rxPyModuleId, "Get" + rName,
                                []( auto&& x ) { x.def( "id", &GetInserterContainerModuleType::iId ); } );
} // function

#endif

}; // namespace libMA
