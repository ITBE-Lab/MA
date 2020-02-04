/**
 * @file get_inserter_container_module.h
 * @brief templates for database inserters
 * @author Markus Schmidt
 */
#pragma once

#include "module/module.h"


template <typename DBCon, typename TableType, typename... InsertTypes> class InserterContainer : public Container
{
    const static size_t BUFFER_SIZE = 500;

  protected:
    using INSERTER_TYPE = decltype( TableType<DBCon>::template getBulkInserter<BUFFER_SIZE> );
    INSERTER_TYPE pInserter;

  public:
    int64_t iId;
    using InsertTypes = InsertTypes;

    InserterContainer( std::shared_ptr<DBCon> pConnection, int64_t iId )
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
 * @todo this and the GetJumpCallerContainerModule could be unified via a templated class...
 */
template <typename InserterContainerType, typename DBCon, typename DBConInit, typename TableType>
class GetInserterContainerModule : public Module<InserterContainerType<DBCon>, false, DBCon>
{
    int64_t iId;

  public:
    using InserterContainerType = InserterContainerType;
    using DBCon = DBCon;
    using DBConInit = DBConInit;
    using TableType = TableType;

    ///@brief creates a new inserter from given arguments
    GetInserterContainerModule( const ParameterSetManager& rParameters, std::shared_ptr<DBConInit> pConnection,
                                TableType::ColumnTypes... xArguments )
        : iId( TableType<DBConInit>( pConnection ).insert( xArguments... ) )
    {} // constructor

    ///@brief creates an inserter for an existing db element
    GetInserterContainerModule( const ParameterSetManager& rParameters, int64_t iId ) : iId( iId )
    {} // constructor

    /// @brief create a new jump inserter container
    std::shared_ptr<InserterContainerType<DBCon>> execute( std::shared_ptr<DBCon> pConnection )
    {
        return std::make_shared<InserterContainerType<DBCon>>( pConnection, iId );
    } // method
}; // class

/**
 * @brief Wraps a jump inserter, so that it can become part of a computational graph.
 */
template <typename DBCon, typename InserterContainerType>
class InserterModule : public Module<Container, false, InserterContainerType, InserterContainerType::InsertTypes>
{
  public:
    InserterModule( const ParameterSetManager& rParameters )
    {} // constructor

    /// @brief insert
    std::shared_ptr<Container> execute( std::shared_ptr<InserterContainerType> pInserter,
                                        std::shared_ptr<InserterContainerType::InsertTypes>... pArgs )
    {
        pInserter->insert( pArgs... );

        return std::make_shared<Container>( );
        // end of score for xGuard
    } // method
}; // class

#ifdef WITH_PYTHON

template <typename GetInserterContainerModuleType>
inline void exportInserterContainer( py::module& rxPyModuleId, const std::string& rName )
{
    // export the templated class
    py::class_<GetInserterContainerModuleType::InserterContainerType, Container,
               std::shared_ptr<GetInserterContainerModuleType::InserterContainerType>>( rxPyModuleId, rName )
        .def( py::init<std::shared_ptr<DBConType>, int64_t>( ) )
        .def( "insert", GetInserterContainerModuleType::InserterContainerType::insert )
        .def( "close", GetInserterContainerModuleType::InserterContainerType::close );

    exportModule<GetInserterContainerModuleType, std::shared_ptr<GetInserterContainerModuleType::DBConInit>,
                 GetInserterContainerModuleType::TableType::ColumnTypes...>(
        rxPyModuleId, "Get" + rName, []( auto&& x ) { x.def( "id", &GetInserterContainerModuleType::iId ) ); } );
} // function

#endif