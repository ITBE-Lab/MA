/**
 * @file pool_container.h
 * @brief implements libMA::PoolContainer a container that holds a connection to a database.
 * @author Markus Schmidt
 */
#pragma once

#include "container/container.h"
#include "db_con_pool.h"

namespace libMA
{
/**
 * @brief container that holds a SQLDBConPool
 * @details
 * The purpose of this container is to hold a database pool, so that the pool
 * can be used in a computational graph.
 */
template <class DBCon> class PoolContainer : public Container
{
  public:
    SQLDBConPool<typename DBCon::DBImplForwarded> xPool;

    PoolContainer( size_t uiPoolSize = 1, const json& jDBConData = json{} ) : xPool( uiPoolSize, jDBConData )
    {} // constructor

    PoolContainer( size_t uiPoolSize, std::string sSchemaName )
        : PoolContainer( uiPoolSize, json{{"SCHEMA", sSchemaName}} )
    {} // constructor

}; // class

} // namespace libMA

#ifdef WITH_PYTHON

void exportPoolContainer( py::module& rxPyModuleId );

#endif
