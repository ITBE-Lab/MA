/**
 * @file pool_container.h
 * @brief implements libMS::PoolContainer a container that holds a connection to a database.
 * @author Markus Schmidt
 */
#ifdef WITH_DB
#pragma once

#include "db_con_pool.h"
#include "ms/container/container.h"
#include "ms/util/export.h"

namespace libMS
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
        : PoolContainer(
              uiPoolSize,
              json{{SCHEMA, {{NAME, sSchemaName}}},
#ifdef USE_PG
                   {CONNECTION, {{HOSTNAME, "localhost"}, {USER, "postgres"}, {PASSWORD, "admin"}, {PORT, 0}}}
#endif
#ifdef USE_MSQL
                   {CONNECTION, {{HOSTNAME, "localhost"}, {USER, "root"}, {PASSWORD, "admin"}, {PORT, 0}}}
#endif
              } )
    {} // constructor

}; // class

} // namespace libMS

#ifdef WITH_PYTHON

void exportPoolContainer( libMS::SubmoduleOrganizer& xOrganizer );

#endif
#endif
