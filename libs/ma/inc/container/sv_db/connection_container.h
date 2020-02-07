/**
 * @file connection_container.h
 * @brief implements libMA::ConnectionContainer a container that holds a connection to a database.
 * @author Markus Schmidt
 */
#pragma once

#include "container/container.h"

namespace libMA
{
/**
 * @brief container that holds a connection
 * @details
 * The purpose of this container is to hold a database connection, so that the connection
 * can be used in a computational graph.
 */
template <class DBCon> class ConnectionContainer : public Container
{
  public:
    std::shared_ptr<DBCon> pConnection;

    ConnectionContainer( std::shared_ptr<DBCon> pConnection ) : pConnection( pConnection )
    {} // constructor

}; // class

} // namespace libMA
