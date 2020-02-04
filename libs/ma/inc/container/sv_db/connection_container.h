/**
 * @file connection_container.h
 * @brief implements libMA::ConnectionContainer a container that holds a connection to a database.
 * @author Markus Schmidt
 */
#pragma once

namespace libMA
{
    template <typename DBCon> class ConnectionContainer : public DBCon, public Container
    {
        using DBCon::DBCon;
    }; // class
} // namespace libMA
