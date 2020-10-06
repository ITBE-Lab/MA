
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include <cstdlib>
#include <iostream>

/* Regarding PostgreSQL:
 * Activation of spatial extensions for database:
 * Requried statement: "create extension postgis;"
 * In Ubuntu for DB postgres and user postgres: psql --user postgres -c 'create extension postgis' postgres
 */

// DB related includes
#include "sql_api.h"
#include "wkb_spatial.h"


auto jDBConfig = json{{SCHEMA,
                       {
                           {NAME, "testing"}
                           , { FLAGS, { DROP_ON_CLOSURE } }
                       }},
                      {CONNECTION, {{HOSTNAME, "127.0.0.1"}}}};

int main( void )
{
    return doNoExcept( [&] {
        std::vector<int64_t> vuiIds;
        std::vector<WKBUint64Rectangle> vxRectangles;

        auto pDatabase = std::make_shared<SQLDB<DBConn>>( jDBConfig );
        {
            SQLTableWithLibIncrPriKey<SQLDB<DBConn>, WKBUint64Rectangle> xTestTable(
                pDatabase,
                {{TABLE_NAME, "rectangle_test"},
                 {TABLE_COLUMNS,
                  {
                      // FIXME: PLACEHOLDER seems to unused - delete
                      {{COLUMN_NAME, "rectangle"}, {PLACEHOLDER, "ST_PolyFromWKB(?, 0)"}},
                  }}} );

            for( size_t uiI = 0; uiI < 100000; uiI++ )
            {
                vxRectangles.emplace_back( WKBUint64Rectangle( geom::Rectangle<uint64_t>( 1, 1, 1, 1 ) ) );
                vxRectangles.emplace_back( WKBUint64Rectangle( geom::Rectangle<uint64_t>( 1, 0, 9, 10 ) ) );
                vxRectangles.emplace_back( WKBUint64Rectangle( geom::Rectangle<uint64_t>( 2, 0, 8, 10 ) ) );
            } // for

            auto pBulkInserter = xTestTable.getBulkInserter<500>( );

            for( auto& xRect : vxRectangles )
            {
                // std::cout << xRect << std::endl;
                vuiIds.push_back( pBulkInserter->insert( xRect ) );
            } // for
        } // scope SQLQuery
        {
            SQLQuery<SQLDB<DBConn>, WKBUint64Rectangle> xQuery(
#ifdef POSTGRESQL
                // Regarding the explicit type cast '::geometry'
                // See:
                pDatabase, "SELECT ST_AsBinary(rectangle::geometry) FROM rectangle_test WHERE id = $1" );
#else // MySQL
                pDatabase, "SELECT ST_AsBinary(rectangle) FROM rectangle_test WHERE id = ?" );
#endif

            for( size_t uiI = 0; uiI < vxRectangles.size( ); uiI++ )
            {
                /*auto xRect =*/ xQuery.execAndGetNthCell<0>( vuiIds[ uiI ] ).getRect( );
                // std::cout << xRect << std::endl;
            } // for
        } // scope SQLQuery
    } ); // doNoExcept
} // main function
