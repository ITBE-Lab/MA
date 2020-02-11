// Remarks: MySQL tests does not work in Visual Studio 2017.
// Reason: https://forums.mysql.com/read.php?168,379738,379738

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include <cstdlib>
#include <iostream>

// DB related includes
#include "MySQL_con.h" // MySQL connector
#include "common.h"
#include "wkb_spatial.h"

auto jDBConfig =
    json{ { SCHEMA, "testing" },
          { CONNECTION, { { HOSTNAME, "127.0.0.1" }, { USER, "root" }, { PASSWORD, "admin" }, { PORT, 0 } } } };

int main( void )
{
    return doNoExcept( [ & ] {
        std::vector<int64_t> vuiIds;
        std::vector<WKBUint64Rectangle> vxRectangles;

        {
            auto pDatabase = std::make_shared<SQLDB<MySQLConDB>>( jDBConfig );
            SQLTableWithAutoPriKey<SQLDB<MySQLConDB>, WKBUint64Rectangle> xTestTable(
                pDatabase,
                { { TABLE_NAME, "rectangle_test" },
                  { TABLE_COLUMNS,
                    {
                        { { COLUMN_NAME, "rectangle" }, { PLACEHOLDER, "ST_PolyFromWKB(?, 0)" } },
                    } } } );

        vxRectangles.emplace_back( WKBUint64Rectangle( geom::Rectangle<uint64_t>( 1, 1, 1, 1 ) ) );
        vxRectangles.emplace_back( WKBUint64Rectangle( geom::Rectangle<uint64_t>( 1, 0, 9, 10 ) ) );
        vxRectangles.emplace_back( WKBUint64Rectangle( geom::Rectangle<uint64_t>( 2, 0, 8, 10 ) ) );

            for( auto& xRect : vxRectangles )
            {
                std::cout << xRect << std::endl;
                vuiIds.push_back( xTestTable.insert( xRect ) );
            } // for
        } // scope database
        {
            auto pDatabase = std::make_shared<SQLDB<MySQLConDB>>( jDBConfig );
            SQLQuery<SQLDB<MySQLConDB>, WKBUint64Rectangle> xQuery(
                pDatabase, "SELECT ST_AsBinary(rectangle) FROM rectangle_test WHERE id = ?" );

            for( size_t uiI = 0; uiI < vxRectangles.size( ); uiI++ )
            {
                auto xRect = xQuery.execAndGetNthCell<0>( vuiIds[ uiI ] ).getRect( );
                std::cout << xRect << std::endl;
            } // for
        } // scope database
    } ); // doNoExcept
} // main function
