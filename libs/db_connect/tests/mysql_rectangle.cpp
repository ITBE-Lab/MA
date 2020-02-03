#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include <cstdlib>
#include <iostream>

// DB related includes
#include "common.h"
// #include "geom.h"
#include "MySQL_con.h" // MySQL connector
#include "wkb_spatial.h"


std::ostream& operator<<( std::ostream& xOS, const geomUtil::Interval<uint64_t>& xInterval )
{
    xOS << "(" << xInterval.start( ) << ", " << xInterval.end( ) << "] ";
    return xOS;
} // operator

std::ostream& operator<<( std::ostream& xOS, const WKBUint64Rectangle& xRectWKB )
{
    const auto xRect = xRectWKB.getRect( );
	xOS << std::dec << "Rectangle: x: " << xRect.xXAxis << " y: " << xRect.xYAxis << std::endl;
    xOS << "WKBPolygon: ";
    for( auto uiI : xRectWKB.aData )
        xOS << std::hex << (int)uiI << " ";
    xOS << std::endl;
    return xOS;
} // operator

int main( void )
{
    std::vector<int64_t> vuiIds;
    std::vector<WKBUint64Rectangle> vxRectangles;
    {
        auto pDatabase = std::make_shared<SQLDB<MySQLConDB>>( );

        SQLTableWithAutoPriKey<SQLDB<MySQLConDB>, WKBUint64Rectangle> xTestTable(
            pDatabase,
            { { TABLE_NAME, "rectangle_test" },
              { TABLE_COLUMNS,
                {
                    { { COLUMN_NAME, "rectangle" }, { PLACEHOLDER, "ST_PolyFromWKB(?, 0)" } },
                } } } );

        vxRectangles.emplace_back( WKBUint64Rectangle( geomUtil::Rectangle<uint64_t>( 1, 1, 1, 1 ) ) );
        vxRectangles.emplace_back( WKBUint64Rectangle( geomUtil::Rectangle<uint64_t>( 1, 0, 9, 10 ) ) );
        vxRectangles.emplace_back( WKBUint64Rectangle( geomUtil::Rectangle<uint64_t>( 2, 0, 8, 10 ) ) );

        for( auto& xRect : vxRectangles )
        {
            std::cout << xRect << std::endl;
            vuiIds.push_back( xTestTable.insert( xRect ) );
        }
    }

    {
        auto pDatabase = std::make_shared<SQLDB<MySQLConDB>>( );
        SQLQuery<SQLDB<MySQLConDB>, WKBUint64Rectangle> xQuery(
            pDatabase, "SELECT ST_AsBinary(rectangle) FROM rectangle_test WHERE id = ?" );

        for( size_t uiI = 0; uiI < vxRectangles.size( ); uiI++ )
        {
            auto xRect = xQuery.execAndGetNthCell<0>( vuiIds[ uiI ] ).getRect( );
            std::cout << xRect << std::endl;
        }
    }

    return EXIT_SUCCESS;
} /// main function
