#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "MySQL_con.h" // MySQL connector
#include "common.h"
#include "geom.h"
#include "spatial.h"
#include <cstdlib>
#include <iostream>


std::ostream& operator<<( std::ostream& xOS, const geomUtil::Interval<uint64_t>& xInterval )
{
    xOS << "(" << xInterval.start( ) << ", " << xInterval.end( ) << "] ";
    return xOS;
}

std::ostream& operator<<( std::ostream& xOS, const geomUtil::Rectangle<uint64_t>& xRect )
{
    xOS << "Rectangle: x: " << xRect.xXAxis << " y: " << xRect.xYAxis << std::endl;
    return xOS;
}

typedef geomUtil::WKBPolygon<geomUtil::Rectangle<uint64_t>::uiSizeWKB> WKBRectangle;
std::ostream& operator<<( std::ostream& xOS, const MyRectangle& xRect )
{
    xOS << "WKBPolygon: ";
    for( auto uiI : xRect.aData )
        xOS << std::hex << (int)uiI << " ";
    xOS << std::endl;
    return xOS;
}

int main( void )
{
    std::vector<int64_t> vuiIds;
    std::vector<WKBRectangle> vxRectanglesWKB;
    std::vector<geomUtil::Rectangle<uint64_t>> vxRectangles;

    vxRectangles.emplace_back( 1, 1, 1, 1 );
    vxRectangles.emplace_back( 1, 0, 9, 10 );
    vxRectangles.emplace_back( 2, 0, 8, 10 );

    {
        auto pDatabase = std::make_shared<SQLDB<MySQLConDB>>( );

        SQLTableWithAutoPriKey<SQLDB<MySQLConDB>, WKBRectangle> xTestTable(
            pDatabase,
            {{TABLE_NAME, "rectangle_test"},
             {TABLE_COLUMNS,
              {
                  {{COLUMN_NAME, "rectangle"}, {PLACEHOLDER, "ST_PolyFromWKB(?, 0)"}},
              }}} );


        for( auto& xRect : vxRectangles )
        {
            vxRectanglesWKB.push_back( xRect.getWKB( ) );
            vuiIds.push_back( xTestTable.insert( vxRectanglesWKB.back( ) ) );
        }
    }

    {
        auto pDatabase = std::make_shared<SQLDB<MySQLConDB>>( );
        SQLQuery<SQLDB<MySQLConDB>, WKBRectangle> xQuery(
            pDatabase, "SELECT ST_AsBinary(rectangle) FROM rectangle_test WHERE id = ?" );

        for( size_t uiI = 0; uiI < vxRectangles.size( ); uiI++ )
        {
            auto xWKB = xQuery.execAndGetNthCell<0>( vuiIds[ uiI ] );
            geomUtil::Rectangle<uint64_t> xRect;
            xRect.fromWKB( xWKB );
            std::cout << std::dec << xRect << std::endl;
            assert( xRect == vxRectangles[ uiI ] );
        }
    }

    return EXIT_SUCCESS;
} /// main function