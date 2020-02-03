/* Authors: Arne Kutzner and Markus Schmidt
 * Created: Jan. 2020
 * MIT License
 * @file spatial.h
 * @brief Support of spatial types in SQL databases.
 */

#pragma once

#include "MySQL_con.h" // NEW DATABASE INTERFACE
#include "geom.h"

typedef geomUtil::WKBPolygon<geomUtil::Rectangle<uint64_t>::uiSizeWKB> MyRectangle;

/* Integration of geomUtil::Rectangle as data-type in the MySQL interface.
 */
// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string MySQLConDB::TypeTranslator::getSQLTypeName<MyRectangle>( )
{
    return "BLOB";
} // specialized method

// Part 2: Input arguments: Set the start of the blob (void *), size of the blob and type of the blob.
template <> inline void MySQLConDB::StmtArg::set( const MyRectangle& rRectangle )
{
    this->uiLength = static_cast<unsigned long>( geomUtil::Rectangle<uint64_t>::uiSizeWKB );
    pMySQLBind->buffer_length = static_cast<unsigned long>( geomUtil::Rectangle<uint64_t>::uiSizeWKB );
    pMySQLBind->buffer_type = MYSQL_TYPE_BLOB; // this type must be equal to the type in Part 3.
    pMySQLBind->buffer = rRectangle.getData( );
} // specialized method

// Part 3: Code for supporting query output:
//         1. Via the third argument of the call of init, set the MySQL datatype for your cell type.
//         2. Using storeVarSizeCel, fetch the blob from the byte-buffer of the cell.
template <> struct /* MySQLConDB:: */ RowCell<MyRectangle> : public /* MySQLConDB::*/ RowCellBase<MyRectangle>
{
    inline void init( MYSQL_BIND* pMySQLBind, MyRectangle* pCellValue, size_t uiColNum )
    {
        RowCellBase<MyRectangle>::init( pMySQLBind, pCellValue, MYSQL_TYPE_BLOB, uiColNum );
    } // method

    // Fetch the blob from the buffer.
    inline void storeVarSizeCell( )
    {
        //assert( this->uiLength == geomUtil::Rectangle<uint64_t>::uiSizeWKB );
        pCellValue->setData( this->pVarLenBuf.get( ) );
    } // method
}; // specialized class

inline std::ostream& operator<<( std::ostream& xOS, const MyRectangle& xRect )
{
    xOS << "WKBPolygon: ";
    for( auto uiI : xRect.aData )
        xOS << std::hex << (int)uiI << " ";
    xOS << std::endl;
    return xOS;
}