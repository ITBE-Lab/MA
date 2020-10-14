/* Authors: Arne Kutzner and Markus Schmidt
 * Created: July. 2020
 * This file is part of the ITBE-Lab code collection.
 * MIT License
 * @file postgre_sql_spatial.h
 * @brief Spatial datatype support for the PostgreSQL engine support.
 */

// Activation of spatial extensions for database:
// create extension postgis;
// In Ubuntu for DB postgres and user postgres: psql --user postgres -c 'create extension postgis' postgres

#include "wkb_spatial.h" // definitions of WKB data-structures

/* Integration of geom::Rectangle as data-type in the MySQL interface. */

// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string PostgreSQLDBCon::TypeTranslator::getSQLTypeName<WKBUint64Rectangle>( )
{
    return "bytea"; // WKB data are passed as BLOB
} // specialized method

// Part1b : Spatial types require an indication that the argument passed at a placeholder's
//          position has the format 'WKB'.
//          https://postgis.net/docs/ST_GeomFromWKB.html
template <>
inline std::string
PostgreSQLDBCon::TypeTranslator::getPlaceholderForType<WKBUint64Rectangle>( const std::string& rsInsertedText )
{
    return "ST_GeomFromWKB(" + rsInsertedText + ", 0)"; // Corresponding MySQL: ST_PolyFromWKB
} // specialized method

// Part 2: Input arguments: Set the start of the blob (void *), size of the blob and type of the blob.
template <> inline void PostgreSQLDBCon::StmtArg::set( const WKBUint64Rectangle& rRectangle )
{
    rpParamValue = (char*)rRectangle.getData( );
    riParamLength = static_cast<int>( WKBUint64Rectangle::uiSizeWKB );
    riParamFormat = PG_BINARY_ARG;

#if 0 // Debugging 
    std::cout << (int)((uint8_t*)pMySQLBind->buffer)[0] << " " << (int)((uint8_t*)pMySQLBind->buffer)[1] << " "
        << (int)((uint8_t*)pMySQLBind->buffer)[2] << " " << (int)((uint8_t*)pMySQLBind->buffer)[3] << " "
        << (int)((uint8_t*)pMySQLBind->buffer)[4] << " len= " << pMySQLBind->buffer_length << std::endl;
#endif
} // specialized method

// Part 3: Code for supporting query output:
//         1. Via the third argument of the call of init, set the MySQL datatype for your cell type.
//         2. Using storeVarSizeCel, fetch the blob from the byte-buffer of the cell.
template <>
struct /* MySQLConDB:: */ PGRowCell<WKBUint64Rectangle> : public /* MySQLConDB::*/ PGRowCellBase<WKBUint64Rectangle>
{
    inline void init( WKBUint64Rectangle* pCellValue, size_t uiColNum )
    {
        PGRowCellBase<WKBUint64Rectangle>::init( pCellValue, uiColNum );
    } // method

    inline void store( const PGresult* pPGRes )
    {
        pCellValue->setData( this->getValPtr( pPGRes ) );
    } // method
}; // specialized class

// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string PostgreSQLDBCon::TypeTranslator::getSQLTypeName<WKBPoint>( )
{
    return "POINT"; // WKB data are passed as BLOB
} // specialized method

// Part1b : Spatial types require an indication that the argument passed at a placeholder's
//          position has the format 'WKB'.
template <>
inline std::string PostgreSQLDBCon::TypeTranslator::getPlaceholderForType<WKBPoint>( const std::string& rsInsertedText )
{
    return "ST_PointFromWKB(" + rsInsertedText + ", 0)";
} // specialized method

// Part 2: Input arguments: Set the start of the blob (void *), size of the blob and type of the blob.
template <> inline void PostgreSQLDBCon::StmtArg::set( const WKBPoint& rPoint )
{
    rpParamValue = (char*)rPoint.getData( );
    riParamLength = static_cast<int>( WKBPoint::uiSizeWKB );
    riParamFormat = PG_BINARY_ARG;
} // specialized method

// Part 3: Code for supporting query output:
//         1. Via the third argument of the call of init, set the MySQL datatype for your cell type.
//         2. Using storeVarSizeCel, fetch the blob from the byte-buffer of the cell.
template <> struct /* PostgreSQLDBCon:: */ PGRowCell<WKBPoint> : public /* PostgreSQLDBCon::*/ PGRowCellBase<WKBPoint>
{
    inline void init( WKBPoint* pCellValue, size_t uiColNum )
    {
        PGRowCellBase<WKBPoint>::init( pCellValue, uiColNum );
    } // method

    inline void store( const PGresult* pPGRes )
    {
        pCellValue->setData( this->getValPtr( pPGRes ) );
    } // method
}; // specialized class