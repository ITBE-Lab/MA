/* Authors: Arne Kutzner and Markus Schmidt
 * Created: Nov. 2019
 * This file is part of the ITBE-Lab code collection.
 * MIT License
 * @file mysql_spatial.h
 * @brief Spatial datatype support for the MySQL engine support.
 *        NOTICE: This file is intended to be included via mysql.h merely.
 */

#include "wkb_spatial.h" // definitions of WKB data-structures

/* Integration of geom::Rectangle as data-type in the MySQL interface. */

// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string MySQLConDB::TypeTranslator::getSQLTypeName<WKBUint64Rectangle>( )
{
    return "BLOB"; // WKB data are passed as BLOB
} // specialized method

// Part1b : Spatial types require an indication that the argument passed at a placeholder's
//          position has the format 'WKB'.
template <> inline std::string MySQLConDB::TypeTranslator::getPlaceholderForType<WKBUint64Rectangle>( )
{
    return "ST_PolyFromWKB(?, 0)";
} // specialized method

// Part 2: Input arguments: Set the start of the blob (void *), size of the blob and type of the blob.
template <> inline void MySQLConDB::StmtArg::set( const WKBUint64Rectangle& rRectangle )
{
    this->uiLength = static_cast<unsigned long>( WKBUint64Rectangle::uiSizeWKB );
    pMySQLBind->buffer_length = static_cast<unsigned long>( WKBUint64Rectangle::uiSizeWKB );
    pMySQLBind->buffer_type = MYSQL_TYPE_BLOB; // this type must be equal to the type in Part 3.
    pMySQLBind->buffer = rRectangle.getData( );

#if 0 // Debugging 
    std::cout << (int)( (uint8_t*)pMySQLBind->buffer )[ 0 ] << " " << (int)( (uint8_t*)pMySQLBind->buffer )[ 1 ] << " "
              << (int)( (uint8_t*)pMySQLBind->buffer )[ 2 ] << " " << (int)( (uint8_t*)pMySQLBind->buffer )[ 3 ] << " "
              << (int)( (uint8_t*)pMySQLBind->buffer )[ 4 ] << " len= " << pMySQLBind->buffer_length << std::endl;
#endif
} // specialized method

// Part 3: Code for supporting query output:
//         1. Via the third argument of the call of init, set the MySQL datatype for your cell type.
//         2. Using storeVarSizeCel, fetch the blob from the byte-buffer of the cell.
template <>
struct /* MySQLConDB:: */ RowCell<WKBUint64Rectangle> : public /* MySQLConDB::*/ RowCellBase<WKBUint64Rectangle>
{
    inline void init( MYSQL_BIND* pMySQLBind, WKBUint64Rectangle* pCellValue, size_t uiColNum )
    {
        RowCellBase<WKBUint64Rectangle>::init( pMySQLBind, pCellValue, MYSQL_TYPE_BLOB, uiColNum );
    } // method

    // Fetch the blob from the buffer.
    inline void storeVarSizeCell( )
    {
        // assert( this->uiLength == geom::Rectangle<uint64_t>::uiSizeWKB );
        // Fill the buffer with WKB values
        pCellValue->setData( this->pVarLenBuf.get( ) );
    } // method
}; // specialized class
  