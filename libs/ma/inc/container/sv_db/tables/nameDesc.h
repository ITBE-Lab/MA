/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "db_base.h"
#include "sql_api.h"

namespace libMA
{
template <typename DBCon>
using NameDescTableType = SQLTableWithAutoPriKey<DBCon,
                                                 std::string, // name
                                                 std::string, // desc
                                                 int64_t // timestamp
                                                 >;

/**
 * @brief template for descriptive table
 * @details
 * template for table that saves name, description, timestamp and id
 */
template <typename DBCon> class NameDescTable : public NameDescTableType<DBCon>
{

    std::shared_ptr<DBCon> pDatabase;
    SQLStatement<DBCon> xDelete;
    SQLQuery<DBCon, int64_t> xGetId;
    SQLQuery<DBCon, std::string, std::string, int64_t> xGetName;
    SQLQuery<DBCon, uint32_t> xNum;
    SQLQuery<DBCon, uint32_t> xExists;
    SQLQuery<DBCon, uint32_t> xNameExists;
    SQLQuery<DBCon, int64_t> xNewestUnique;

  public:

    NameDescTable( std::shared_ptr<DBCon> pDatabase, std::string sTableName )
        : NameDescTableType<DBCon>(
              pDatabase, // the database where the table resides
              json{
                  {TABLE_NAME, sTableName},
                  {TABLE_COLUMNS,
                   {{{COLUMN_NAME, "name"}},
                    {{COLUMN_NAME, "_desc_"}}, // The column name was originally "desc", which is a keyword in MySQL
                    {{COLUMN_NAME, "time_stamp"}}}},
              } ),
          pDatabase( pDatabase ),
          xDelete( pDatabase, ( "DELETE FROM " + sTableName + " WHERE name = ?" ).c_str( ) ),
          xGetId( pDatabase,
                  ( "SELECT id FROM " + sTableName + " WHERE name = ? ORDER BY time_stamp ASC LIMIT 1" ).c_str( ) ),
          xGetName( pDatabase, ( "SELECT name, _desc_, time_stamp FROM " + sTableName + " WHERE id = ?" ).c_str( ) ),
          xNum( pDatabase, ( "SELECT COUNT(*) FROM " + sTableName ).c_str( ) ),
          xExists( pDatabase, ( "SELECT COUNT(*) FROM " + sTableName + " WHERE id = ?" ).c_str( ) ),
          xNameExists( pDatabase, ( "SELECT COUNT(*) FROM " + sTableName + " WHERE name = ?" ).c_str( ) ),
          // FIXED: outer and inner are reserved words in MySQL
          xNewestUnique(
              pDatabase,
              ( "SELECT id FROM " + sTableName + " AS _outer_ WHERE ( SELECT COUNT(*) FROM " + sTableName +
                " AS _inner_ WHERE _inner_.name = _outer_.name AND _inner_.time_stamp >= _outer_.time_stamp ) < ?" )
                  .c_str( ) )
    {} // default constructor

    inline void deleteName( std::string& rS )
    {
        xDelete.execAndBind( rS );
        // vDump( std::cout );
    } // method

    inline int64_t getId( std::string& rS )
    {
        return xGetId.scalar( rS );
    } // method

    inline bool exists( int64_t iId )
    {
        return xExists.scalar( iId ) > 0;
    } // method

    inline bool nameExists( std::string sName )
    {
        return xNameExists.scalar( sName ) > 0;
    } // method

    inline std::string getName( int64_t iId )
    {
        // return std::get<0>(xGetName.vExecuteAndReturnIterator(iId).get());
        return xGetName.execAndGetNthCell<0>( iId );
    } // method

    inline std::string getDesc( int64_t iId )
    {
        // return std::get<1>( xGetName.vExecuteAndReturnIterator( iId ).get( ) );
        return xGetName.execAndGetNthCell<1>( iId );
    } // method

    inline std::string getDate( int64_t iId )
    {
        // auto now_c = (std::time_t)std::get<2>( xGetName.vExecuteAndReturnIterator( iId ).get( ) );
        auto now_c = ( std::time_t )( xGetName.execAndGetNthCell<2>( iId ) );
        std::stringstream ss;
#ifdef _MSC_VER
#pragma warning( suppress : 4996 ) // @todo find another way to do this
        ss << std::put_time( std::localtime( &now_c ), "%c" );
#else
        ss << std::put_time( std::localtime( &now_c ), "%c" );
#endif
        return ss.str( );
    } // method

    inline uint32_t size( )
    {
        return xNum.scalar( );
    } // method

    using ColTypesForw = TypePack<std::string, std::string>; // redefine this to match the insert function
    inline int64_t insert( std::string sName, std::string sDesc )
    {
        // return this->xInsertRow(
        return NameDescTableType<DBCon>::insert(
            sName, sDesc, (int64_t)std::chrono::system_clock::to_time_t( std::chrono::system_clock::now( ) ) );
    } // method

    inline std::vector<int64_t> getNewestUnique( uint32_t uiNum )
    {
        // return xNewestUnique.executeAndStoreInVector<0>(uiNum);
        throw std::runtime_error( "executeAndStoreInVector not implemented yet." );
        return std::vector<int64_t>( );
    } // method
}; // class

template <typename DBCon> class SvJumpRunTable : public NameDescTable<DBCon>
{
  public:
    SvJumpRunTable( std::shared_ptr<DBCon> pDatabase )
        : NameDescTable<DBCon>( pDatabase, "sv_jump_run_table" )
    {} // default constructor
}; // class

} // namespace libMA
