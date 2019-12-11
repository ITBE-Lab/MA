/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "db_sql.h"

namespace libMA
{

typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string, // name
                                                 std::string, // desc
                                                 int64_t // timestamp
                                                 >
    TP_NAME_DESC_TABLE;
/**
 * @brief template for descriptive table
 * @details
 * template for table that saves name, description, timestamp and id
 */
class NameDescTable : public TP_NAME_DESC_TABLE
{
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    const std::string sTableName;
    CppSQLiteExtQueryStatement<int64_t> xDelete;
    CppSQLiteExtQueryStatement<int64_t> xGetId;
    CppSQLiteExtQueryStatement<std::string, std::string, int64_t> xGetName;
    CppSQLiteExtQueryStatement<uint32_t> xNum;
    CppSQLiteExtQueryStatement<uint32_t> xExists;
    CppSQLiteExtQueryStatement<uint32_t> xNameExists;
    CppSQLiteExtQueryStatement<int64_t> xNewestUnique;

  public:
    NameDescTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase, const std::string sTableName )
        : TP_NAME_DESC_TABLE( *pDatabase, // the database where the table resides
                              sTableName, // name of the table in the database
                              // column definitions of the table
                              std::vector<std::string>{"name", "desc", "time_stamp"} ),
          pDatabase( pDatabase ),
          sTableName( sTableName ),
          xDelete( *pDatabase, ( "DELETE FROM " + sTableName + " WHERE name == ?" ).c_str( ) ),
          xGetId( *pDatabase,
                  ( "SELECT id FROM " + sTableName + " WHERE name == ? ORDER BY time_stamp ASC LIMIT 1" ).c_str( ) ),
          xGetName( *pDatabase, ( "SELECT name, desc, time_stamp FROM " + sTableName + " WHERE id == ?" ).c_str( ) ),
          xNum( *pDatabase, ( "SELECT COUNT(*) FROM " + sTableName ).c_str( ) ),
          xExists( *pDatabase, ( "SELECT COUNT(*) FROM " + sTableName + " WHERE id == ?" ).c_str( ) ),
          xNameExists( *pDatabase, ( "SELECT COUNT(*) FROM " + sTableName + " WHERE name == ?" ).c_str( ) ),
          xNewestUnique( *pDatabase,
                         ( "SELECT id FROM " + sTableName + " AS outer WHERE ( SELECT COUNT(*) FROM " + sTableName +
                           " AS inner WHERE inner.name = outer.name AND inner.time_stamp >= outer.time_stamp ) < ?" )
                             .c_str( ) )
    {} // default constructor

    inline void deleteName( std::string& rS )
    {
        xDelete.bindAndExecQuery<>( rS );
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
        return std::get<0>( xGetName.vExecuteAndReturnIterator( iId ).get( ) );
    } // method

    inline std::string getDesc( int64_t iId )
    {
        return std::get<1>( xGetName.vExecuteAndReturnIterator( iId ).get( ) );
    } // method

    inline std::string getDate( int64_t iId )
    {
        auto now_c = (std::time_t)std::get<2>( xGetName.vExecuteAndReturnIterator( iId ).get( ) );
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

    inline int64_t insert( std::string sName, std::string sDesc )
    {
        return this->xInsertRow(
            sName, sDesc, (int64_t)std::chrono::system_clock::to_time_t( std::chrono::system_clock::now( ) ) );
    } // method

    inline std::vector<int64_t> getNewestUnique( uint32_t uiNum )
    {
        return xNewestUnique.executeAndStoreInVector<0>( uiNum );
    } // method
}; // class

} // namespace libMA