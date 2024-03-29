/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "sql_api.h"

namespace libMSV
{

template <typename DBCon>
using SvCallerRunTableType = SQLTableWithAutoPriKey<DBCon,
                                                    std::string, // name
                                                    std::string, // desc
                                                    int64_t, // timestamp
                                                    PriKeyDefaultType // sv_jump_run_id
                                                    >;

const json jSvCallerRunTableDef = {
    { TABLE_NAME, "sv_caller_run_table" },
    { TABLE_COLUMNS,
      { { { COLUMN_NAME, "_name_" } },
        { { COLUMN_NAME, "_desc_" } }, // The column name was originally "desc", which is a keyword in MySQL
        { { COLUMN_NAME, "time_stamp" } },
        { { COLUMN_NAME, "sv_jump_run_id" } } } }
    // @todo ask arne about inserting NULL
    /*,{FOREIGN_KEY, {{COLUMN_NAME, "sv_jump_run_id"}, {REFERENCES, "sv_jump_run_table(id)"}}}*/ };

template <typename DBCon> class SvCallerRunTable : public SvCallerRunTableType<DBCon>
{

    std::shared_ptr<DBCon> pDatabase;
    SQLStatement<DBCon> xDelete; // Discuss Markus: Shouldn't this be a statement?
    SQLQuery<DBCon, int64_t> xGetId;
    SQLQuery<DBCon, int64_t> xHasId;
    SQLQuery<DBCon, int64_t> xGetIds;
    SQLQuery<DBCon, std::string, std::string, int64_t, int64_t> xGetName;
    SQLQuery<DBCon, uint64_t> xNum;
    SQLQuery<DBCon, uint64_t> xExists;
    SQLQuery<DBCon, uint64_t> xNameExists;
    SQLQuery<DBCon, int64_t> xNewestUnique;

  public:
    SvCallerRunTable( std::shared_ptr<DBCon> pDB )
        : SvCallerRunTableType<DBCon>( pDB, // the database where the table resides
                                       jSvCallerRunTableDef ),
          pDatabase( pDB ),
          xDelete( pDB, "DELETE FROM sv_caller_run_table WHERE _name_ = ?" ),
          xGetId( pDB, "SELECT id FROM sv_caller_run_table WHERE _name_ = ? ORDER BY time_stamp DESC LIMIT 1" ),
          xHasId( pDB, "SELECT COUNT(*) FROM sv_caller_run_table WHERE _name_ = ?" ),
          xGetIds( pDB, "SELECT id FROM sv_caller_run_table" ),
          xGetName( pDB, "SELECT _name_, _desc_, time_stamp, "
                         "CASE WHEN sv_jump_run_id is NULL THEN -1 ELSE sv_jump_run_id END AS v1 "
                         "FROM sv_caller_run_table WHERE id = ?" ),
          xNum( pDB, "SELECT COUNT(*) FROM sv_caller_run_table " ),
          xExists( pDB, "SELECT COUNT(*) FROM sv_caller_run_table WHERE id = ?" ),
          xNameExists( pDB, "SELECT COUNT(*) FROM sv_caller_run_table WHERE _name_ = ?" ),
          xNewestUnique(
              pDB,
              "SELECT id FROM sv_caller_run_table AS _outer_ WHERE ( SELECT COUNT(*) FROM sv_caller_run_table AS "
              "_inner_ WHERE _inner_._name_ = _outer_._name_ AND _inner_.time_stamp >= _outer_.time_stamp ) < ? "
              "AND _desc_ = ? " )
    {} // default constructor

    inline void deleteName( std::string& rS )
    {
        xDelete.exec( rS );
        // vDump( std::cout );
    } // method

    inline int64_t getId( std::string& rS )
    {
        // int64_t xId = xGetId.scalar( rS ); // FIXED: unused
        return xGetId.scalar( rS );
    } // method

    inline bool hasName( std::string& rS )
    {
        return xHasId.scalar( rS ) > 0;
    } // method

    inline std::vector<int64_t> getIds( )
    {
        return xGetIds.template executeAndStoreInVector<0>( );
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
        return xGetName.template execAndGetNthCell<0>( iId );
    } // method

    inline std::string getDesc( int64_t iId )
    {
        return xGetName.template execAndGetNthCell<1>( iId );
    } // method

    inline int64_t getSvJumpRunId( int64_t iId )
    {
        return xGetName.template execAndGetNthCell<3>( iId );
    } // method

    inline std::string getDate( int64_t iId )
    {
        auto now_c = ( std::time_t )( xGetName.template execAndGetNthCell<2>( iId ) );
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

    using ColTypesForw = TypePack<std::string, std::string, int64_t>; // redefine this to match the insert function

    inline PriKeyDefaultType insert( std::string sName, std::string sDesc, int64_t uiJumpRunId )
    {
        auto iTimeNow = (int64_t)std::chrono::system_clock::to_time_t( std::chrono::system_clock::now( ) );
        if( uiJumpRunId < 0 )
        {
            // Insert NULL at the position of uiJumpRunId
            return SvCallerRunTableType<DBCon>::insertNonSafe( sName, sDesc, iTimeNow, nullptr );
        }
        return SvCallerRunTableType<DBCon>::insert( sName, sDesc, iTimeNow, uiJumpRunId );
    } // method

    inline std::vector<int64_t> getNewestUnique( uint32_t uiNum, std::string sDesc )
    {
        return xNewestUnique.template executeAndStoreInVector<0>( uiNum, sDesc );
    } // method
}; // class

} // namespace libMSV
