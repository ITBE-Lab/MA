/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "common.h" // NEW DATABASE INTERFACE
#include "db_sql.h"

namespace libMA
{

typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int64_t, // sequencer id (foreign key)
                                                 std::string, // read name
                                                 NucSeqSql // read sequence
                                                 >
    TP_READ_TABLE;
/**
 * @brief this table saves reads
 */
class ReadTable : public TP_READ_TABLE
{
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    bool bDoDuplicateWarning = true;

  public:
    CppSQLiteExtQueryStatement<int32_t> xGetReadId;
    CppSQLiteExtQueryStatement<NucSeqSql, std::string> xGetRead;

    ReadTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
        : TP_READ_TABLE( *pDatabase, // the database where the table resides
                         "read_table", // name of the table in the database
                         // column definitions of the table
                         std::vector<std::string>{"sequencer_id", "name", "sequence"},
                         // constraints for table
                         std::vector<std::string>{"FOREIGN KEY (sequencer_id) REFERENCES sequencer_table(id) "} ),
          pDatabase( pDatabase ),
          xGetReadId( *pDatabase, "SELECT id FROM read_table WHERE sequencer_id == ? AND name == ? " ),
          xGetRead( *pDatabase, "SELECT sequence, name FROM read_table WHERE id == ? " )
    {} // default constructor

    inline int64_t insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pRead )
    {
        return xInsertRow( uiSequencerId, pRead->sName, NucSeqSql( pRead ) );
    } // method

    inline std::shared_ptr<NucSeq> getRead( int64_t iId )
    {
        auto xTuple = xGetRead.vExecuteAndReturnIterator( iId ).get( );
        std::get<0>( xTuple ).pNucSeq->iId = iId;
        std::get<0>( xTuple ).pNucSeq->sName = std::get<1>( xTuple );
        return std::get<0>( xTuple ).pNucSeq;
    } // method
}; // class

/* NEW DATABASE INTERFACE */

template <typename DBCon>
using ReadTableType = SQLTableWithAutoPriKey<DBCon,
                                             int64_t, // sequencer id (foreign key)
                                             std::string, // read name
                                             NucSeqSql // read sequence
                                             >;
const json jReadTableDef = {
    {TABLE_NAME, "read_table"},
    {TABLE_COLUMNS, {{{COLUMN_NAME, "sequencer_id"}}, {{COLUMN_NAME, "name"}}, {{COLUMN_NAME, "sequence"}}}},
    {FOREIGN_KEY, {{COLUMN_NAME, "sequencer_id"}, {REFERENCES, "sequencer_table(id)"}}}};
/**
 * @brief this table saves reads
 */
template <typename DBCon> class _ReadTable : public ReadTableType<DBCon>
{
    bool bDoDuplicateWarning = true;

  public:
    SQLQuery<DBCon, int32_t> xGetReadId;
    SQLQuery<DBCon, NucSeqSql, std::string> xGetRead;

    _ReadTable( std::shared_ptr<SQLDB<DBCon>> pDB )
        : ReadTableType<DBCon>( pDB, jReadTableDef ),
          xGetReadId( pDB, "SELECT id FROM read_table WHERE sequencer_id = ? AND name = ? " ),
          xGetRead( pDB, "SELECT sequence, name FROM read_table WHERE id = ? " )
    {} // default constructor


    inline int64_t insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pRead )
    {
        return ReadTableType<DBCon>::insert( uiSequencerId, pRead->sName, NucSeqSql( pRead ) );
    } // method

    inline std::shared_ptr<NucSeq> getRead( int64_t iId )
    {
        if( !xGetRead.execAndFetch( iId ) )
            throw std::runtime_error( "Read with id " + std::to_string( iId ) +
                                      " could not be found in the database." );
        auto xTuple = xGetRead.get( );
        std::get<0>( xTuple ).pNucSeq->iId = iId;
        std::get<0>( xTuple ).pNucSeq->sName = std::get<1>( xTuple );
        assert( !xGetRead.next( ) );
        return std::get<0>( xTuple ).pNucSeq;
    } // method
}; // class

} // namespace libMA
