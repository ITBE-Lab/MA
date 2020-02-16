/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "common.h" // NEW DATABASE INTERFACE
#include "container/svJump.h"
#include "db_sql.h"

namespace libMA
{


template <typename DBCon>
using ReadTableType = SQLTableWithLibIncrPriKey<DBCon,
                                             PriKeyDefaultType, // sequencer id (foreign key)
                                             std::string, // read name
                                             std::shared_ptr<CompressedNucSeq> // read sequence
                                             >;
const json jReadTableDef = {
    {TABLE_NAME, "read_table"},
    {TABLE_COLUMNS, {{{COLUMN_NAME, "sequencer_id"}}, {{COLUMN_NAME, "name"}}, {{COLUMN_NAME, "sequence"}}}},
    {FOREIGN_KEY, {{COLUMN_NAME, "sequencer_id"}, {REFERENCES, "sequencer_table(id)"}}}};
/**
 * @brief this table saves reads
 */
template <typename DBCon> class ReadTable : public ReadTableType<DBCon>
{
    bool bDoDuplicateWarning = true;

  public:
    SQLQuery<DBCon, int32_t> xGetReadId;
    SQLQuery<DBCon, std::shared_ptr<CompressedNucSeq>, std::string> xGetRead;

    ReadTable( std::shared_ptr<DBCon> pDB )
        : ReadTableType<DBCon>( pDB, jReadTableDef ),
          xGetReadId( pDB, "SELECT id FROM read_table WHERE sequencer_id = ? AND name = ? " ),
          xGetRead( pDB, "SELECT sequence, name FROM read_table WHERE id = ? " )
    {} // default constructor


    inline int64_t insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pRead )
    {
        return ReadTableType<DBCon>::insert( uiSequencerId, pRead->sName, makeSharedCompNucSeq( *pRead ) );
    } // method

    inline std::shared_ptr<NucSeq> getRead( int64_t iId )
    {
        if( !xGetRead.execAndFetch( iId ) )
            throw std::runtime_error( "Read with id " + std::to_string( iId ) +
                                      " could not be found in the database." );
        auto xTuple = xGetRead.get( );
        std::get<0>( xTuple )->pUncomNucSeq->iId = iId;
        std::get<0>( xTuple )->pUncomNucSeq->sName = std::get<1>( xTuple );
        assert( !xGetRead.next( ) );
        return std::get<0>( xTuple )->pUncomNucSeq;
    } // method
}; // class

} // namespace libMA
