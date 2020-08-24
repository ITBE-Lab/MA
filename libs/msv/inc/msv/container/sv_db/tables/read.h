/**
 * @file read.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "msv/container/svJump.h"
#include "sql_api.h" // NEW DATABASE INTERFACE

namespace libMSV
{


template <typename DBCon>
using ReadTableType = SQLTableWithLibIncrPriKey<DBCon,
                                                PriKeyDefaultType, // sequencer id (foreign key)
                                                std::string, // read name
                                                std::shared_ptr<CompressedNucSeq> // read sequence
                                                >;
const json jReadTableDef = {
    {TABLE_NAME, "read_table"},
    {TABLE_COLUMNS, {{{COLUMN_NAME, "sequencer_id"}}, {{COLUMN_NAME, "_name_"}}, {{COLUMN_NAME, "sequence"}}}},
    {FOREIGN_KEY, {{COLUMN_NAME, "sequencer_id"}, {REFERENCES, "sequencer_table(id)"}}}};
/**
 * @brief this table saves reads
 */
template <typename DBCon> class ReadTable : public ReadTableType<DBCon>
{

  public:
    // std::shared_ptr<DBCon> pDB;
    SQLQuery<DBCon, PriKeyDefaultType> xGetReadId;
    SQLQuery<DBCon, uint64_t> xGetNumMatchingReads;
    SQLQuery<DBCon, std::shared_ptr<CompressedNucSeq>, std::string> xGetRead;
    SQLQuery<DBCon, PriKeyDefaultType> xGetSeqId;

    ReadTable( std::shared_ptr<DBCon> pDB )
        : ReadTableType<DBCon>( pDB, jReadTableDef ),

          // pDB( pDB ),
          xGetReadId( pDB, "SELECT id FROM read_table WHERE sequencer_id = ? AND _name_ = ? " ),
          xGetNumMatchingReads( pDB, "SELECT COUNT(*) FROM read_table WHERE sequencer_id = ? AND _name_ = ? " ),
          xGetRead( pDB, "SELECT sequence, _name_ FROM read_table WHERE id = ? " ),
          xGetSeqId( pDB, "SELECT sequencer_id FROM read_table WHERE id = ? " )
    {} // default constructor


    inline int64_t insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pRead )
    {
        return ReadTableType<DBCon>::insert( uiSequencerId, pRead->sName, makeSharedCompNucSeq( *pRead ) );
    } // method

    inline std::shared_ptr<NucSeq> getRead( int64_t iId )
    {
        if( !xGetRead.execAndFetch( (PriKeyDefaultType)iId ) )
            throw std::runtime_error( "Read with id " + std::to_string( iId ) +
                                      " could not be found in the database." );
        auto xTuple = xGetRead.get( );
        std::get<0>( xTuple )->pUncomNucSeq->iId = iId;
        std::get<0>( xTuple )->pUncomNucSeq->sName = std::get<1>( xTuple );
        if( xGetRead.next( ) )
            assert( false ); // can never find two reads with same ID
        return std::get<0>( xTuple )->pUncomNucSeq;
    } // method
    inline int64_t getSeqId( int64_t iReadId )
    {
        return xGetSeqId.scalar( iReadId );
    } // method

    inline int64_t getReadId( int64_t iSeqId, std::string sName )
    {
        // char* pcEscapedName = PQescapeLiteral( this->pDB->pPGConn, sName.c_str( ), sName.size( ) );
        // std::string sEscapedName( pcEscapedName );
        if( xGetNumMatchingReads.scalar( (PriKeyDefaultType)iSeqId, sName ) == 0 )
            return -1;
        return xGetReadId.scalar( (PriKeyDefaultType)iSeqId, sName );
    } // method

    inline std::vector<std::shared_ptr<NucSeq>> getUsedReads( std::shared_ptr<DBCon> pDB )
    {
        SQLQuery<DBCon, std::shared_ptr<CompressedNucSeq>, std::string> xGetAllUsedReads(
            pDB, "SELECT sequence, _name_ FROM read_table WHERE id IN (SELECT DISTINCT read_id FROM sv_jump_table)" );
        std::vector<std::shared_ptr<NucSeq>> vRet;
        xGetAllUsedReads.execAndForAll( [&]( std::shared_ptr<CompressedNucSeq> pComp, std::string sName ) {
            vRet.push_back( pComp->pUncomNucSeq );
            vRet.back( )->sName = sName;
        } );
        return vRet;
    } // method
}; // class

} // namespace libMSV
