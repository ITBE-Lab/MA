/**
 * @file kMerFilter.h
 * @details
 * Stores a set of kMers that act as filter
 */
#pragma once

#include "ma/container/nucSeq.h"
#include "msv/module/count_k_mers.h"
#include "sql_api.h" // NEW DATABASE INTERFACE

namespace libMSV
{


template <typename DBCon>
using KMerFilterTableType = SQLTable<DBCon,
                                     PriKeyDefaultType, // sequencer id (foreign key)
                                     NucSeqSql, // k_mer
                                     uint32_t // num_occ
                                     >;
const json jKMerFilterDef = {
    {TABLE_NAME, "k_mer_filter_table"},
    {TABLE_COLUMNS, {{{COLUMN_NAME, "sequencer_id"}}, {{COLUMN_NAME, "k_mer"}}, {{COLUMN_NAME, "num_occ"}}}},
    {FOREIGN_KEY, {{COLUMN_NAME, "sequencer_id"}, {REFERENCES, "sequencer_table(id)"}}}};
/**
 * @brief this table saves k-mers
 */
template <typename DBCon> class KMerFilterTable : public KMerFilterTableType<DBCon>
{
  public:
    SQLQuery<DBCon, NucSeqSql, uint32_t> xGetAll;
    KMerFilterTable( std::shared_ptr<DBCon> pDB )
        : KMerFilterTableType<DBCon>( pDB, jKMerFilterDef ),
          xGetAll( pDB, "SELECT k_mer, num_occ FROM k_mer_filter_table WHERE sequencer_id = ? " )
    {} // default constructor

    void insert_counter_set( PriKeyDefaultType uiSequencerId, std::shared_ptr<KMerCounter> pCounter, size_t uiT )
    {
        auto pBulkInserter = this->template getBulkInserter<500>( );
        pCounter->iterate( [&]( const NucSeq& xNucSeq, size_t uiCnt ) {
            if( uiCnt > uiT )
            {
                auto pInsert = std::make_shared<NucSeq>( xNucSeq );
                pBulkInserter->insert( uiSequencerId, NucSeqSql( pInsert ), uiCnt );
            } // if
        } );
    } // method

    std::shared_ptr<KMerCounter> getCounter( PriKeyDefaultType uiSeqId, size_t uiW )
    {
        auto pCnt = std::make_shared<KMerCounter>( uiW );
        xGetAll.execAndForAll(
            [&]( NucSeqSql xNucSeq, uint32_t uiNumOcc ) { pCnt->addKMer( *xNucSeq.pNucSeq, uiNumOcc ); }, uiSeqId );
        return pCnt;
    }
}; // class

} // namespace libMSV
