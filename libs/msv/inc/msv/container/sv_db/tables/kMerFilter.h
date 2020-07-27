/**
 * @file kMerFilter.h
 * @details
 * Stores a set of kMers that act as filter
 */
#pragma once

#include "ma/container/nucSeq.h"
#include "sql_api.h" // NEW DATABASE INTERFACE

namespace libMSV
{


template <typename DBCon>
using KMerFilterTableType = SQLTable<DBCon,

                                     PriKeyDefaultType, // sequencer id (foreign key)
                                     std::shared_ptr<libMA::CompressedNucSeq>, // k_mer
                                     uint32_t // num_occ
                                     >;
const json jKMerFilterDef = {
    {TABLE_NAME, "k_mer_filter_table"},
    {TABLE_COLUMNS, {{{COLUMN_NAME, "sequencer_id"}}, {{COLUMN_NAME, "k_mer"}}, {{COLUMN_NAME, "num_occ"}}}},
    {CONSTRAINTS, "UNIQUE(sequencer_id, k_mer)"},
    {FOREIGN_KEY, {{COLUMN_NAME, "sequencer_id"}, {REFERENCES, "sequencer_table(id)"}}}};
/**
 * @brief this table saves k-mers
 */
template <typename DBCon> class KMerFilterTable : public KMerFilterTableType<DBCon>
{
  public:
    KMerFilterTable( std::shared_ptr<DBCon> pDB ) : KMerFilterTableType<DBCon>( pDB, jKMerFilterDef )
    {} // default constructor

    // override
    virtual std::string makeInsertStmt( const size_t uiNumVals = 1 ) const
    {
        auto sStmt = KMerFilterTableType<DBCon>::makeInsertStmt( uiNumVals );
        sStmt.append( " ON DUPLICATE KEY UPDATE num_occ = num_occ + 1 " );
        return sStmt;
    } // function
}; // class

} // namespace libMSV
