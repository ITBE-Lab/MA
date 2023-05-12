/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "sql_api.h" // NEW DATABASE INTERFACE

namespace libMSV
{


class CoverageCollector : public libMS::Container
{
public:
    const size_t uiBinSize;
    std::vector<std::atomic<size_t>> vCoverage;

    CoverageCollector(libMA::nucSeqIndex uiGenomeSize) :
     uiBinSize(pGlobalParams->xCoverageBinSize->get()),
     vCoverage(uiGenomeSize / uiBinSize)
    {}

    inline void addRead(std::set<int64_t> xCoverage)
    {
        for(size_t uiPos : xCoverage)
            vCoverage[uiPos]++; // atomic so no parralelization needed
    }

};


template <typename DBCon>
using CoverageTableType = SQLTableWithAutoPriKey<DBCon,
                                                   int64_t, // run_id
                                                   int64_t, // bin_id
                                                   int64_t // bin_coverage
                                                   >;
const json jCoverageTableDef = { { TABLE_NAME, "coverage_table" },
                                  { TABLE_COLUMNS, { { { COLUMN_NAME, "run_id" } }, 
                                                     { { COLUMN_NAME, "bin_id" } }, 
                                                     { { COLUMN_NAME, "bin_coverage" } } } } };

/**
 * @brief contains the name of a sequencer run
 * @details
 * In order to create sv jumps with different parameters for different read types,
 * each read type has to have a entry in the sequencer_table.
 * Then jumps can be computed separated for each sequencer_table entry.
 * Although the jumps are created separated, they can all be used together in the line sweep.
 */
template <typename DBCon> class CoverageTable : public CoverageTableType<DBCon>
{
    SQLQuery<DBCon, uint64_t> xQueryCoverage;
    SQLStatement<DBCon> xUpdateCoverage;


  public:
    CoverageTable( std::shared_ptr<DBCon> pDB ) : 
          CoverageTableType<DBCon>( pDB, jCoverageTableDef ),
          xQueryCoverage( pDB, "SELECT bin_coverage FROM coverage_table WHERE run_id = ? AND bin_id = ? LIMIT 1" ),
          xUpdateCoverage( pDB,
                     "UPDATE coverage_table "
                     "SET bin_coverage = bin_coverage + ? "
                     "WHERE run_id = ? " 
                     "AND   bin_id = ? " )
    {} // default constructor
    
    inline void createIndices( )
    {
        this->addIndex(
            json{ { INDEX_NAME, "bin_id_index" },
                  { INDEX_COLUMNS, "run_id, bin_id" } } );
    } // method

    inline void dropIndices( )
    {
        this->dropIndex( json{ { INDEX_NAME, "bin_id_index" } } );
    }

    inline uint32_t getCoverage( int64_t jump_run_id, int64_t uiBinId )
    {
        return xQueryCoverage.scalar( jump_run_id, uiBinId );
    } // method
    
    inline void incCoverage( int64_t jump_run_id, int64_t uiBinId, uint32_t uiAmount )
    {
        xUpdateCoverage.exec( uiAmount, jump_run_id, uiBinId );
    } // method

    inline void initCoverage( int64_t jump_run_id, int64_t uiNumBins )
    {
        for(size_t uiI = 0; uiI < uiNumBins; uiI++)
            CoverageTableType<DBCon>::insert( jump_run_id, uiI, 0 );
    }

    inline void initCoverageFromColl( int64_t jump_run_id, CoverageCollector& rColl )
    {
        for(size_t uiI = 0; uiI < rColl.vCoverage.size(); uiI++)
            CoverageTableType<DBCon>::insert( jump_run_id, uiI, rColl.vCoverage[uiI] );
    }
}; // class

} // namespace libMSV
