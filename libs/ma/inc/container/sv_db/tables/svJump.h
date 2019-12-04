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

typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int64_t, // sv_jump_run_id (foreign key)
                                                 int64_t, // read_id (foreign key)
                                                 int64_t, // sort_pos_start
                                                 int64_t, // sort_pos_end
                                                 uint32_t, // from_pos
                                                 uint32_t, // to_pos
                                                 uint32_t, // query_from
                                                 uint32_t, // query_to
                                                 uint32_t, // num_supporting_nt
                                                 bool, // from_forward @todo save space by compressing booleans?
                                                 bool, // to_forward
                                                 bool // from_seed_start
                                                 >
    TP_SV_JUMP_TABLE;
class SvJumpTable : public TP_SV_JUMP_TABLE
{
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    CppSQLiteExtQueryStatement<uint32_t> xQuerySize;
    CppSQLiteExtQueryStatement<int64_t> xDeleteRun;

  public:
    SvJumpTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
        : TP_SV_JUMP_TABLE( *pDatabase, // the database where the table resides
                            "sv_jump_table", // name of the table in the database
                            // column definitions of the table
                            std::vector<std::string>{"sv_jump_run_id", "read_id", "sort_pos_start", "sort_pos_end",
                                                     "from_pos", "to_pos", "query_from", "query_to",
                                                     "num_supporting_nt", "from_forward", "to_forward",
                                                     "from_seed_start"},
                            // constraints for table
                            std::vector<std::string>{
                                "FOREIGN KEY (sv_jump_run_id) REFERENCES sv_jump_run_table(id) ON DELETE CASCADE",
                                "FOREIGN KEY (read_id) REFERENCES read_table(id)"} ),
          pDatabase( pDatabase ),
          xQuerySize( *pDatabase, "SELECT COUNT(*) FROM sv_jump_table" ),
          xDeleteRun( *pDatabase, "DELETE FROM sv_jump_table WHERE sv_jump_run_id IN ( SELECT id FROM "
                                  "sv_jump_run_table WHERE name == ?)" )
    {} // default constructor

    inline void createIndices( int64_t uiRun )
    {
        // https://www.sqlite.org/queryplanner.html -> 3.2. Searching And Sorting With A Covering Index
        // index intended for the sweep over the start of all sv-rectangles
        // interestingly sv_jump_run_id needs to be part of the index even if it's in the condition...
        pDatabase->execDML( ( "CREATE INDEX IF NOT EXISTS sv_jump_table_sort_index_start_" + std::to_string( uiRun ) +
                              " ON sv_jump_table"
                              "(sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward,"
                              " to_forward, from_seed_start, num_supporting_nt, id, read_id, sv_jump_run_id) "
                              "WHERE sv_jump_run_id == " +
                              std::to_string( uiRun ) )
                                .c_str( ) );
        // index intended for the sweep over the end of all sv-rectangles
        pDatabase->execDML( ( "CREATE INDEX IF NOT EXISTS sv_jump_table_sort_index_end_" + std::to_string( uiRun ) +
                              " ON sv_jump_table"
                              "(sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward,"
                              " to_forward, from_seed_start, num_supporting_nt, id, read_id, sv_jump_run_id) "
                              "WHERE sv_jump_run_id == " +
                              std::to_string( uiRun ) )
                                .c_str( ) );
    } // method

    inline uint32_t numJumps( )
    {
        return xQuerySize.scalar( );
    } // method

    inline void deleteRun( std::string& rS )
    {
        xDeleteRun.bindAndExecQuery<>( rS );
    } // method
}; // class

} // namespace libMA