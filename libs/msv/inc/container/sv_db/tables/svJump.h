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
using SvJumpTableType = SQLTableWithLibIncrPriKey<DBCon,
                                                  PriKeyDefaultType, // sv_jump_run_id (foreign key)
                                                  PriKeyDefaultType, // read_id (foreign key)
                                                  int64_t, // sort_pos_start
                                                  int64_t, // sort_pos_end
                                                  uint32_t, // from_pos
                                                  uint32_t, // to_pos
                                                  uint32_t, // query_from
                                                  uint32_t, // query_to
                                                  uint32_t, // num_supporting_nt
                                                  bool, // from_forward
                                                  bool, // to_forward
                                                  bool // from_seed_start
                                                  >;


template <typename DBCon> class SvJumpTable : public SvJumpTableType<DBCon>
{
    std::shared_ptr<DBCon> pDatabase;
    SQLQuery<DBCon, uint32_t> xQuerySize;
    SQLStatement<DBCon> xDeleteRun;

  public:
    json jSvCallTableDef( )
    {
        return json{
            {TABLE_NAME, "sv_jump_table"},
            {TABLE_COLUMNS,
             {{{COLUMN_NAME, "sv_jump_run_id"}},
              {{COLUMN_NAME, "read_id"}},
              {{COLUMN_NAME, "sort_pos_start"}},
              {{COLUMN_NAME, "sort_pos_end"}},
              {{COLUMN_NAME, "from_pos"}},
              {{COLUMN_NAME, "to_pos"}},
              {{COLUMN_NAME, "query_from"}},
              {{COLUMN_NAME, "query_to"}},
              {{COLUMN_NAME, "num_supporting_nt"}},
              {{COLUMN_NAME, "from_forward"}},
              {{COLUMN_NAME, "to_forward"}},
              {{COLUMN_NAME, "from_seed_start"}}}},
            {FOREIGN_KEY, {{COLUMN_NAME, "sv_jump_run_id"}, {REFERENCES, "sv_jump_run_table(id) ON DELETE CASCADE"}}},
            {FOREIGN_KEY, {{COLUMN_NAME, "read_id"}, {REFERENCES, "read_table(id)"}}}};
    } // method
#if DEBUG_LEVEL == 0
    // increase the bulk inserter size on this table
    using uiBulkInsertSize = std::integral_constant<size_t, 5000>;
#endif

    SvJumpTable( std::shared_ptr<DBCon> pDatabase )
        : SvJumpTableType<DBCon>( pDatabase, // the database where the table resides
                                  jSvCallTableDef( ) ),
          pDatabase( pDatabase ),
          xQuerySize( pDatabase, "SELECT COUNT(*) FROM sv_jump_table WHERE sv_jump_run_id = ?" ),
          xDeleteRun( pDatabase, "DELETE FROM sv_jump_table WHERE sv_jump_run_id IN ( SELECT id FROM "
                                 "sv_jump_run_table WHERE name = ?)" )
    {} // default constructor


    // @todo make the jumps rectangles as well and then use an r-tree for the sweep?
    inline void createIndices( int64_t uiRun )
    {
        // https://www.sqlite.org/queryplanner.html -> 3.2. Searching And Sorting With A Covering Index

        // index intended for the sweep over the start of all sv-rectangles
        // interestingly sv_jump_run_id needs to be part of the index even if it's in the condition...
        this->addIndex(
            json{{INDEX_NAME, "sort_start"},
                 {INDEX_COLUMNS, "sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward,"
                                 " to_forward, from_seed_start, num_supporting_nt, id, read_id, sv_jump_run_id"}} );

        // index intended for the sweep over the end of all sv-rectangles
        this->addIndex(
            json{{INDEX_NAME, "sort_end"},
                 {INDEX_COLUMNS, "sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward,"
                                 " to_forward, from_seed_start, num_supporting_nt, id, read_id, sv_jump_run_id"}} );
    } // method

    inline void dropIndices( int64_t uiRun )
    {
        this->dropIndex( json{{INDEX_NAME, "sort_start"}} );

        // index intended for the sweep over the end of all sv-rectangles
        this->dropIndex( json{{INDEX_NAME, "sort_end"}} );
    }

    inline uint32_t numJumps( int64_t jump_run_id )
    {
        return xQuerySize.scalar( jump_run_id );
    } // method

    inline void deleteRun( std::string& rS )
    {
        xDeleteRun.execAndBind( rS );
    } // method
}; // class


} // namespace libMSV
