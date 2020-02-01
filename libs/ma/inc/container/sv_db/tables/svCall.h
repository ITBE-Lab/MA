/**
 * @file svCall.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "common.h"
#include "system.h"
#include <csignal>
#include <string>

// using DBCon = MySQLConDB; // FIXME: Implement as template parameter

namespace libMA
{
template <typename DBCon>
using SvCallTableType = SQLTableWithAutoPriKey<DBCon, // DB connector type
                                               int64_t, // sv_caller_run_id (foreign key)
                                               uint32_t, // from_pos (geometry)
                                               uint32_t, // to_pos (geometry)
                                               uint32_t, // from_size (geometry)
                                               uint32_t, // to_size (geometry)
                                               bool, // switch_strand
                                               NucSeqSql, // inserted_sequence
                                               uint32_t, // supporting_reads
                                               uint32_t, // reference_ambiguity
                                               int64_t // regex_id
                                               >;
#ifdef WITHOUT_R_TREE
// Emergency definition of the R-tree table for structural variations.
template <typename DBCon>
using SvCallRTreeType = SQLTable<DBCon, // DB connector type
                                 uint32_t, // id (geometry)
                                 uint32_t, // run_id_a (geometry)
                                 uint32_t, // run_id_b (geometry)
                                 uint32_t, // minX (geometry)
                                 uint32_t, // maxX (geometry)
                                 uint32_t, // minY (geometry)
                                 uint32_t // maxY (geometry)
                                 >;

const json jSvCallRTreeDef = { { TABLE_NAME, "sv_call_r_tree" },
                               { TABLE_COLUMNS,
                                 { { { COLUMN_NAME, "id" } },
                                   { { COLUMN_NAME, "run_id_a" } },
                                   { { COLUMN_NAME, "run_id_b" } },
                                   { { COLUMN_NAME, "minX" } },
                                   { { COLUMN_NAME, "maxX" } },
                                   { { COLUMN_NAME, "minY" } },
                                   { { COLUMN_NAME, "maxY" } } } } };
#endif

template <typename DBCon> class SvCallTable : private SvCallTableType<DBCon>
{
    std::shared_ptr<DBCon> pDatabase;

#ifndef WITHOUT_R_TREE
    class RTreeIndex
    {
      public:
        RTreeIndex( std::shared_ptr<DBCon> pDatabase )
        {
            // create a R*tree index
            if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
            {
                /* We drop the table in the case that it exists already. */
                pDatabase->execDML( "DROP TABLE IF EXISTS sv_call_r_tree" );
                pDatabase->execDML( "CREATE VIRTUAL TABLE sv_call_r_tree USING rtree_i32( "
                                    "       id, " // key from sv_call_table
                                    "       run_id_a, run_id_b, " // run id -> has to be a dimension for it to work
                                    "       minX, maxX, " // start & end of from positions
                                    "       minY, maxY " // start & end of to positions
                                    "   )" );

            } // if
        } // constructor
    }; // class
#else
    class RTreeIndex : public SvCallRTreeType<DBCon>
    {
      public:
        /** @brief Construct the R-tree table in given database. */
        RTreeIndex( std::shared_ptr<DBCon> pDB ) : SvCallRTreeType<DBCon>( pDB, jSvCallRTreeDef )
        {} // constructor
    }; // class
#endif // WITHOUT_R_TREE

  public:
    /**
     * @brief SQL code for score computation of sv call
     */
    static std::string getSqlForCallScore( std::string sTableName = "" )
    {
        if( !sTableName.empty( ) )
            sTableName.append( "." );
        return " ( " + sTableName + "supporting_reads * 1.0 ) / " + sTableName + "reference_ambiguity ";
    } // method

  private:
    class OverlapCache
    {
        typedef SQLTable<DBCon, // DB connector type
                         double, // score_a
                         int64_t, // blur_min
                         int64_t // blur_max
                         >
            OverlapCacheTableType;
        /**
         * @brief this speeds up the overlap analysis
         * @details
         * For the overlap analysis we have to compute the number of overlapping calls multiple times with different
         * minimal scores.
         * For huge databases, this takes forever due to a lot of searches in the R*Tree. However, even if we analyze
         * only once, we do not need to search the R*Tree multiple times if we cache all overlapping calls (no matter
         * what score) This table implements this.
         * @note if calls are changed, deleted, or new ones are inserted using this results in undefined behavior.
         */
        class OverlapCacheTable : public OverlapCacheTableType
        {
          public:
            json jOverlapCacheTableDef( std::string sTableName )
            {
                return json{ { TABLE_NAME, sTableName },
                             { TABLE_COLUMNS,
                               { { { COLUMN_NAME, "score_a" } },
                                 { { COLUMN_NAME, "blur_min" } },
                                 { { COLUMN_NAME, "blur_max" } } } } };
            } // method

            std::string sTableName;
            SQLQuery<DBCon, int64_t, double, uint32_t, uint32_t, uint32_t, uint32_t, bool> xNumOverlaps;
            SQLQuery<DBCon, uint32_t> xNumOverlapsCache;
            // DEL: SQLStatement<DBCon> xCreateIndex;

            OverlapCacheTable( std::shared_ptr<DBCon> pDatabase, int64_t iCallerRunIdA, int64_t iCallerRunIdB,
                               bool bDoCreate = false )
                : OverlapCacheTableType( pDatabase, // the database where the table resides
                                                    // name of the table in the database
                                         jOverlapCacheTableDef( "overlap_cache_table_" +
                                                                std::to_string( iCallerRunIdA ) + "_" +
                                                                std::to_string( iCallerRunIdB ) ) ),
                  sTableName( "overlap_cache_table_" + std::to_string( iCallerRunIdA ) + "_" +
                              std::to_string( iCallerRunIdB ) ),
                  xNumOverlaps( pDatabase,
                                // each inner call can overlap an outer call at most once
                                "SELECT id, " + SvCallTable::getSqlForCallScore( ) +
                                    ", from_pos, from_size, "
                                    "       to_pos, to_size, switch_strand "
                                    "FROM sv_call_table "
                                    "WHERE sv_caller_run_id = ? "
                                    "AND " +
                                    SvCallTable::getSqlForCallScore( ) + " >= ? " ),
                  xNumOverlapsCache( pDatabase,
                                     "SELECT COUNT(*) "
                                     "FROM " +
                                         sTableName +
                                         " WHERE score_a >= ? "
                                         "AND blur_min <= ? "
                                         "AND blur_max >= ? " )
            // create an index so that we can query the results quickly
            // DEL: xCreateIndex( pDatabase,
            // DEL:               "CREATE INDEX IF NOT EXISTS " + sTableName + "_index ON " + sTableName + " (score_a) "
            // DEL: )
            {} // default constructor
        }; // class (OverlapCacheTable)

        typedef SQLTable<DBCon, // DB connector type
                         int64_t, // sv_caller_run_id_a
                         int64_t, // sv_caller_run_id_b
                         int64_t // blur
                         >
            OverlapCacheComputedTableType;

        class OverlapCacheComputedTable : public OverlapCacheComputedTableType
        {
          public:
            json jOverlapCacheComputedTableDef( )
            {
                return json{ { TABLE_NAME, "overlap_cache_computed_table" },
                             { TABLE_COLUMNS,
                               {
                                   { { COLUMN_NAME, "sv_caller_run_id_a" } },
                                   { { COLUMN_NAME, "sv_caller_run_id_b" } },
                                   { { COLUMN_NAME, "blur" } },
                               } } };
            } // method

            SQLQuery<DBCon, int64_t> xIsComputed;
            SQLStatement<DBCon> xDeleteFromTable;

            OverlapCacheComputedTable( std::shared_ptr<DBCon> pDatabase )
                : OverlapCacheComputedTableType( pDatabase, // the database where the table resides
                                                 jOverlapCacheComputedTableDef( ) ), // table definition
                  xIsComputed( pDatabase,
                               "SELECT COUNT(*) "
                               "FROM overlap_cache_computed_table "
                               "WHERE sv_caller_run_id_a = ? "
                               "AND sv_caller_run_id_b = ? "
                               "AND blur >= ? " ),
                  xDeleteFromTable( pDatabase,
                                    "DELETE FROM overlap_cache_computed_table "
                                    "WHERE sv_caller_run_id_a = ? "
                                    "AND sv_caller_run_id_b = ? " )
            {} // default constructor

        }; // class

        std::shared_ptr<DBCon> pDatabase;
        std::shared_ptr<std::mutex> pWriteLock;
        std::string sDBName;
        std::shared_ptr<OverlapCacheComputedTable> pHelper;

      public:
        OverlapCache( std::shared_ptr<DBCon> pDatabase, std::shared_ptr<std::mutex> pWriteLock,
                      std::string sDBName )
            : pDatabase( pDatabase ),
              pWriteLock( pWriteLock ),
              sDBName( sDBName ),
              pHelper( std::make_shared<OverlapCacheComputedTable>( pDatabase ) )
        {} // constructor

        // iCallerRunIdA used as ground truth
        void createCache( int64_t iCallerRunIdA, int64_t iCallerRunIdB, int64_t iAllowedDist )
        {
            if( pHelper->xIsComputed.scalar( iCallerRunIdA, iCallerRunIdB, iAllowedDist ) > 0 )
                return;

            // FIXME: Not implemented yet:
            // CppSQLiteExtImmediateTransactionContext xTransactionContext( *pDatabase );
            pHelper->xDeleteFromTable.exec( iCallerRunIdA, iCallerRunIdB );

            // create the cache table
            OverlapCacheTable xCacheTable( pDatabase, iCallerRunIdA, iCallerRunIdB, true );

            // simply get all calls in iCallerRunIdB
            auto vResults = xCacheTable.xNumOverlaps.executeAndStoreAllInVector( iCallerRunIdB, 0 );

#if 0 // Code replacement is in #else branch.
            {
                ThreadPool xPool( 32 ); // FIXME: Read this from parameters

                for( size_t uiJobId = 0; uiJobId < 32; uiJobId++ )
                    xPool.enqueue(
                        [ & ]( size_t, size_t uiJobId_ ) {
                            SQLDB<DBCon> xNewConn( sDBName );
                            SQLQuery<DBCon, int64_t> xNumOverlapsHelper1(
                                xNewConn,
                                // make sure that inner overlaps the outer:
                                "SELECT outer.id "
                                "FROM sv_call_table AS outer, sv_call_r_tree AS idx_outer "
                                "WHERE outer.id == idx_outer.id "
                                "AND idx_outer.run_id_b >= ? " // dim 1
                                "AND idx_outer.run_id_a <= ? " // dim 1
                                "AND idx_outer.maxX >= ? " // dim 2
                                "AND idx_outer.minX <= ? " // dim 2
                                "AND idx_outer.maxY >= ? " // dim 3
                                "AND idx_outer.minY <= ? " // dim 3
                                "AND outer.switch_strand = ? "
                                "LIMIT 1 " );
                            SQLQuery<DBCon, int64_t> xNumOverlapsHelper2(
                                xNewConn,
                                // make sure that inner does not overlap with any other call with higher score
                                "SELECT inner2.id "
                                "FROM sv_call_table AS inner2, sv_call_r_tree AS idx_inner2 "
                                "WHERE inner2.id == idx_inner2.id "
                                "AND idx_inner2.id != ? "
                                // tuple comparison here, so that overlapping calls with equal score always have
                                // one call take priority (in this case always the one inserted after)
                                "AND (" +
                                    SvCallTable::getSqlForCallScore( "inner2" ) +
                                    ", inner2.id) > (?, ?) "
                                    "AND idx_inner2.run_id_b >= ? " // dim 1
                                    "AND idx_inner2.run_id_a <= ? " // dim 1
                                    "AND idx_inner2.maxX >= ? " // dim 2
                                    "AND idx_inner2.minX <= ? " // dim 2
                                    "AND idx_inner2.maxY >= ? " // dim 3
                                    "AND idx_inner2.minY <= ? " // dim 3
                                    "AND inner2.switch_strand == ? "
                                    "LIMIT 1 " );
                            for( size_t uiTupId = uiJobId_; uiTupId < vResults.size( ); uiTupId += 32 )
                            {
                                auto& rTup = vResults[ uiTupId ];
                                int64_t iId = std::get<0>( rTup );
                                double dScore = std::get<1>( rTup );
                                uint32_t uiFromStart = std::get<2>( rTup );
                                uint32_t uiFromSize = std::get<3>( rTup );
                                uint32_t uiToStart = std::get<4>( rTup );
                                uint32_t uiToSize = std::get<5>( rTup );
                                bool bSwitchStrand = std::get<6>( rTup );

                                // get the maximal blur for which the current call does not overlap with another one
                                int64_t iBlurMax = 0;
                                while( true )
                                {
                                    // check if current call overlaps a call with higher score of the same run
                                    auto xIt2 = xNumOverlapsHelper2.vExecuteAndReturnIterator(
                                        iId, dScore, iId, iCallerRunIdB, iCallerRunIdB, uiFromStart - iBlurMax,
                                        uiFromStart + uiFromSize + iBlurMax, uiToStart - iBlurMax,
                                        uiToStart + uiToSize + iBlurMax, bSwitchStrand );
                                    if( !xIt2.eof( ) )
                                        break;
                                    iBlurMax++;
                                    if( iBlurMax >= iAllowedDist )
                                        break;
                                } // while

                                // call immediately overlaps with call with higher score
                                if( iBlurMax == 0 )
                                    continue;

                                // get the minimal blur for which the current call does overlap a ground truth
                                int64_t iBlurMin = iAllowedDist;
                                while( true )
                                {
                                    auto xIt1 = xNumOverlapsHelper1.vExecuteAndReturnIterator(
                                        iCallerRunIdA, iCallerRunIdA, uiFromStart - iBlurMin,
                                        uiFromStart + uiFromSize + iBlurMin, uiToStart - iBlurMin,
                                        uiToStart + uiToSize + iBlurMin, bSwitchStrand );
                                    if( xIt1.eof( ) )
                                        break;
                                    iBlurMin--;
                                    if( iBlurMin < 0 )
                                        break;
                                } // while
                                iBlurMin++;

                                // call does not overlap with ground truth at all
                                if( iBlurMin > iAllowedDist )
                                    continue;

                                // fill the cache with the connections
                                std::lock_guard<std::mutex> xGuard( *pWriteLock );
                                xCacheTable.xInsertRow( dScore, iBlurMin, iBlurMax );
                            } // for
                        },
                        uiJobId );
            } // scope of xPool
#else // code replacement or new SQL API:
            {
                const size_t uiNumThreads = 1; // In MySQL, one thread allowed only.
                ThreadPool xPool( uiNumThreads ); // DESIGN: We require a database connection pool ...

                for( size_t uiJobId = 0; uiJobId < uiNumThreads; uiJobId++ )
                    xPool.enqueue(
                        [ & ]( size_t, size_t uiJobId_ ) {
                            // SQLDB<DBCon> xNewConn(sDBName);
                            SQLQuery<DBCon, int64_t> xNumOverlapsHelper1(
                                pDatabase,
                                // make sure that inner overlaps the outer:
                                "SELECT outer.id "
                                "FROM sv_call_table AS outer, sv_call_r_tree AS idx_outer "
                                "WHERE outer.id == idx_outer.id "
                                "AND idx_outer.run_id_b >= ? " // dim 1
                                "AND idx_outer.run_id_a <= ? " // dim 1
                                "AND idx_outer.maxX >= ? " // dim 2 -> spatial geometry dim 1
                                "AND idx_outer.minX <= ? " // dim 2 -> spatial geometry dim 1
                                "AND idx_outer.maxY >= ? " // dim 3 -> spatial geometry dim 2
                                "AND idx_outer.minY <= ? " // dim 3 -> spatial geometry dim 2
                                "AND outer.switch_strand = ? "
                                "LIMIT 1 " );
                            SQLQuery<DBCon, int64_t> xNumOverlapsHelper2(
                                pDatabase,
                                // make sure that inner does not overlap with any other call with higher score
                                "SELECT inner2.id "
                                "FROM sv_call_table AS inner2, sv_call_r_tree AS idx_inner2 "
                                "WHERE inner2.id == idx_inner2.id "
                                "AND idx_inner2.id != ? "
                                // tuple comparison here, so that overlapping calls with equal score always have
                                // one call take priority (in this case always the one inserted after)
                                "AND (" +
                                    SvCallTable::getSqlForCallScore( "inner2" ) +
                                    ", inner2.id) > (?, ?) "
                                    "AND idx_inner2.run_id_b >= ? " // dim 1
                                    "AND idx_inner2.run_id_a <= ? " // dim 1
                                    "AND idx_inner2.maxX >= ? " // dim 2 -> spatial geometry dim 1
                                    "AND idx_inner2.minX <= ? " // dim 2 -> spatial geometry dim 1
                                    "AND idx_inner2.maxY >= ? " // dim 3 -> spatial geometry dim 2
                                    "AND idx_inner2.minY <= ? " // dim 3 -> spatial geometry dim 2
                                    "AND inner2.switch_strand = ? "
                                    "LIMIT 1 " );
                            for( size_t uiTupId = uiJobId_; uiTupId < vResults.size( ); uiTupId += uiNumThreads )
                            {
                                auto& rTup = vResults[ uiTupId ];
                                int64_t iId = std::get<0>( rTup );
                                double dScore = std::get<1>( rTup );
                                uint32_t uiFromStart = std::get<2>( rTup );
                                uint32_t uiFromSize = std::get<3>( rTup );
                                uint32_t uiToStart = std::get<4>( rTup );
                                uint32_t uiToSize = std::get<5>( rTup );
                                bool bSwitchStrand = std::get<6>( rTup );

                                // get the maximal blur for which the current call does not overlap with another one
                                int64_t iBlurMax = 0;

                                while( true )
                                {
                                    // check if current call overlaps a call with higher score of the same run
                                    // DEL:auto xIt2 = xNumOverlapsHelper2.vExecuteAndReturnIterator(
                                    // DEL:    iId, dScore, iId, iCallerRunIdB, iCallerRunIdB, uiFromStart - iBlurMax,
                                    // DEL:    uiFromStart + uiFromSize + iBlurMax, uiToStart - iBlurMax,
                                    // DEL:    uiToStart + uiToSize + iBlurMax, bSwitchStrand );
                                    xNumOverlapsHelper2.execAndFetch(
                                        iId, dScore, iId, iCallerRunIdB, iCallerRunIdB, uiFromStart - iBlurMax,
                                        uiFromStart + uiFromSize + iBlurMax, uiToStart - iBlurMax,
                                        uiToStart + uiToSize + iBlurMax, bSwitchStrand );

                                    // DEL: if( !xIt2.eof( ) )
                                    if( !xNumOverlapsHelper2.eof( ) )
                                        break;
                                    iBlurMax++;
                                    if( iBlurMax >= iAllowedDist )
                                        break;
                                } // while

                                // call immediately overlaps with call with higher score
                                if( iBlurMax == 0 )
                                    continue;

                                // get the minimal blur for which the current call does overlap a ground truth
                                int64_t iBlurMin = iAllowedDist;
                                while( true )
                                {
                                    // DEL:auto xIt1 = xNumOverlapsHelper1.vExecuteAndReturnIterator(
                                    // DEL:    iCallerRunIdA, iCallerRunIdA, uiFromStart - iBlurMin,
                                    // DEL:    uiFromStart + uiFromSize + iBlurMin, uiToStart - iBlurMin,
                                    // DEL:    uiToStart + uiToSize + iBlurMin, bSwitchStrand );
                                    xNumOverlapsHelper1.execAndFetch(
                                        iCallerRunIdA, iCallerRunIdA, uiFromStart - iBlurMin,
                                        uiFromStart + uiFromSize + iBlurMin, uiToStart - iBlurMin,
                                        uiToStart + uiToSize + iBlurMin, bSwitchStrand );

                                    // DEL: if( xIt1.eof( ) )
                                    if( xNumOverlapsHelper1.eof( ) )
                                        break;
                                    iBlurMin--;
                                    if( iBlurMin < 0 )
                                        break;
                                } // while
                                iBlurMin++;

                                // call does not overlap with ground truth at all
                                if( iBlurMin > iAllowedDist )
                                    continue;

                                // fill the cache with the connections
                                std::lock_guard<std::mutex> xGuard( *pWriteLock );
                                // DEL: xCacheTable.xInsertRow( dScore, iBlurMin, iBlurMax );
								xCacheTable.insert(dScore, iBlurMin, iBlurMax);
                            } // for
                        },
                        uiJobId );
            } // scope of xPool
#endif
#if 0
// py code

import sqlite3
conn = sqlite3.connect("/mnt/ssd1/sv_datasets2/del_human/svs.db")
c = conn.cursor()
c.execute("DELETE FROM overlap_cache_computed_table")
c.execute("DROP TABLE overlap_cache_table_9_1")
c.execute("DELETE FROM overlap_cache_table")
conn.commit()
conn.close()
exit()

#endif
            // remember that we have computed the cache
            pHelper->insert( iCallerRunIdA, iCallerRunIdB, iAllowedDist ); // xInsertRow->insert
            // DEL: xCacheTable.xCreateIndex.exec( );
            xCacheTable.addIndex(
                json{ { INDEX_NAME, xCacheTable.sTableName + "_index" }, { INDEX_COLUMNS, "score_a" } } );
        } // method
        /**
         * @brief helper function; see libMA::SvCallTable::numOverlaps
         * @details
         * iCallerRunIdA used as ground truth
         */
        inline uint32_t numOverlaps( int64_t iCallerRunIdA, int64_t iCallerRunIdB, double dMinScore,
                                     int64_t iAllowedDist )
        {
            metaMeasureAndLogDuration<false>( "createCache", [ & ]( ) { //
                createCache( iCallerRunIdA, iCallerRunIdB, iAllowedDist ); //
            } );
            uint32_t uiRet;
            metaMeasureAndLogDuration<false>( "xNumOverlapsCache.scalar", [ & ]( ) { //
                OverlapCacheTable xCacheTable( pDatabase, iCallerRunIdA, iCallerRunIdB );
                // xNumOverlapsCache.bindAndExplain(iCallerRunIdB, iCallerRunIdA, dMinScore, iAllowedDist);
                uiRet = xCacheTable.xNumOverlapsCache.scalar( dMinScore, iAllowedDist, iAllowedDist );
            } );
            return uiRet;
        } // method

        /**
         * @brief returns how many calls are invalid because they overlap another call with higher score
         */
        inline uint32_t numInvalidCalls( int64_t iCallerRunIdA, double dMinScore, int64_t iAllowedDist )
        {
            uint32_t uiRet = 0;
#if 0
            xNumOverlaps.vExecuteAndForAllRowsUnpackedDo(
                [&]( int64_t iId, double dScore, uint32_t uiFromStart, uint32_t uiFromSize, uint32_t uiToStart,
                     uint32_t uiToSize, bool bSwitchStrand ) {
                    if( xNumOverlapsHelper2
                            .vExecuteAndReturnIterator(
                                iCallerRunIdA, dScore, iCallerRunIdA, uiFromStart - iAllowedDist,
                                uiFromStart + uiFromSize + iAllowedDist, uiToStart - iAllowedDist,
                                uiToStart + uiToSize + iAllowedDist, bSwitchStrand )
                            .eof( ) )
                        return;
                    uiRet += 1;
                },
                iCallerRunIdA, dMinScore );
#endif

            return uiRet;
        } // method
    }; // class

#ifndef WITHOUT_R_TREE
    RTreeIndex xIndex; // merely required for calling the constructor
#else
    std::shared_ptr<RTreeIndex> xIndex;
#endif
    // Was originally CppSQLiteExtInsertStatement statement.
    // FIXME: Add insert statement SQLStatement<DBCon, int64_t, int64_t, int64_t, uint32_t, uint32_t, uint32_t,
    // uint32_t> xInsertRTree;
    SQLQuery<DBCon, uint32_t> xQuerySize;
    SQLQuery<DBCon, uint32_t> xQuerySizeSpecific;
    SQLQuery<DBCon, int64_t> xCallArea;
    SQLQuery<DBCon, double> xMaxScore;
    SQLQuery<DBCon, double> xMinScore;
    SQLQuery<DBCon, int64_t, bool, uint32_t, uint32_t, NucSeqSql, uint32_t> xNextCallForwardContext;
    SQLQuery<DBCon, int64_t, bool, uint32_t, uint32_t, NucSeqSql, uint32_t> xNextCallBackwardContext;
    SQLStatement<DBCon> xSetCoverageForCall;
    SQLStatement<DBCon> xDeleteCall1, xDeleteCall2;
    SQLStatement<DBCon> xUpdateCall, xUpdateRTree;
    std::shared_ptr<OverlapCache> pOverlapCache;

  public:
    // Consider: Place the table on global level
    json jSvCallTableDef( )
    {
        return json{
            { TABLE_NAME, "sv_call_table" },
            { TABLE_COLUMNS,
              {
                  { { COLUMN_NAME, "sv_caller_run_id" } },
                  { { COLUMN_NAME, "from_pos" } },
                  { { COLUMN_NAME, "to_pos" } },
                  { { COLUMN_NAME, "from_size" } },
                  { { COLUMN_NAME, "to_size" } },
                  { { COLUMN_NAME, "switch_strand" } },
                  { { COLUMN_NAME, "inserted_sequence" } },
                  { { COLUMN_NAME, "supporting_reads" } },
                  { { COLUMN_NAME, "reference_ambiguity" } },
                  { { COLUMN_NAME, "regex_id" } },
              } }, //
            { FOREIGN_KEY,
              { { COLUMN_NAME, "sv_caller_run_id" }, { REFERENCES, "sv_caller_run_table(id) ON DELETE CASCADE" } } },
            { FOREIGN_KEY,
              { { COLUMN_NAME, "regex_id" }, { REFERENCES, "sv_call_reg_ex_table(id) ON DELETE SET NULL" } } } };
    }; // method

    SvCallTable( std::shared_ptr<DBCon> pDatabase, std::shared_ptr<std::mutex> pWriteLock,
                 std::string sDBName )
        : SvCallTableType<DBCon>( pDatabase, // the database where the table resides
                                  jSvCallTableDef( ) ), // table definition

          pDatabase( pDatabase ),
#ifndef WITHOUT_R_TREE
          // pReconstructionTable( std::make_shared<ReconstructionTable>( pDatabase ) ),
          xIndex( pDatabase ), // create R-tree (R-tree table view)
#else
          xIndex( std::make_shared<RTreeIndex>( pDatabase ) ),
#endif

          // FIXME: Add insert statement xInsertRTree( *pDatabase, "sv_call_r_tree", false ),
          xQuerySize( pDatabase, "SELECT COUNT(*) FROM sv_call_table" ),

          xQuerySizeSpecific( pDatabase, "SELECT COUNT(*) FROM sv_call_table, sv_call_r_tree "
                                         "WHERE sv_call_table.id = sv_call_r_tree.id "
                                         "AND sv_call_r_tree.run_id_a >= ? " // dim 1
                                         "AND sv_call_r_tree.run_id_b <= ? " // dim 1
                                         "AND " +
                                             getSqlForCallScore( ) + " >= ? " ),
          xCallArea( pDatabase,
                     "SELECT SUM( from_size * to_size ) FROM sv_call_table, sv_call_r_tree "
                     "WHERE sv_call_table.id = sv_call_r_tree.id "
                     "AND sv_call_r_tree.run_id_a >= ? " // dim 1
                     "AND sv_call_r_tree.run_id_b <= ? " // dim 1
                     "AND " +
                         getSqlForCallScore( ) + " >= ? " ),
          xMaxScore( pDatabase,
                     "SELECT " + getSqlForCallScore( ) +
                         " FROM sv_call_table, sv_call_r_tree "
                         "WHERE sv_call_table.id = sv_call_r_tree.id "
                         "AND sv_call_r_tree.run_id_a >= ? " // dim 1
                         "AND sv_call_r_tree.run_id_b <= ? " // dim 1
                         "ORDER BY " +
                         getSqlForCallScore( ) + " DESC LIMIT 1 " ),
          xMinScore( pDatabase,
                     "SELECT " + getSqlForCallScore( ) +
                         " FROM sv_call_table, sv_call_r_tree "
                         "WHERE sv_call_table.id = sv_call_r_tree.id "
                         "AND sv_call_r_tree.run_id_a >= ? " // dim 1
                         "AND sv_call_r_tree.run_id_b <= ? " // dim 1
                         "ORDER BY " +
                         getSqlForCallScore( ) + " ASC LIMIT 1 " ),

          xNextCallForwardContext(
              pDatabase,
              "SELECT sv_call_table.id, switch_strand, to_pos, to_size, inserted_sequence, from_pos + from_size "
              "FROM sv_call_table "
              "WHERE sv_call_table.sv_caller_run_id = ? " // dim 1
              "AND sv_call_table.from_pos >= ? " // dim 2
#if 0
                  "AND NOT EXISTS ( "
                  "  SELECT * "
                  "  FROM reconstruction_table "
                  "  WHERE sv_call_table.id == reconstruction_table.call_id "
                  ") "
#endif
              "ORDER BY sv_call_table.from_pos ASC "
              "LIMIT 1 " ),
          xNextCallBackwardContext( pDatabase,
                                    "SELECT sv_call_table.id, switch_strand, from_pos, from_size, "
                                    "       inserted_sequence, to_pos "
                                    "FROM sv_call_table "
                                    "WHERE sv_call_table.sv_caller_run_id = ? " // dim 1
                                    "AND sv_call_table.to_pos <= ? " // dim 2
#if 0
                                        "AND NOT EXISTS ( "
                                        "  SELECT * "
                                        "  FROM reconstruction_table "
                                        "  WHERE sv_call_table.id == reconstruction_table.call_id "
                                        ") "
#endif
                                    "ORDER BY sv_call_table.to_pos DESC "
                                    "LIMIT 1 " ),
          xSetCoverageForCall( pDatabase,
                               "UPDATE sv_call_table "
                               "SET reference_ambiguity = ? "
                               "WHERE id = ?" ),
          xDeleteCall1( pDatabase,
                        "DELETE FROM sv_call_r_tree "
                        "WHERE id = ? " ),
          xDeleteCall2( pDatabase,
                        "DELETE FROM sv_call_table "
                        "WHERE id = ? " ),
          xUpdateCall( pDatabase,
                       "UPDATE sv_call_table "
                       "SET from_pos = ?, "
                       "    to_pos = ?, "
                       "    from_size = ?, "
                       "    to_size = ?, "
                       "    switch_strand = ?, "
                       "    inserted_sequence = ?, "
                       "    supporting_reads = ?, "
                       "    reference_ambiguity = ? "
                       "WHERE id = ? " ),
          xUpdateRTree( pDatabase,
                        "UPDATE sv_call_r_tree "
                        "SET run_id_a = ?, "
                        "    run_id_b = ?, "
                        "    minX = ?, "
                        "    maxX = ?, "
                        "    minY = ?, "
                        "    maxY = ? "
                        "WHERE id = ? " ),
          pOverlapCache( std::make_shared<OverlapCache>( pDatabase, pWriteLock, sDBName ) )
    {} // default constructor

    inline void addScoreIndex( int64_t iCallerRunId )
    {
        // Discuss Markus: This index definition is somehow defect ...
        // DEL:SQLStatement<DBCon>( pDatabase,
        // DEL:                     ( "CREATE INDEX IF NOT EXISTS sv_call_table_score_index_" +
        // DEL:                       std::to_string( iCallerRunId ) + " ON sv_call_table ( " + getSqlForCallScore( ) +
        // DEL:                       " )                         " +
        // DEL:                       " WHERE sv_caller_run_id = " + std::to_string( iCallerRunId ) ) DEL
        // DEL:                     :.c_str( ) )
        // DEL:    .exec( );
#if 0 // Must be fixed together with Markus		
        this->addIndex( json{ { INDEX_NAME, "sv_call_table_score_index_" + std::to_string(iCallerRunId) },
                              { INDEX_COLUMNS, getSqlForCallScore( ) },
                              { WHERE, "sv_caller_run_id = " + std::to_string( iCallerRunId ) } } );
#endif
    } // method

    inline uint32_t numCalls( )
    {
        return xQuerySize.scalar( );
    } // method

    inline uint32_t numCalls( int64_t iCallerRunId, double dMinScore )
    {
        return xQuerySizeSpecific.scalar( iCallerRunId, iCallerRunId, dMinScore );
    } // method

    inline void updateCoverage( SvCall& rCall )
    {
        xSetCoverageForCall.exec( (uint32_t)rCall.uiReferenceAmbiguity, rCall.iId );
    } // method

    inline void deleteCall( int64_t iCallId )
    {
        xDeleteCall1.exec( iCallId );
        xDeleteCall2.exec( iCallId );
    } // method

    inline void deleteCall( SvCall& rCall )
    {
        deleteCall( rCall.iId );
    } // method

    inline int64_t insertCall( int64_t iSvCallerRunId, SvCall& rCall )
    {
        // DEL: int64_t iCallId = //
        // DEL:     this->xInsertRow( iSvCallerRunId, (uint32_t)rCall.uiFromStart, (uint32_t)rCall.uiToStart,
        // DEL:                       (uint32_t)rCall.uiFromSize, (uint32_t)rCall.uiToSize,
        // DEL:                       rCall.bSwitchStrand, // NucSeqSql can deal
        // DEL:                       with nullpointers NucSeqSql( rCall.pInsertedSequence ), //
        // DEL:                       (uint32_t)rCall.uiNumSuppReads, (uint32_t)rCall.uiReferenceAmbiguity, -1 );
        int64_t iCallId = this->insert( iSvCallerRunId, //
                                        (uint32_t)rCall.uiFromStart, //
                                        (uint32_t)rCall.uiToStart, //
                                        (uint32_t)rCall.uiFromSize, (uint32_t)rCall.uiToSize, rCall.bSwitchStrand,
                                        // NucSeqSql can deal with nullpointers
                                        NucSeqSql( rCall.pInsertedSequence ), (uint32_t)rCall.uiNumSuppReads,
                                        (uint32_t)rCall.uiReferenceAmbiguity, -1 );

        rCall.iId = iCallId;

        // DEL: xInsertRTree( iCallId, iSvCallerRunId, iSvCallerRunId, (uint32_t)rCall.uiFromStart,
        // DEL:               (uint32_t)rCall.uiFromStart + (uint32_t)rCall.uiFromSize, (uint32_t)rCall.uiToStart,
        // DEL:               (uint32_t)rCall.uiToStart + (uint32_t)rCall.uiToSize );
        // DEBUG: std::cout << "R-Tree insert:" << std::endl;
        xIndex->insert( static_cast<int32_t>( iCallId ), //
                        static_cast<int32_t>( iSvCallerRunId ), //
                        (uint32_t)iSvCallerRunId, (uint32_t)rCall.uiFromStart,
                        (uint32_t)rCall.uiFromStart + (uint32_t)rCall.uiFromSize, (uint32_t)rCall.uiToStart,
                        (uint32_t)rCall.uiToStart + (uint32_t)rCall.uiToSize );
        return iCallId;
    } // method

    inline int64_t updateCall( int64_t iSvCallerRunId, SvCall& rCall )
    {
        xUpdateCall.exec( (uint32_t)rCall.uiFromStart, (uint32_t)rCall.uiToStart, (uint32_t)rCall.uiFromSize,
                          (uint32_t)rCall.uiToSize, rCall.bSwitchStrand,
                          // NucSeqSql can deal with nullpointers
                          NucSeqSql( rCall.pInsertedSequence ), (uint32_t)rCall.uiNumSuppReads,
                          (uint32_t)rCall.uiReferenceAmbiguity, rCall.iId );
        xUpdateRTree.exec( iSvCallerRunId, iSvCallerRunId, (uint32_t)rCall.uiFromStart,
                           (uint32_t)rCall.uiFromStart + (uint32_t)rCall.uiFromSize, (uint32_t)rCall.uiToStart,
                           (uint32_t)rCall.uiToStart + (uint32_t)rCall.uiToSize, rCall.iId );
        return rCall.iId;
    } // method

    inline int64_t callArea( int64_t iCallerRunId, double dMinScore )
    {
        return xCallArea.scalar( iCallerRunId, iCallerRunId, dMinScore );
    } // method

    inline double maxScore( int64_t iCallerRunId )
    {
        return xMaxScore.scalar( iCallerRunId, iCallerRunId );
    } // method

    inline double minScore( int64_t iCallerRunId )
    {
        return xMinScore.scalar( iCallerRunId, iCallerRunId );
    } // method

    /**
     * @brief returns how many calls of run A are overlapped by a call in run B
     * @details
     * Only considers calls of run A with score >= to dMinScore.
     * Calls that are no further away than iAllowedDist are considered overlapping (can be used to add some fuzziness).
     * If two calls in run A overlap, only the one with higher score counts; If both have the same score the one with
     * the higher id is kept.
     */
    inline uint32_t numOverlaps( int64_t iCallerRunIdA, int64_t iCallerRunIdB, double dMinScore, int64_t iAllowedDist )
    {
        return pOverlapCache->numOverlaps( iCallerRunIdA, iCallerRunIdB, dMinScore, iAllowedDist );
    } // method

    /**
     * @brief returns the average distance of class from the overlapped (due to fuzziness) SV
     */
    inline double blurOnOverlaps( int64_t iCallerRunIdA, int64_t iCallerRunIdB, double dMinScore, int64_t iAllowedDist )
    {
        int64_t uiSum = 0;
        int64_t uiCount = 0;
        for( int64_t iI = 0; iI <= iAllowedDist; iI++ )
        {
            uint32_t uiAmount = numOverlaps( iCallerRunIdA, iCallerRunIdB, dMinScore, iI );
            uiSum += uiAmount * iI;
            uiCount += uiAmount;
        } // for
        return uiSum / (double)uiCount;
    } // method

    /**
     * @brief returns how many calls are invalid because they overlap another call with higher score
     */
    inline uint32_t numInvalidCalls( int64_t iCallerRunIdA, double dMinScore, int64_t iAllowedDist )
    {
        return pOverlapCache->numInvalidCalls( iCallerRunIdA, dMinScore, iAllowedDist );
    } // method

    /// @brief returns call id, jump start pos, next context, next from position, jump end position
    inline std::tuple<int64_t, uint32_t, bool, NucSeqSql, uint32_t> getNextCall( int64_t iCallerRun, //
                                                                                 uint32_t uiFrom, //
                                                                                 bool bForwardContext )
    {
        std::tuple<int64_t, uint32_t, bool, NucSeqSql, uint32_t> xRet;
        std::get<0>( xRet ) = -1;
        std::get<2>( xRet ) = bForwardContext; // does nothing...
        if( bForwardContext )
        {
            xNextCallForwardContext.execAndFetch( iCallerRun, uiFrom );
            if( !xNextCallForwardContext.eof( ) )
            {
                auto xQ = xNextCallForwardContext.get( );
                std::get<0>( xRet ) = std::get<0>( xQ );
                std::get<1>( xRet ) = std::get<5>( xQ );
                std::get<2>( xRet ) = !std::get<1>( xQ );
                std::get<3>( xRet ) = std::get<4>( xQ );
                std::get<4>( xRet ) = std::get<1>( xQ ) ? std::get<2>( xQ ) + std::get<3>( xQ ) : std::get<2>( xQ );
            } // if
        } // if
        else
        {
            xNextCallBackwardContext.execAndFetch( iCallerRun, uiFrom );
            if( !xNextCallForwardContext.eof( ) )
            {
                auto xQ = xNextCallForwardContext.get( );
                std::get<0>( xRet ) = std::get<0>( xQ );
                std::get<1>( xRet ) = std::get<5>( xQ );
                std::get<2>( xRet ) = std::get<1>( xQ );
                std::get<3>( xRet ) = std::get<4>( xQ );
                std::get<4>( xRet ) = !std::get<1>( xQ ) ? std::get<2>( xQ ) + std::get<3>( xQ ) : std::get<2>( xQ );
            } // if
        } // else
        return xRet;
    } // method

    inline std::shared_ptr<Pack> reconstructSequencedGenome( std::shared_ptr<Pack> pRef, int64_t iCallerRun )
    {
        {
            // Discuss Markus: Why from_pos + from_size -> Results in SQL error
            // DEL:SQLStatement<DBCon>( pDatabase,
            // DEL:                     ( "CREATE INDEX IF NOT EXISTS tmp_reconstruct_seq_index_1_" +
            // DEL:                       std::to_string( iCallerRun ) +
            // DEL:                       " ON sv_call_table (from_pos, id, switch_strand, to_pos, to_size, "
            // DEL:                       "                  inserted_sequence, from_pos + from_size) "
            // DEL:                       "WHERE sv_caller_run_id = " +
            // DEL:                       std::to_string( iCallerRun ) )
            // DEL:                         .c_str( ) )
            // DEL:    .exec( );

            // Discuss Markus: inserted_sequence is of type NucSeqSql and so a LONGBLOB. This kind of unlimited
            // BLOB's or text can not be used as part of an index in MySQL. See:
            // https://stackoverflow.com/questions/1827063/mysql-error-key-specification-without-a-key-length
            this->addIndex( json{ { INDEX_NAME, "tmp_reconstruct_seq_index_1_" + std::to_string( iCallerRun ) },
                                  { INDEX_COLUMNS, "from_pos, id, switch_strand, to_pos, to_size, "
                                                   "from_size" },
                                  { WHERE, "sv_caller_run_id = " + std::to_string( iCallerRun ) } } );

            // DEL:SQLStatement<DBCon>( pDatabase,
            // DEL:                     ( "CREATE INDEX IF NOT EXISTS tmp_reconstruct_seq_index_2_" +
            // DEL:                       std::to_string( iCallerRun ) +
            // DEL:                       " ON sv_call_table (to_pos, id, switch_strand, from_pos, from_size, "
            // DEL:                       "                  inserted_sequence) "
            // DEL:                       "WHERE sv_caller_run_id = " +
            // DEL:                       std::to_string( iCallerRun ) )
            // DEL:                         .c_str( ) )
            // DEL:    .exec( );
            this->addIndex( json{ { INDEX_NAME, "tmp_reconstruct_seq_index_2_" + std::to_string( iCallerRun ) },
                                  { INDEX_COLUMNS, "to_pos, id, switch_strand, from_pos, from_size" },
                                  { WHERE, "sv_caller_run_id = " + std::to_string( iCallerRun ) } } );
        } // scope for CppSQLiteExtStatement

        // @todo at the moment this does not deal with jumped over sequences
        // @todo at the moment this does not check the regex (?)
        auto pRet = std::make_shared<Pack>( );

        std::set<int64_t> xVisitedCalls;

        NucSeq xCurrChrom;
        uint32_t uiCurrPos = 0;
        uint32_t uiContigCnt = 1;
        bool bForwContext = true;
        while( true )
        {
            // get the next call
            std::tuple<int64_t, uint32_t, bool, NucSeqSql, uint32_t> tNextCall;
            uint32_t uiIntermediatePos = uiCurrPos;
            do
            {
                // search for the next call that we have not visited yet...
                metaMeasureAndLogDuration<false>(
                    "SQL", [ & ]( ) { tNextCall = this->getNextCall( iCallerRun, uiIntermediatePos, bForwContext ); } );
                if( std::get<0>( tNextCall ) == -1 || // there is no next call
                                                      // we have not visited the next call
                    xVisitedCalls.find( std::get<0>( tNextCall ) ) == xVisitedCalls.end( ) )
                    break;
                // we have visited the next call and need to search again

                // @todo this is extremely inefficient (if we have cylces in our graph which we do not at the moment)
                uiIntermediatePos += bForwContext ? 1 : -1;
            } while( true );
#if 0
            std::cout << "id: " << std::get<0>( tNextCall ) << " from: " << std::get<1>( tNextCall )
                      << " to: " << std::get<4>( tNextCall ) << ( std::get<2>( tNextCall ) ? " forward" : " rev-comp" )
                      << std::endl;
#endif
            if( std::get<0>( tNextCall ) == -1 ) // if there are no more calls
            {
                metaMeasureAndLogDuration<false>( "seq copy final", [ & ]( ) {
                    // for jumps to the end of the genome we do not want to extract the last contig...
                    // this check becomes necessary since the with the current index system,
                    // we would either extract the last nucleotide of the genome twice or extract the
                    // reverse complement of the last contig...
                    if( pRef->uiUnpackedSizeForwardStrand != uiCurrPos )
                        pRef->vExtractContext( uiCurrPos, xCurrChrom, true, bForwContext );

                    pRet->vAppendSequence( "unnamed_contig_" + std::to_string( uiContigCnt++ ), "no_description_given",
                                           xCurrChrom );
                    xCurrChrom.vClear( );

                    // part of the if above...
                    if( pRef->uiUnpackedSizeForwardStrand == uiCurrPos )
                        return;

                    /*
                     * for this we make use of the id system of contigs.
                     * the n forwards contigs have the ids: x*2 | 0 <= x <= n
                     * the n reverse complement contigs have the ids: x*2+1 | 0 <= x <= n
                     */
                    for( int64_t uiI = pRef->uiSequenceIdForPositionOrRev( uiCurrPos ) + ( bForwContext ? 2 : -1 );
                         uiI < (int64_t)pRef->uiNumContigs( ) * 2 && uiI >= 0;
                         uiI += ( bForwContext ? 2 : -2 ) )
                    {
                        pRef->vExtractContig( uiI, xCurrChrom, true );
                        pRet->vAppendSequence( "unnamed_contig_" + std::to_string( uiContigCnt++ ),
                                               "no_description_given", xCurrChrom );
                        xCurrChrom.vClear( );
                    } // for
                } );
                break;
            } // if

            // we reach this point if there are more calls, so tNextCall is set properly here
            metaMeasureAndLogDuration<false>( "seq copy", [ & ]( ) {
                // if the next call is in a different chromosome
                while( pRef->bridgingPositions( uiCurrPos, std::get<1>( tNextCall ) ) )
                {
                    // extract the remaining chromosome into xCurrChrom
                    uiCurrPos = (uint32_t)pRef->vExtractContext( uiCurrPos, xCurrChrom, true, bForwContext );
                    // append xCurrChrom to the pack
                    pRet->vAppendSequence( "unnamed_contig_" + std::to_string( uiContigCnt++ ), "no_description_given",
                                           xCurrChrom );
                    // clear xCurrChrom
                    xCurrChrom.vClear( );
                    // if the next call is several chromosomes over this loops keeps going
                } // while
                // the call is in the current chromosome / we have appended all skipped chromosomes
                if( bForwContext )
                    pRef->vExtractSubsectionN( uiCurrPos, std::get<1>( tNextCall ) + 1, xCurrChrom, true );
                else
                    pRef->vExtractSubsectionN( pRef->uiPositionToReverseStrand( uiCurrPos ),
                                               pRef->uiPositionToReverseStrand( std::get<1>( tNextCall ) ) + 1, //
                                               xCurrChrom,
                                               true );
                // append the skipped over sequence
                if( std::get<3>( tNextCall ).pNucSeq != nullptr )
                    xCurrChrom.vAppend( std::get<3>( tNextCall ).pNucSeq->pxSequenceRef,
                                        std::get<3>( tNextCall ).pNucSeq->length( ) );

                metaMeasureAndLogDuration<false>( "xInsertRow", [ & ]( ) {
                    // remember that we used this call
                    // pReconstructionTable->xInsertRow( std::get<0>( tNextCall ) );
                    xVisitedCalls.insert( std::get<0>( tNextCall ) );
                    bForwContext = std::get<2>( tNextCall );
                    uiCurrPos = std::get<4>( tNextCall );
                } );
            } );
        } // while

        // clear up the temp table
        // pReconstructionTable->clearTable( );
        // @todo hmm eventually the extra indices should be cleared up?
        //       however, these do only exist for call sets where we reconstruct the genome
        //       i.e. the ground truth data set
        // CppSQLiteExtStatement( *pDatabase, "DROP INDEX tmp_reconstruct_seq_index_1 " ).execDML( );
        // CppSQLiteExtStatement( *pDatabase, "DROP INDEX tmp_reconstruct_seq_index_2 " ).execDML( );

        return pRet;
    } // method
}; // namespace libMA

} // namespace libMA
