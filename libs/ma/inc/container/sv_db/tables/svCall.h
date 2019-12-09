/**
 * @file sequencer.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "db_sql.h"
#include <csignal>

namespace libMA
{

typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int64_t, // sv_caller_run_id (foreign key)
                                                 uint32_t, // from_pos
                                                 uint32_t, // to_pos
                                                 uint32_t, // from_size
                                                 uint32_t, // to_size
                                                 bool, // switch_strand
                                                 NucSeqSql, // inserted_sequence
                                                 uint32_t, // supporting_nt
                                                 uint32_t, // coverage
                                                 int64_t // regex_id
                                                 >
    TP_SV_CALL_TABLE;
class SvCallTable : private TP_SV_CALL_TABLE
{

#if 0
        typedef CppSQLiteExtTable<int64_t // call_id
                                  >
            TP_RECONSTRUCTION_TABLE;
        class ReconstructionTable : public TP_RECONSTRUCTION_TABLE
        {
            std::shared_ptr<CppSQLiteDBExtended> pDatabase;

          public:
            ReconstructionTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
                : TP_RECONSTRUCTION_TABLE( *pDatabase, // the database where the table resides
                                           "reconstruction_table", // name of the table in the database
                                           // column definitions of the table
                                           std::vector<std::string>{"call_id"},
                                           false ),
                  pDatabase( pDatabase )
            {
                pDatabase->execDML( "CREATE INDEX IF NOT EXISTS reconstruction_table_index "
                                    "ON reconstruction_table (call_id) " );
            } // default constructor

        }; // class
#endif

    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    // std::shared_ptr<ReconstructionTable> pReconstructionTable;
    class RTreeIndex
    {
      public:
        RTreeIndex( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
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
    RTreeIndex xIndex; // merely required for calling the constructor
    CppSQLiteExtInsertStatement<int64_t, int64_t, int64_t, uint32_t, uint32_t, uint32_t, uint32_t> xInsertRTree;
    CppSQLiteExtQueryStatement<uint32_t> xQuerySize;
    CppSQLiteExtQueryStatement<uint32_t> xQuerySizeSpecific;
    CppSQLiteExtQueryStatement<int64_t, double, uint32_t, uint32_t, uint32_t, uint32_t, bool> xNumOverlaps;
    CppSQLiteExtQueryStatement<int64_t> xNumOverlapsHelper1;
    CppSQLiteExtQueryStatement<int64_t> xNumOverlapsHelper2;
    CppSQLiteExtQueryStatement<int64_t> xCallArea;
    CppSQLiteExtQueryStatement<double> xMaxScore;
    CppSQLiteExtQueryStatement<double> xMinScore;
    CppSQLiteExtQueryStatement<int64_t, bool, uint32_t, uint32_t, NucSeqSql, uint32_t> xNextCallForwardContext;
    CppSQLiteExtQueryStatement<int64_t, bool, uint32_t, uint32_t, NucSeqSql, uint32_t> xNextCallBackwardContext;
    CppSQLiteExtStatement xSetCoverageForCall;
    CppSQLiteExtStatement xDeleteCall1, xDeleteCall2;
    CppSQLiteExtStatement xUpdateCall, xUpdateRTree;

  public:
    SvCallTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
        : TP_SV_CALL_TABLE( *pDatabase, // the database where the table resides
                            "sv_call_table", // name of the table in the database
                            // column definitions of the table
                            std::vector<std::string>{"sv_caller_run_id", "from_pos", "to_pos", "from_size", "to_size",
                                                     "switch_strand", "inserted_sequence", "supporting_nt", "coverage",
                                                     "regex_id"},
                            // constraints for table
                            std::vector<std::string>{
                                "FOREIGN KEY (sv_caller_run_id) REFERENCES sv_caller_run_table(id) ON DELETE CASCADE",
                                "FOREIGN KEY (regex_id) REFERENCES sv_call_reg_ex_table(id) ON DELETE SET NULL"} ),
          pDatabase( pDatabase ),
          // pReconstructionTable( std::make_shared<ReconstructionTable>( pDatabase ) ),
          xIndex( pDatabase ),
          xInsertRTree( *pDatabase, "sv_call_r_tree", false ),
          xQuerySize( *pDatabase, "SELECT COUNT(*) FROM sv_call_table" ),
          xQuerySizeSpecific( *pDatabase, "SELECT COUNT(*) FROM sv_call_table, sv_call_r_tree "
                                          "WHERE sv_call_table.id == sv_call_r_tree.id "
                                          "AND sv_call_r_tree.run_id_a >= ? " // dim 1
                                          "AND sv_call_r_tree.run_id_b <= ? " // dim 1
                                          "AND (supporting_nt*1.0)/coverage >= ? " ),
          xNumOverlaps( *pDatabase,
                        // each inner call can overlap an outer call at most once
                        "SELECT id, supporting_nt*1.0/coverage, from_pos, from_size, to_pos, to_size, "
                        "       switch_strand "
                        "FROM sv_call_table "
                        "WHERE sv_caller_run_id = ? "
                        "AND supporting_nt*1.0/coverage >= ? " ),
          xNumOverlapsHelper1( *pDatabase,
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
                               "AND outer.switch_strand == ? "
                               "LIMIT 1 " ),
          xNumOverlapsHelper2( *pDatabase,
                               // make sure that inner does not overlap with any other call with higher score
                               "SELECT inner2.id "
                               "FROM sv_call_table AS inner2, sv_call_r_tree AS idx_inner2 "
                               "WHERE inner2.id == idx_inner2.id "
                               "AND idx_inner2.id != ? "
                               "AND (inner2.supporting_nt*1.0)/inner2.coverage >= ? "
                               "AND idx_inner2.run_id_b >= ? " // dim 1
                               "AND idx_inner2.run_id_a <= ? " // dim 1
                               "AND idx_inner2.maxX >= ? " // dim 2
                               "AND idx_inner2.minX <= ? " // dim 2
                               "AND idx_inner2.maxY >= ? " // dim 3
                               "AND idx_inner2.minY <= ? " // dim 3
                               "AND inner2.switch_strand == ? "
                               "LIMIT 1 " ),
          xCallArea( *pDatabase,
                     "SELECT SUM( from_size * to_size ) FROM sv_call_table, sv_call_r_tree "
                     "WHERE sv_call_table.id == sv_call_r_tree.id "
                     "AND sv_call_r_tree.run_id_a >= ? " // dim 1
                     "AND sv_call_r_tree.run_id_b <= ? " // dim 1
                     "AND (supporting_nt*1.0)/coverage >= ? " ),
          xMaxScore( *pDatabase,
                     "SELECT supporting_nt*1.0/coverage FROM sv_call_table, sv_call_r_tree "
                     "WHERE sv_call_table.id == sv_call_r_tree.id "
                     "AND sv_call_r_tree.run_id_a >= ? " // dim 1
                     "AND sv_call_r_tree.run_id_b <= ? " // dim 1
                     "ORDER BY (supporting_nt*1.0)/coverage DESC LIMIT 1 " ),
          xMinScore( *pDatabase,
                     "SELECT (supporting_nt*1.0)/coverage FROM sv_call_table, sv_call_r_tree "
                     "WHERE sv_call_table.id == sv_call_r_tree.id "
                     "AND sv_call_r_tree.run_id_a >= ? " // dim 1
                     "AND sv_call_r_tree.run_id_b <= ? " // dim 1
                     "ORDER BY (supporting_nt*1.0)/coverage ASC LIMIT 1 " ),
          xNextCallForwardContext(
              *pDatabase,
              "SELECT sv_call_table.id, switch_strand, to_pos, to_size, inserted_sequence, from_pos + from_size "
              "FROM sv_call_table "
              "WHERE sv_call_table.sv_caller_run_id == ? " // dim 1
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
          xNextCallBackwardContext( *pDatabase,
                                    "SELECT sv_call_table.id, switch_strand, from_pos, from_size, "
                                    "       inserted_sequence, to_pos "
                                    "FROM sv_call_table "
                                    "WHERE sv_call_table.sv_caller_run_id == ? " // dim 1
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
          xSetCoverageForCall( *pDatabase,
                               "UPDATE sv_call_table "
                               "SET coverage = ? "
                               "WHERE id == ?" ),
          xDeleteCall1( *pDatabase,
                        "DELETE FROM sv_call_r_tree "
                        "WHERE id == ? " ),
          xDeleteCall2( *pDatabase,
                        "DELETE FROM sv_call_table "
                        "WHERE id == ? " ),
          xUpdateCall( *pDatabase,
                       "UPDATE sv_call_table "
                       "SET from_pos = ?, "
                       "    to_pos = ?, "
                       "    from_size = ?, "
                       "    to_size = ?, "
                       "    switch_strand = ?, "
                       "    inserted_sequence = ?, "
                       "    supporting_nt = ?, "
                       "    coverage = ? "
                       "WHERE id == ? " ),
          xUpdateRTree( *pDatabase,
                        "UPDATE sv_call_r_tree "
                        "SET run_id_a = ?, "
                        "    run_id_b = ?, "
                        "    minX = ?, "
                        "    maxX = ?, "
                        "    minY = ?, "
                        "    maxY = ? "
                        "WHERE id == ? " )
    {} // default constructor

    inline void addScoreIndex( int64_t iCallerRunId )
    {
        CppSQLiteExtStatement( *pDatabase,
                               ( "CREATE INDEX IF NOT EXISTS sv_call_table_score_index_" +
                                 std::to_string( iCallerRunId ) +
                                 " ON sv_call_table ((supporting_nt*1.0)/coverage) "
                                 "WHERE sv_caller_run_id == " +
                                 std::to_string( iCallerRunId ) )
                                   .c_str( ) )
            .execDML( );
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
        xSetCoverageForCall.bindAndExecute( (uint32_t)rCall.uiCoverage, rCall.iId );
    } // method

    inline void deleteCall( int64_t iCallId )
    {
        xDeleteCall1.bindAndExecute( iCallId );
        xDeleteCall2.bindAndExecute( iCallId );
    } // method

    inline void deleteCall( SvCall& rCall )
    {
        deleteCall( rCall.iId );
    } // method

    inline int64_t insertCall( int64_t iSvCallerRunId, SvCall& rCall )
    {
        int64_t iCallId = this->xInsertRow( iSvCallerRunId, (uint32_t)rCall.uiFromStart, (uint32_t)rCall.uiToStart,
                                            (uint32_t)rCall.uiFromSize, (uint32_t)rCall.uiToSize, rCall.bSwitchStrand,
                                            // NucSeqSql can deal with nullpointers
                                            NucSeqSql( rCall.pInsertedSequence ), (uint32_t)rCall.uiNumSuppNt,
                                            (uint32_t)rCall.uiCoverage, -1 );
        rCall.iId = iCallId;
        xInsertRTree( iCallId, iSvCallerRunId, iSvCallerRunId, (uint32_t)rCall.uiFromStart,
                      (uint32_t)rCall.uiFromStart + (uint32_t)rCall.uiFromSize, (uint32_t)rCall.uiToStart,
                      (uint32_t)rCall.uiToStart + (uint32_t)rCall.uiToSize );
        return iCallId;
    } // method

    inline int64_t updateCall( int64_t iSvCallerRunId, SvCall& rCall )
    {
        xUpdateCall.bindAndExecute( (uint32_t)rCall.uiFromStart, (uint32_t)rCall.uiToStart, (uint32_t)rCall.uiFromSize,
                                    (uint32_t)rCall.uiToSize, rCall.bSwitchStrand,
                                    // NucSeqSql can deal with nullpointers
                                    NucSeqSql( rCall.pInsertedSequence ), (uint32_t)rCall.uiNumSuppNt,
                                    (uint32_t)rCall.uiCoverage, rCall.iId );
        xUpdateRTree.bindAndExecute( iSvCallerRunId, iSvCallerRunId, (uint32_t)rCall.uiFromStart,
                                     (uint32_t)rCall.uiFromStart + (uint32_t)rCall.uiFromSize,
                                     (uint32_t)rCall.uiToStart, (uint32_t)rCall.uiToStart + (uint32_t)rCall.uiToSize,
                                     rCall.iId );
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
     * returns how many calls of run A are overlapped by a call in run B
     */
    inline uint32_t numOverlaps( int64_t iCallerRunIdA, int64_t iCallerRunIdB, double dMinScore, int64_t iAllowedDist )
    {
        // uint32_t uiNumCalls = numCalls( iCallerRunIdA, 0 ) * 3;
        uint32_t uiRet = 0;
        auto vResults = xNumOverlaps.executeAndStoreAllInVector( iCallerRunIdB, dMinScore );
        for( auto xTup : vResults )
        {
            int64_t iId = std::get<0>( xTup );
            double dScore = std::get<1>( xTup );
            uint32_t uiFromStart = std::get<2>( xTup );
            uint32_t uiFromSize = std::get<3>( xTup );
            uint32_t uiToStart = std::get<4>( xTup );
            uint32_t uiToSize = std::get<5>( xTup );
            bool bSwitchStrand = std::get<6>( xTup );

            if( xNumOverlapsHelper1
                    .vExecuteAndReturnIterator( iCallerRunIdA, iCallerRunIdA, uiFromStart - iAllowedDist,
                                                uiFromStart + uiFromSize + iAllowedDist, uiToStart - iAllowedDist,
                                                uiToStart + uiToSize + iAllowedDist, bSwitchStrand )
                    .eof( ) )
                continue;
            if( !xNumOverlapsHelper2
                     .vExecuteAndReturnIterator( iId, dScore, iCallerRunIdB, iCallerRunIdB, uiFromStart - iAllowedDist,
                                                 uiFromStart + uiFromSize + iAllowedDist, uiToStart - iAllowedDist,
                                                 uiToStart + uiToSize + iAllowedDist, bSwitchStrand )
                     .eof( ) )
                continue;
            uiRet += 1;
        } // for

        return uiRet;
    } // method

    /**
     * returns the average distance of class from the overlapped (due to fuzziness) SV
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
     * returns how many calls are invalid because they overlap another call with higher score
     */
    inline uint32_t numInvalidCalls( int64_t iCallerRunIdA, double dMinScore, int64_t iAllowedDist )
    {
        uint32_t uiRet = 0;
        xNumOverlaps.vExecuteAndForAllRowsUnpackedDo(
            [&]( int64_t iId, double dScore, uint32_t uiFromStart, uint32_t uiFromSize, uint32_t uiToStart,
                 uint32_t uiToSize, bool bSwitchStrand ) {
                if( xNumOverlapsHelper2
                        .vExecuteAndReturnIterator( iCallerRunIdA, dScore, uiFromStart - iAllowedDist,
                                                    uiFromStart + uiFromSize + iAllowedDist, uiToStart - iAllowedDist,
                                                    uiToStart + uiToSize + iAllowedDist, bSwitchStrand )
                        .eof( ) )
                    return;
                uiRet += 1;
            },
            iCallerRunIdA, dMinScore );

        return uiRet;
    } // method

    // returns call id, jump start pos, next context, next from position, jump end position
    inline std::tuple<int64_t, uint32_t, bool, NucSeqSql, uint32_t> getNextCall( int64_t iCallerRun, //
                                                                                 uint32_t uiFrom, //
                                                                                 bool bForwardContext )
    {
        std::tuple<int64_t, uint32_t, bool, NucSeqSql, uint32_t> xRet;
        std::get<0>( xRet ) = -1;
        std::get<2>( xRet ) = bForwardContext; // does nothing...
        if( bForwardContext )
        {
            auto itQ = xNextCallForwardContext.vExecuteAndReturnIterator( iCallerRun, uiFrom );
            if( !itQ.eof( ) )
            {
                auto xQ = itQ.get( );
                std::get<0>( xRet ) = std::get<0>( xQ );
                std::get<1>( xRet ) = std::get<5>( xQ );
                std::get<2>( xRet ) = !std::get<1>( xQ );
                std::get<3>( xRet ) = std::get<4>( xQ );
                std::get<4>( xRet ) = std::get<1>( xQ ) ? std::get<2>( xQ ) + std::get<3>( xQ ) : std::get<2>( xQ );
            } // if
        } // if
        else
        {
            auto itQ = xNextCallBackwardContext.vExecuteAndReturnIterator( iCallerRun, uiFrom );
            if( !itQ.eof( ) )
            {
                auto xQ = itQ.get( );
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
            CppSQLiteExtStatement( *pDatabase,
                                   ( "CREATE INDEX IF NOT EXISTS tmp_reconstruct_seq_index_1_" +
                                     std::to_string( iCallerRun ) +
                                     " ON sv_call_table (from_pos, id, switch_strand, to_pos, to_size, "
                                     "                  inserted_sequence, from_pos + from_size) "
                                     "WHERE sv_caller_run_id == " +
                                     std::to_string( iCallerRun ) )
                                       .c_str( ) )
                .execDML( );
            CppSQLiteExtStatement( *pDatabase,
                                   ( "CREATE INDEX IF NOT EXISTS tmp_reconstruct_seq_index_2_" +
                                     std::to_string( iCallerRun ) +
                                     " ON sv_call_table (to_pos, id, switch_strand, from_pos, from_size, "
                                     "                  inserted_sequence) "
                                     "WHERE sv_caller_run_id == " +
                                     std::to_string( iCallerRun ) )
                                       .c_str( ) )
                .execDML( );
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
                    "SQL", [&]( ) { tNextCall = this->getNextCall( iCallerRun, uiIntermediatePos, bForwContext ); } );
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
                metaMeasureAndLogDuration<false>( "seq copy final", [&]( ) {
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
            metaMeasureAndLogDuration<false>( "seq copy", [&]( ) {
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

                metaMeasureAndLogDuration<false>( "xInsertRow", [&]( ) {
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
}; // class

} // namespace libMA