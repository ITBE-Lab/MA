/**
 * @file svCall.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "container/pack.h"
#include "container/svJump.h"
#include "geom.h"
#include "sql_api.h"
#include "system.h"
#include "threadPool.h"
#include "wkb_spatial.h"
#include <csignal>
#include <string>

namespace libMA
{
template <typename DBCon>
using SvCallTableType = SQLTableWithLibIncrPriKey<DBCon, // DB connector type
                                                  PriKeyDefaultType, // sv_caller_run_id (foreign key)
                                                  uint32_t, // from_pos (geometry)
                                                  uint32_t, // to_pos (geometry)
                                                  uint32_t, // from_size (geometry)
                                                  uint32_t, // to_size (geometry)
                                                  bool, // switch_strand
                                                  std::shared_ptr<CompressedNucSeq>, // inserted_sequence
                                                  uint32_t, // supporting_reads
                                                  uint32_t, // reference_ambiguity
                                                  int64_t, // regex_id
                                                  WKBUint64Rectangle // rectangle (geometry)
                                                  >;

template <typename DBCon> class SvCallTable : public SvCallTableType<DBCon>
{

    std::shared_ptr<DBCon> pConnection;
    SQLQuery<DBCon, uint32_t> xQuerySize;
    SQLQuery<DBCon, uint32_t> xQuerySizeSpecific;
    SQLQuery<DBCon, int64_t> xCallArea;
    SQLQuery<DBCon, double> xMaxScore;
    SQLQuery<DBCon, double> xMinScore;
    SQLStatement<DBCon> xSetCoverageForCall;
    SQLStatement<DBCon> xDeleteCall;
    SQLStatement<DBCon> xUpdateCall;
    SQLQuery<DBCon, WKBUint64Rectangle, bool, double, PriKeyDefaultType> xNumOverlaps;
    SQLQuery<DBCon, uint32_t> xHelperIntersecCallWHigherScore;
    SQLQuery<DBCon, uint32_t> xHelperIntersectingCall;
    // std::shared_ptr<OverlapCache> pOverlapCache;

  public:
    // Consider: Place the table on global level
    json jSvCallTableDef( )
    {
        return json{
            {TABLE_NAME, "sv_call_table"},
            {TABLE_COLUMNS,
             {{{COLUMN_NAME, "sv_caller_run_id"}},
              {{COLUMN_NAME, "from_pos"}},
              {{COLUMN_NAME, "to_pos"}},
              {{COLUMN_NAME, "from_size"}},
              {{COLUMN_NAME, "to_size"}},
              {{COLUMN_NAME, "switch_strand"}},
              {{COLUMN_NAME, "inserted_sequence"}},
              {{COLUMN_NAME, "supporting_reads"}},
              {{COLUMN_NAME, "reference_ambiguity"}},
              {{COLUMN_NAME, "regex_id"}},
              {{COLUMN_NAME, "rectangle"}, {CONSTRAINTS, "NOT NULL"}}}},
            {GENERATED_COLUMNS,
             {{{COLUMN_NAME, "score"}, {TYPE, "DOUBLE"}, {AS, "(supporting_reads * 1.0) / reference_ambiguity"}}}},
            // @todo how to insert NULL into sv_caller_run_id ?
            // @todo ask arne about inserting NULL
            /*{FOREIGN_KEY,
             {{COLUMN_NAME, "sv_caller_run_id"}, {REFERENCES, "sv_caller_run_table(id) ON DELETE CASCADE"}}},*/
            // @todo regex_id is unused at the moment
            /*{FOREIGN_KEY, {{COLUMN_NAME, "regex_id"}, {REFERENCES, "sv_call_reg_ex_table(id) ON DELETE SET NULL"}}}*/
        };
    }; // method

    SvCallTable( std::shared_ptr<DBCon> pConnection )
        : SvCallTableType<DBCon>( pConnection, // the database where the table resides
                                  jSvCallTableDef( ) ), // table definition
          pConnection( pConnection ),
          xQuerySize( pConnection, "SELECT COUNT(*) FROM sv_call_table" ),

          xQuerySizeSpecific( pConnection, "SELECT COUNT(*) FROM sv_call_table "
                                           "WHERE sv_caller_run_id = ? "
                                           "AND score >= ? " ),
          xCallArea( pConnection,
                     "SELECT SUM( from_size * to_size ) FROM sv_call_table "
                     "WHERE sv_caller_run_id = ? "
                     "AND score >= ? " ),
          xMaxScore( pConnection,
                     "SELECT score "
                     " FROM sv_call_table "
                     "WHERE sv_caller_run_id = ? "
                     "ORDER BY score DESC LIMIT 1 " ),
          xMinScore( pConnection,
                     "SELECT score "
                     " FROM sv_call_table "
                     "WHERE sv_caller_run_id = ? "
                     "ORDER BY score ASC LIMIT 1 " ),
          xSetCoverageForCall( pConnection,
                               "UPDATE sv_call_table "
                               "SET reference_ambiguity = ? "
                               "WHERE id = ? " ),
          xDeleteCall( pConnection,
                       "DELETE FROM sv_call_table "
                       "WHERE id = ? " ),
          xUpdateCall( pConnection,
                       "UPDATE sv_call_table "
                       "SET from_pos = ?, "
                       "    to_pos = ?, "
                       "    from_size = ?, "
                       "    to_size = ?, "
                       "    switch_strand = ?, "
                       "    inserted_sequence = ?, "
                       "    supporting_reads = ?, "
                       "    reference_ambiguity = ?, "
                       "    rectangle = ST_PolyFromWKB(?, 0) "
                       "WHERE id = ? " ),
          xNumOverlaps( pConnection,
                        "SELECT ST_AsBinary(rectangle), switch_strand, score, id "
                        "FROM sv_call_table "
                        "WHERE sv_caller_run_id = ? "
                        "AND score >= ? ",
                        json{}, "SvCallTable::xNumOverlaps" ),
          xHelperIntersecCallWHigherScore( pConnection,
                                           "SELECT COUNT(*) "
                                           "FROM sv_call_table "
                                           "WHERE sv_caller_run_id = ? "
                                           "AND MBRIntersects(rectangle, ST_PolyFromWKB(?, 0)) "
                                           "AND switch_strand = ? "
                                           "AND (score, id) > (?, ?) "
                                           "LIMIT 1 ",
                                           json{}, "SvCallTable::xHelperIntersecCallWHigherScore" ),
          xHelperIntersectingCall( pConnection,
                                   "SELECT COUNT(*) "
                                   "FROM sv_call_table "
                                   "WHERE sv_caller_run_id = ? "
                                   "AND MBRIntersects(rectangle, ST_PolyFromWKB(?, 0)) "
                                   "AND switch_strand = ? "
                                   "LIMIT 1 ",
                                   json{}, "SvCallTable::xHelperIntersectingCall" )
    //,pOverlapCache( std::make_shared<OverlapCache>( pDatabase, pWriteLock, sDBName ) ) //@todo delete this..
    {} // default constructor


    /**
     * @brief returns how many calls of run A are overlapped by a call in run B
     * @details
     * Only considers calls of run A with score >= to dMinScore.
     * Calls that are no further away than iAllowedDist are considered overlapping (can be used to add some
     * fuzziness). If two calls in run A overlap, only the one with higher score counts; If both have the same score
     * the one with the higher id is kept.
     * UAAAAGGGH mysql is too stupid to use the rectangle spatial index if the queries in here are combined...
     * splitting them up results in the desired behaviour.
     * Even with FORCE INDEX hints the optimizer insists on making a full table scan
     */
    inline uint32_t numOverlaps( int64_t iCallerRunIdA, int64_t iCallerRunIdB, double dMinScore, int64_t iAllowedDist )
    {
        uint32_t uiRet = 0;
        // @todo parallelize this
        xNumOverlaps.execAndForAll(
            [this, &uiRet, iCallerRunIdA, iCallerRunIdB, iAllowedDist]( WKBUint64Rectangle xWKB, bool bSwitchStrand,
                                                                        double fScore, PriKeyDefaultType iId ) {
                auto xRect = xWKB.getRect( );
                xRect.resize( iAllowedDist );
                WKBUint64Rectangle xWKBResized( xRect );
                if( xHelperIntersecCallWHigherScore.scalar( iCallerRunIdA, xWKBResized, bSwitchStrand, fScore, iId ) ==
                    0 )
                    if( xHelperIntersectingCall.scalar( iCallerRunIdB, xWKBResized, bSwitchStrand ) > 0 )
                        uiRet++;
            },
            iCallerRunIdA, dMinScore );
        return uiRet;
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
        uint32_t uiRet = 0;
        xNumOverlaps.execAndForAll(
            [this, &uiRet, iCallerRunIdA, iAllowedDist]( WKBUint64Rectangle xWKB, bool bSwitchStrand, double fScore,
                                                         PriKeyDefaultType iId ) {
                auto xRect = xWKB.getRect( );
                xRect.resize( iAllowedDist );
                WKBUint64Rectangle xWKBResized( xRect );
                if( xHelperIntersecCallWHigherScore.scalar( iCallerRunIdA, xWKBResized, bSwitchStrand, fScore, iId ) >
                    0 )
                    uiRet++;
            },
            iCallerRunIdA, dMinScore );
        return uiRet;
    } // method

    inline void genIndices( int64_t iCallerRunId )
    {
        this->addIndex( json{{INDEX_NAME, "rectangle"}, {INDEX_COLUMNS, "rectangle"}, {INDEX_TYPE, "SPATIAL"}} );

        // see: https://dev.mysql.com/doc/refman/5.7/en/create-table-generated-columns.html
        // and: https://dev.mysql.com/doc/refman/5.7/en/create-table-secondary-indexes.html
        this->addIndex( json{{INDEX_NAME, "runId_score"}, {INDEX_COLUMNS, "sv_caller_run_id, score"}} );
    } // method

    inline void dropIndices( int64_t iCallerRunId )
    {
        this->dropIndex( json{{INDEX_NAME, "rectangle"}} );
        this->dropIndex( json{{INDEX_NAME, "runId_score"}} );
    } // method

    inline uint32_t numCalls( )
    {
        return xQuerySize.scalar( );
    } // method

    inline uint32_t numCalls( int64_t iCallerRunId, double dMinScore )
    {
        return xQuerySizeSpecific.scalar( iCallerRunId, dMinScore );
    } // method
    inline uint32_t numCalls_py( int64_t iCallerRunId, double dMinScore )
    {
        return numCalls( iCallerRunId, dMinScore );
    } // method

    inline void updateCoverage( SvCall& rCall )
    {
        xSetCoverageForCall.exec( (uint32_t)rCall.uiReferenceAmbiguity, rCall.iId );
    } // method

    inline void deleteCall( int64_t iCallId )
    {
        xDeleteCall.exec( iCallId );
    } // method

    inline void deleteCall( SvCall& rCall )
    {
        deleteCall( rCall.iId );
    } // method

    inline int64_t insertCall( int64_t iSvCallerRunId, SvCall& rCall )
    {
        auto xRectangle = WKBUint64Rectangle( rCall );
        int64_t iCallId =
            this->insert( iSvCallerRunId, //
                          (uint32_t)rCall.xXAxis.start( ), //
                          (uint32_t)rCall.xYAxis.start( ), //
                          (uint32_t)rCall.xXAxis.size( ), //
                          (uint32_t)rCall.xYAxis.size( ), //
                          rCall.bSwitchStrand,
                          // can deal with nullpointers
                          makeSharedCompNucSeq( rCall.pInsertedSequence ), //
                          (uint32_t)rCall.uiNumSuppReads, (uint32_t)rCall.uiReferenceAmbiguity, -1, xRectangle );
        rCall.iId = iCallId;

        return iCallId;
    } // method

    inline int64_t updateCall( int64_t iSvCallerRunId, SvCall& rCall )
    {
        auto xRectangle = WKBUint64Rectangle( rCall );
        xUpdateCall.exec( (uint32_t)rCall.xXAxis.start( ), (uint32_t)rCall.xYAxis.start( ),
                          (uint32_t)rCall.xXAxis.size( ), (uint32_t)rCall.xYAxis.size( ), rCall.bSwitchStrand,
                          // can deal with nullpointers
                          makeSharedCompNucSeq( rCall.pInsertedSequence ), (uint32_t)rCall.uiNumSuppReads,
                          (uint32_t)rCall.uiReferenceAmbiguity, xRectangle, rCall.iId );
        return rCall.iId;
    } // method

    inline int64_t callArea( int64_t iCallerRunId, double dMinScore )
    {
        return xCallArea.scalar( iCallerRunId, dMinScore );
    } // method

    inline double maxScore( int64_t iCallerRunId )
    {
        return xMaxScore.scalar( iCallerRunId );
    } // method

    inline double minScore( int64_t iCallerRunId )
    {
        return xMinScore.scalar( iCallerRunId );
    } // method


    using NextCallType =
        SQLQuery<DBCon, int64_t, bool, uint32_t, uint32_t, std::shared_ptr<CompressedNucSeq>, uint32_t>;

    /** @brief returns call id, jump start pos, next context, next from position, jump end position
     *  @details helper function for reconstructSequencedGenome */
    inline std::tuple<int64_t, uint32_t, bool, std::shared_ptr<NucSeq>, uint32_t>
    getNextCall( int64_t iCallerRun, //
                 uint32_t uiFrom, //
                 bool bForwardContext,
                 NextCallType& xNextCallForwardContext,
                 NextCallType& xNextCallBackwardContext )
    {
        std::tuple<int64_t, uint32_t, bool, std::shared_ptr<NucSeq>, uint32_t> xRet;
        std::get<0>( xRet ) = -1;
        std::get<2>( xRet ) = bForwardContext; // does nothing...
        if( bForwardContext && xNextCallForwardContext.execAndFetch( uiFrom ) )
        {
            auto xQ = xNextCallForwardContext.get( );
            std::get<0>( xRet ) = std::get<0>( xQ );
            std::get<1>( xRet ) = std::get<5>( xQ );
            std::get<2>( xRet ) = !std::get<1>( xQ );
            std::get<3>( xRet ) = std::get<4>( xQ ) == nullptr ? nullptr : std::get<4>( xQ )->pUncomNucSeq;
            std::get<4>( xRet ) = std::get<1>( xQ ) ? std::get<2>( xQ ) + std::get<3>( xQ ) : std::get<2>( xQ );
        } // if
        else if( !bForwardContext && xNextCallBackwardContext.execAndFetch( uiFrom ) )
        {
            auto xQ = xNextCallBackwardContext.get( );
            std::get<0>( xRet ) = std::get<0>( xQ );
            std::get<1>( xRet ) = std::get<5>( xQ );
            std::get<2>( xRet ) = std::get<1>( xQ );
            std::get<3>( xRet ) = std::get<4>( xQ ) == nullptr ? nullptr : std::get<4>( xQ )->pUncomNucSeq;
            std::get<4>( xRet ) = !std::get<1>( xQ ) ? std::get<2>( xQ ) + std::get<3>( xQ ) : std::get<2>( xQ );
        } // if
        return xRet;
    } // method

    /** @brief reconstruct a sequenced genome from a reference and the calls of the run with id iCallerRun.
     *  @details
     *  @todo at the moment this does not deal with jumped over sequences
     *  @todo at the moment this does not check the regex (?)
     *  Creates a reconstruction_table that is filled with all unused calls from iCallerRun and then deletes the calls
     *  one by one until the sequenced genome is reconstructed
     */
    inline std::shared_ptr<Pack> reconstructSequencedGenome( std::shared_ptr<Pack> pRef, PriKeyDefaultType iCallerRun )
    {
        auto pRet = std::make_shared<Pack>( );
        auto pTransaction = this->pConnection->uniqueGuardedTrxn( );
        SQLTable<DBCon, PriKeyDefaultType, uint32_t, uint32_t> xReconstructionTable(
            this->pConnection,
            json{{TABLE_NAME, "reconstruction_table"},
                 {CPP_EXTRA, "DROP ON DESTRUCTION"},
                 {TABLE_COLUMNS,
                  {
                      {{COLUMN_NAME, "call_id"},
                       {REFERENCES, "sv_call_table(id)"},
                       {CONSTRAINTS, "NOT NULL PRIMARY KEY"}},
                      {{COLUMN_NAME, "from_pos"}},
                      {{COLUMN_NAME, "to_pos"}}
                      //
                  }}} );

        SQLStatement<DBCon> xInsert( this->pConnection, "INSERT INTO reconstruction_table (call_id, from_pos, to_pos) "
                                                        "SELECT id, from_pos, to_pos "
                                                        "FROM sv_call_table "
                                                        "WHERE sv_caller_run_id = ? " );
        SQLStatement<DBCon> xDelete( this->pConnection, "DELETE FROM reconstruction_table "
                                                        "WHERE call_id = ? " );

        metaMeasureAndLogDuration<false>( "fill reconstruction table", [&]( ) { xInsert.exec( iCallerRun ); } );

        metaMeasureAndLogDuration<false>( "create indices on reconstruction table", [&]( ) {
            xReconstructionTable.addIndex( json{{INDEX_NAME, "tmp_rct_from"}, {INDEX_COLUMNS, "from_pos, call_id"}} );
            xReconstructionTable.addIndex( json{{INDEX_NAME, "tmp_rct_to"}, {INDEX_COLUMNS, "to_pos, call_id"}} );
        } );

        NextCallType xNextCallForwardContext(
            this->pConnection,
            "SELECT id, switch_strand, sv_call_table.to_pos, to_size, inserted_sequence, "
            "       sv_call_table.from_pos + from_size "
            "FROM sv_call_table "
            "INNER JOIN reconstruction_table ON reconstruction_table.call_id = sv_call_table.id "
            "WHERE reconstruction_table.from_pos >= ? "
            "ORDER BY reconstruction_table.from_pos ASC "
            "LIMIT 1 " );
        NextCallType xNextCallBackwardContext(
            this->pConnection,
            "SELECT id, switch_strand, sv_call_table.from_pos, from_size, "
            "       inserted_sequence, sv_call_table.to_pos "
            "FROM sv_call_table "
            "INNER JOIN reconstruction_table ON reconstruction_table.call_id = sv_call_table.id "
            "WHERE reconstruction_table.to_pos <= ? "
            "ORDER BY reconstruction_table.to_pos DESC "
            "LIMIT 1 " );

#if DEBUG_LEVEL > 0
        std::set<int64_t> xVisitedCalls;
#endif

        NucSeq xCurrChrom;
        uint32_t uiCurrPos = 0;
        uint32_t uiContigCnt = 1;
        bool bForwContext = true;
        while( true )
        {
            // get the next call
            std::tuple<int64_t, uint32_t, bool, std::shared_ptr<NucSeq>, uint32_t> tNextCall;
            uint32_t uiIntermediatePos = uiCurrPos;
            do
            {
                // search for the next call that we have not visited yet...
                metaMeasureAndLogDuration<false>( "SQL", [&]( ) {
                    tNextCall = this->getNextCall( iCallerRun, uiIntermediatePos, bForwContext, xNextCallForwardContext,
                                                   xNextCallBackwardContext );
                } );
#if DEBUG_LEVEL > 0
                if( std::get<0>( tNextCall ) == -1 || // there is no next call
                                                      // we have not visited the next call
                    xVisitedCalls.find( std::get<0>( tNextCall ) ) == xVisitedCalls.end( ) )
#endif
                    break;
                // we have visited the next call and need to search again
                std::cout << "SHOULD NEVER REACH THIS PRINT?" << std::endl;
                assert( false );

                // this is extremely inefficient (if we have cylces in our graph which we do not at the
                // moment)
                uiIntermediatePos += bForwContext ? 1 : -1;
            } while( true );
#if 0
            std::cout << "id: " << std::get<0>( tNextCall ) << " from: " << std::get<1>( tNextCall )
                      << " to: " << std::get<4>( tNextCall ) << ( std::get<2>( tNextCall ) ? " forward" : " rev-comp" )
                      << " inserted_seq: "
                      << ( std::get<3>( tNextCall ) == nullptr
                               ? "nullptr"
                               : ( std::get<3>( tNextCall )->uiSize == 0
                                       ? "empty"
                                       : std::get<3>( tNextCall )->toString( ) ) )
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
                if( std::get<3>( tNextCall ) != nullptr && std::get<3>( tNextCall )->uiSize > 0 )
                    xCurrChrom.vAppend( std::get<3>( tNextCall )->pxSequenceRef, std::get<3>( tNextCall )->length( ) );

                metaMeasureAndLogDuration<false>( "xInsertRow", [&]( ) {
                    // remember that we used this call
                    xDelete.exec( std::get<0>( tNextCall ) );
#if DEBUG_LEVEL > 0
                    xVisitedCalls.insert( std::get<0>( tNextCall ) );
#endif
                    bForwContext = std::get<2>( tNextCall );
                    uiCurrPos = std::get<4>( tNextCall );
                } );
            } );
        } // while

        return pRet;
    } // method
}; // namespace libMA

} // namespace libMA
