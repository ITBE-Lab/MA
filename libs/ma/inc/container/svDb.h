/**
 * @file svDb.h
 * @details
 * The database interface for the structural variant caller
 */

#include "container/container.h"
#include "container/nucSeq.h"
#include "container/pack.h"
#include "container/soc.h"
#include "container/svJump.h"
#include "module/module.h"
#include "util/exception.h"
#include "util/sqlite3.h"
#include <chrono>
#include <ctime>
#include <iomanip>

namespace libMA
{

class SV_DB : public Container
{
  private:
    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string // sequencer name
                                                     >
        TP_SEQUENCER_TABLE;
    class SequencerTable : public TP_SEQUENCER_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;

      public:
        SequencerTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SEQUENCER_TABLE( *pDatabase, // the database where the table resides
                                  "sequencer_table", // name of the table in the database
                                  // column definitions of the table
                                  std::vector<std::string>{"name"},
                                  // constraints for table
                                  std::vector<std::string>{"UNIQUE (name)"} ),
              pDatabase( pDatabase )
        {} // default constructor

        inline void createIndices( )
        {
            pDatabase->execDML( "CREATE INDEX sequencer_id_index ON sequencer_table (id)" );
        } // method

        inline int64_t insertSequencer( std::string& sSequencerName )
        {
            return xInsertRow( sSequencerName );
        } // method
    }; // class

    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int64_t, // sequencer id (foreign key)
                                                     std::string, // read name
                                                     NucSeqSql // read sequence
                                                     >
        TP_READ_TABLE;
    class ReadTable : public TP_READ_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        bool bDoDuplicateWarning = true;

      public:
        CppSQLiteExtQueryStatement<int32_t> xGetReadId;

        ReadTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_READ_TABLE( *pDatabase, // the database where the table resides
                             "read_table", // name of the table in the database
                             // column definitions of the table
                             std::vector<std::string>{"sequencer_id", "name", "sequence"},
                             // constraints for table
                             std::vector<std::string>{"UNIQUE (sequencer_id, name)", //
                                                      "FOREIGN KEY (sequencer_id) REFERENCES sequencer_table(id)"} ),
              pDatabase( pDatabase ),
              xGetReadId( *pDatabase, "SELECT id FROM read_table WHERE sequencer_id = ? AND name = ?" )
        {} // default constructor

        inline void createIndices( )
        {
            pDatabase->execDML( "CREATE INDEX read_id_index ON read_table (id)" );
        } // method

        inline int64_t insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pRead )
        {
            for( size_t i = 0; i < 5; i++ ) // at most do 5 tries
            {
                try
                {
                    return xInsertRow( uiSequencerId, pRead->sName, NucSeqSql( pRead ) );
                } // try
                catch( CppSQLite3Exception& xException )
                {
                    if( bDoDuplicateWarning )
                    {
                        std::cerr << "WARNING: " << xException.errorMessage( ) << std::endl;
                        std::cerr << "Does your data contain duplicate reads? Current read name: " << pRead->sName
                                  << std::endl;
                        std::cerr << "Changing read name to: " << pRead->sName + "_2" << std::endl;
                        std::cerr << "This warning is only displayed once" << std::endl;
                        bDoDuplicateWarning = false;
                    } // if
                    pRead->sName += "_2";
                } // catch
            } // for
            throw AnnotatedException( "Could not insert read after 5 tries" );
        } // method
    }; // class

    typedef CppSQLiteExtTable<int64_t, // first read (foreign key)
                              int64_t // second read (foreign key)
                              >
        TP_PAIRED_READ_TABLE;
    class PairedReadTable : public TP_PAIRED_READ_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        std::shared_ptr<ReadTable> pReadTable;

      public:
        PairedReadTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase, std::shared_ptr<ReadTable> pReadTable );

        ~PairedReadTable( )
        {} // deconstructor

        inline std::pair<int64_t, int64_t> insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pReadA,
                                                       std::shared_ptr<NucSeq> pReadB )
        {
            int64_t uiReadIdA = pReadTable->insertRead( uiSequencerId, pReadA );
            int64_t uiReadIdB = pReadTable->insertRead( uiSequencerId, pReadB );
            xInsertRow( uiReadIdA, uiReadIdB );
            return std::make_pair( uiReadIdA, uiReadIdB );
        } // method
    }; // class

    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string, // name
                                                     std::string, // desc
                                                     int64_t // timestamp
                                                     >
        TP_SV_CALLER_RUN_TABLE;
    class SvCallerRunTable : public TP_SV_CALLER_RUN_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        CppSQLiteExtQueryStatement<int64_t> xDeleteRun;
        CppSQLiteExtQueryStatement<int64_t> xGetRunId;
        CppSQLiteExtQueryStatement<std::string, std::string, int64_t> xGetRunName;
        CppSQLiteExtQueryStatement<uint32_t> xNumRuns;
        CppSQLiteExtQueryStatement<uint32_t> xRunExists;
        CppSQLiteExtQueryStatement<uint32_t> xNameExists;

      public:
        SvCallerRunTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SV_CALLER_RUN_TABLE( *pDatabase, // the database where the table resides
                                      "sv_caller_run_table", // name of the table in the database
                                      // column definitions of the table
                                      std::vector<std::string>{"name", "desc", "time_stamp"} ),
              pDatabase( pDatabase ),
              xDeleteRun( *pDatabase, "DELETE FROM sv_caller_run_table WHERE name == ?" ),
              xGetRunId( *pDatabase,
                         "SELECT id FROM sv_caller_run_table WHERE name == ? ORDER BY time_stamp ASC LIMIT 1" ),
              xGetRunName( *pDatabase, "SELECT name, desc, time_stamp FROM sv_caller_run_table WHERE id == ?" ),
              xNumRuns( *pDatabase, "SELECT COUNT(*) FROM sv_caller_run_table" ),
              xRunExists( *pDatabase, "SELECT COUNT(*) FROM sv_caller_run_table WHERE id == ?" ),
              xNameExists( *pDatabase, "SELECT COUNT(*) FROM sv_caller_run_table WHERE name == ?" )
        {} // default constructor

        inline void createIndices( )
        {
            pDatabase->execDML( "CREATE INDEX sv_caller_run_table_id_index ON sv_caller_run_table (id)" );
        } // method

        inline void dropIndices( )
        {
            pDatabase->execDML( "DROP INDEX IF EXISTS sv_caller_run_table_id_index" );
        } // method

        inline void deleteRun( std::string& rS )
        {
            xDeleteRun.bindAndExecQuery<>( rS );
            vDump( std::cout );
        } // method

        inline int64_t getRunId( std::string& rS )
        {
            return xGetRunId.scalar( rS );
        } // method

        inline bool runExists( int64_t iId )
        {
            return xRunExists.scalar( iId ) > 0;
        } // method

        inline bool nameExists( std::string sName )
        {
            return xNameExists.scalar( sName ) > 0;
        } // method

        inline std::string getRunName( int64_t iId )
        {
            return std::get<0>( xGetRunName.vExecuteAndReturnIterator( iId ).get( ) );
        } // method

        inline std::string getRunDesc( int64_t iId )
        {
            return std::get<1>( xGetRunName.vExecuteAndReturnIterator( iId ).get( ) );
        } // method

        inline std::string getRunDate( int64_t iId )
        {
            auto now_c = (std::time_t)std::get<2>( xGetRunName.vExecuteAndReturnIterator( iId ).get( ) );
            std::stringstream ss;
#ifdef _MSC_VER
#pragma warning( suppress : 4996 ) // @todo find another way to do this
            ss << std::put_time( std::localtime( &now_c ), "%c" );
#else
            ss << std::put_time( std::localtime( &now_c ), "%c" );
#endif
            return ss.str( );
        } // method

        inline uint32_t size( )
        {
            return xNumRuns.scalar( );
        } // method

        inline int64_t insertRun( std::string sName, std::string sDesc )
        {
            return this->xInsertRow(
                sName, sDesc, (int64_t)std::chrono::system_clock::to_time_t( std::chrono::system_clock::now( ) ) );
        } // method
    }; // class

    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int64_t, // sv_caller_run_id (foreign key)
                                                     int64_t, // read_id (foreign key)
                                                     int64_t, // sort_pos_start
                                                     int64_t, // sort_pos_end
                                                     uint32_t, // from_pos
                                                     uint32_t, // to_pos
                                                     uint32_t, // query_from
                                                     uint32_t, // query_to
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
            : TP_SV_JUMP_TABLE(
                  *pDatabase, // the database where the table resides
                  "sv_jump_table", // name of the table in the database
                  // column definitions of the table
                  std::vector<std::string>{"sv_caller_run_id", "read_id", "sort_pos_start", "sort_pos_end", "from_pos",
                                           "to_pos", "query_from", "query_to", "from_forward", "to_forward",
                                           "from_seed_start"},
                  // constraints for table
                  std::vector<std::string>{
                      "FOREIGN KEY (sv_caller_run_id) REFERENCES sv_caller_run_table(id) ON DELETE CASCADE",
                      "FOREIGN KEY (read_id) REFERENCES read_table(id)"} ),
              pDatabase( pDatabase ),
              xQuerySize( *pDatabase, "SELECT COUNT(*) FROM sv_jump_table" ),
              xDeleteRun( *pDatabase, "DELETE FROM sv_jump_table WHERE sv_caller_run_id IN ( SELECT id FROM "
                                      "sv_caller_run_table WHERE name == ?)" )
        {} // default constructor

        inline void createIndices( )
        {
            // https://www.sqlite.org/queryplanner.html -> 3.2. Searching And Sorting With A Covering Index
            // index intended for the sweep over the start of all sv-rectangles
            pDatabase->execDML(
                "CREATE INDEX sv_jump_table_sort_index_start ON sv_jump_table"
                "(sv_caller_run_id, sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward,"
                " to_forward, from_seed_start, id, read_id)" );
            // index intended for the sweep over the end of all sv-rectangles
            pDatabase->execDML( "CREATE INDEX sv_jump_table_sort_index_end ON sv_jump_table"
                                "(sv_caller_run_id, sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward,"
                                " to_forward, from_seed_start, id, read_id)" );
        } // method

        inline void dropIndices( )
        {
            pDatabase->execDML( "DROP INDEX IF EXISTS sv_jump_table_sort_index_start" );
            pDatabase->execDML( "DROP INDEX IF EXISTS sv_jump_table_sort_index_end" );
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

    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string, // regex
                                                     uint32_t // state
                                                     >
        TP_SV_CALL_REG_EX_TABLE;
    class SvCallRegExTable : public TP_SV_CALL_REG_EX_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;

      public:
        SvCallRegExTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SV_CALL_REG_EX_TABLE( *pDatabase, // the database where the table resides
                                       "sv_call_reg_ex_table", // name of the table in the database
                                       // column definitions of the table
                                       std::vector<std::string>{"regex", "state"} ),
              pDatabase( pDatabase )
        {} // default constructor
    }; // class

    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int64_t, // sv_caller_run_id (foreign key)
                                                     uint32_t, // from_pos
                                                     uint32_t, // to_pos
                                                     uint32_t, // from_size
                                                     uint32_t, // to_size
                                                     bool, // switch_strand
                                                     NucSeqSql, // inserted_sequence
                                                     double, // score
                                                     int64_t // regex_id
                                                     >
        TP_SV_CALL_TABLE;
    class SvCallTable : public TP_SV_CALL_TABLE
    {

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
            {} // default constructor

        }; // class

        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        std::shared_ptr<ReconstructionTable> pReconstructionTable;
        CppSQLiteExtQueryStatement<uint32_t> xQuerySize;
        CppSQLiteExtQueryStatement<uint32_t> xQuerySizeSpecific;
        CppSQLiteExtQueryStatement<uint32_t> xNumOverlaps;
        CppSQLiteExtQueryStatement<int64_t> xCallArea;
        CppSQLiteExtQueryStatement<int64_t> xDeleteRun;
        CppSQLiteExtQueryStatement<int64_t, bool, uint32_t, uint32_t, NucSeqSql, uint32_t> xNextCallForwardContext;
        CppSQLiteExtQueryStatement<int64_t, bool, uint32_t, uint32_t, NucSeqSql, uint32_t> xNextCallBackwardContext;

      public:
        SvCallTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SV_CALL_TABLE(
                  *pDatabase, // the database where the table resides
                  "sv_call_table", // name of the table in the database
                  // column definitions of the table
                  std::vector<std::string>{"sv_caller_run_id", "from_pos", "to_pos", "from_size", "to_size",
                                           "switch_strand", "inserted_sequence", "score", "regex_id"},
                  // constraints for table
                  std::vector<std::string>{
                      "FOREIGN KEY (sv_caller_run_id) REFERENCES sv_caller_run_table(id) ON DELETE CASCADE",
                      "FOREIGN KEY (regex_id) REFERENCES sv_call_reg_ex_table(id) ON DELETE SET NULL"} ),
              pDatabase( pDatabase ),
              pReconstructionTable( std::make_shared<ReconstructionTable>( pDatabase ) ),
              xQuerySize( *pDatabase, "SELECT COUNT(*) FROM sv_call_table" ),
              xQuerySizeSpecific( *pDatabase, "SELECT COUNT(*) FROM sv_call_table WHERE sv_caller_run_id == ? "
                                              "AND score >= ?" ),
              xNumOverlaps( *pDatabase,
                            "SELECT COUNT(*) "
                            "FROM sv_call_table AS outer "
                            "WHERE sv_caller_run_id == ? "
                            "AND score >= ? "
                            "AND EXISTS ( "
                            "   SELECT * "
                            "   FROM sv_call_table AS inner "
                            "   WHERE inner.sv_caller_run_id == ? "
                            "   AND score >= ?"
                            "   AND outer.from_pos <= inner.from_pos + inner.from_size + ? "
                            "   AND outer.from_pos + outer.from_size + ? >= inner.from_pos "
                            "   AND outer.to_pos <= inner.to_pos + inner.to_size + ? "
                            "   AND outer.to_pos + outer.to_size + ? >= inner.to_pos "
                            "   AND outer.switch_strand == inner.switch_strand "
                            ")" ),
              xCallArea( *pDatabase,
                         "SELECT SUM( from_size * to_size ) FROM sv_call_table WHERE sv_caller_run_id == ? "
                         "AND score >= ?" ),
              xDeleteRun( *pDatabase, "DELETE FROM sv_call_table WHERE sv_caller_run_id IN ( SELECT id FROM "
                                      "sv_caller_run_table WHERE name == ?)" ),
              xNextCallForwardContext(
                  *pDatabase,
                  "SELECT id, switch_strand, to_pos, to_size, inserted_sequence, from_pos + from_size "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id = ? "
                  "AND from_pos >= ? "
                  "AND id NOT IN ( "
                  "  SELECT call_id "
                  "  FROM reconstruction_table "
                  ") "
                  "ORDER BY from_pos ASC "
                  "LIMIT 1 " ),
              xNextCallBackwardContext( *pDatabase,
                                        "SELECT id, switch_strand, from_pos, from_size, inserted_sequence, to_pos "
                                        "FROM sv_call_table "
                                        "WHERE sv_caller_run_id = ? "
                                        "AND to_pos <= ? "
                                        "AND id NOT IN ( "
                                        "  SELECT call_id "
                                        "  FROM reconstruction_table "
                                        ") "
                                        "ORDER BY to_pos DESC "
                                        "LIMIT 1 " )
        {} // default constructor

        inline uint32_t numCalls( )
        {
            return xQuerySize.scalar( );
        } // method

        inline uint32_t numCalls( int64_t iCallerRunId, double dMinScore )
        {
            return xQuerySizeSpecific.scalar( iCallerRunId, dMinScore );
        } // method

        inline int64_t insertCall( int64_t iSvCallerRunId, SvCall& rCall )
        {
            int64_t iCallId =
                this->xInsertRow( iSvCallerRunId, (uint32_t)rCall.uiFromStart, (uint32_t)rCall.uiToStart,
                                  (uint32_t)rCall.uiFromSize, (uint32_t)rCall.uiToSize, rCall.bSwitchStrand,
                                  // NucSeqSql can deal with nullpointers
                                  NucSeqSql( rCall.pInsertedSequence ), rCall.dScore, -1 );
            rCall.iId = iCallId;
            return iCallId;
        } // method

        inline void deleteRun( std::string& rS )
        {
            xDeleteRun.bindAndExecQuery<>( rS );
        } // method

        inline int64_t callArea( int64_t iCallerRunId, double dMinScore )
        {
            return xCallArea.scalar( iCallerRunId, dMinScore );
        } // method

        /**
         * returns how many calls of run A are overlapped by a call in run B
         */
        inline uint32_t numOverlaps( int64_t iCallerRunIdA, int64_t iCallerRunIdB, double dMinScore,
                                     int64_t iAllowedDist )
        {
            return xNumOverlaps.scalar( iCallerRunIdA, dMinScore, iCallerRunIdB, dMinScore, iAllowedDist, iAllowedDist,
                                        iAllowedDist, iAllowedDist );
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
                    std::get<4>( xRet ) =
                        !std::get<1>( xQ ) ? std::get<2>( xQ ) + std::get<3>( xQ ) : std::get<2>( xQ );
                } // if
            } // else
            return xRet;
        } // method

        inline std::shared_ptr<Pack> reconstructSequencedGenome( std::shared_ptr<Pack> pRef, int64_t iCallerRun )
        {
            // @todo at the moment this does not deal with jumped over sequences
            auto pRet = std::make_shared<Pack>( );

            NucSeq xCurrChrom;
            uint32_t uiCurrPos = 0;
            bool bForwContext = true;
            while( true )
            {
                // get the next call
                auto tNextCall = this->getNextCall( iCallerRun, uiCurrPos, bForwContext );
#if 0
                std::cout << "id: " << std::get<0>( tNextCall ) << " from: " << std::get<1>( tNextCall )
                          << " to: " << std::get<4>( tNextCall )
                          << ( std::get<2>( tNextCall ) ? " forward" : " rev-comp" ) << std::endl;
#endif
                if( std::get<0>( tNextCall ) == -1 ) // if this was the last call
                {
                    pRef->vExtractContext( uiCurrPos, xCurrChrom, true, bForwContext );
                    pRet->vAppendSequence( "unnamed_contig", "no_description_given", xCurrChrom );
                    xCurrChrom.vClear( );
                    /*
                     * for this we make use of the id system of contigs.
                     * the n forwards contigs have the ids: x*2 | 0 <= x <= n
                     * the n reverse complement contigs have the ids: x*2+1 | 0 <= x <= n
                     */
                    for( int64_t uiI = pRef->uiSequenceIdForPositionOrRev( uiCurrPos ) + ( bForwContext ? 2 : -1 );
                         uiI < (int64_t)pRef->uiNumContigs( ) * 2 && uiI >= 0;
                         uiI += bForwContext ? 2 : -2 )
                    {
                        pRef->vExtractContig( uiI, xCurrChrom, true );
                        pRet->vAppendSequence( "unnamed_contig", "no_description_given", xCurrChrom );
                        xCurrChrom.vClear( );
                    } // for
                    break;
                } // if

                if( pRef->bridgingPositions( uiCurrPos, std::get<1>( tNextCall ) ) )
                {
                    pRef->vExtractContext( uiCurrPos, xCurrChrom, true, bForwContext );
                    pRet->vAppendSequence( "unnamed_contig", "no_description_given", xCurrChrom );
                    xCurrChrom.vClear( );
                    pRef->vExtractContext( std::get<1>( tNextCall ), xCurrChrom, true, !bForwContext );
                } // if
                else if( bForwContext )
                    pRef->vExtractSubsectionN( uiCurrPos, std::get<1>( tNextCall ), xCurrChrom, true );
                else
                    pRef->vExtractSubsectionN( pRef->uiPositionToReverseStrand( uiCurrPos ) + 1,
                                               pRef->uiPositionToReverseStrand( std::get<1>( tNextCall ) ) + 1, //
                                               xCurrChrom,
                                               true );
                if( std::get<3>( tNextCall ).pNucSeq != nullptr )
                    xCurrChrom.vAppend( std::get<3>( tNextCall ).pNucSeq->pxSequenceRef,
                                        std::get<3>( tNextCall ).pNucSeq->length( ) );

                // remember that we used this call
                pReconstructionTable->xInsertRow( std::get<0>( tNextCall ) );
                uiCurrPos = std::get<4>( tNextCall );
                bForwContext = std::get<2>( tNextCall );
            } // while

            // clear up the temp table
            pReconstructionTable->clearTable( );

            return pRet;
        } // method
    }; // class

    typedef CppSQLiteExtTable<int64_t, // call_id (foreign key)
                              int64_t // jump_id (foreign key)
                              >
        TP_SV_CALL_SUPPORT_TABLE;
    class SvCallSupportTable : public TP_SV_CALL_SUPPORT_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        CppSQLiteExtQueryStatement<int64_t> xDeleteRun;

      public:
        SvCallSupportTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SV_CALL_SUPPORT_TABLE(
                  *pDatabase, // the database where the table resides
                  "sv_call_support_table", // name of the table in the database
                  // column definitions of the table
                  std::vector<std::string>{"call_id", "jump_id"},
                  false,
                  // constraints for table
                  std::vector<std::string>{"FOREIGN KEY (call_id) REFERENCES sv_call_table(id) ON DELETE CASCADE",
                                           "FOREIGN KEY (jump_id) REFERENCES sv_jump_table(id) ON DELETE CASCADE"} ),
              pDatabase( pDatabase ),
              xDeleteRun( *pDatabase, "DELETE FROM sv_call_support_table WHERE call_id IN ( SELECT id FROM "
                                      "sv_call_table WHERE sv_caller_run_id IN ( SELECT id FROM "
                                      "sv_caller_run_table WHERE name == ?))" )
        {} // default constructor

        inline void deleteRun( std::string& rS )
        {
            xDeleteRun.bindAndExecQuery<>( rS );
        } // method
    }; // class


    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    std::shared_ptr<SequencerTable> pSequencerTable;
    std::shared_ptr<ReadTable> pReadTable;
    std::shared_ptr<PairedReadTable> pPairedReadTable;
    std::shared_ptr<SvCallerRunTable> pSvCallerRunTable;
    std::shared_ptr<SvJumpTable> pSvJumpTable;
    std::shared_ptr<SvCallRegExTable> pSvCallRegExTable;
    std::shared_ptr<SvCallTable> pSvCallTable;
    std::shared_ptr<SvCallSupportTable> pSvCallSupportTable;

    friend class NucSeqFromSql;
    friend class AllNucSeqFromSql;
    friend class PairedNucSeqFromSql;
    friend class PairedReadTable;
    friend class SortedSvJumpFromSql;
    friend class SvCallerRunsFromDb;
    friend class SvCallsFromDb;

  public:
    SV_DB( std::string sName, enumSQLite3DBOpenMode xMode )
        : pDatabase( std::make_shared<CppSQLiteDBExtended>( "", sName, xMode ) ),
          pSequencerTable( std::make_shared<SequencerTable>( pDatabase ) ),
          pReadTable( std::make_shared<ReadTable>( pDatabase ) ),
          pPairedReadTable( std::make_shared<PairedReadTable>( pDatabase, pReadTable ) ),
          pSvCallerRunTable( std::make_shared<SvCallerRunTable>( pDatabase ) ),
          pSvJumpTable( std::make_shared<SvJumpTable>( pDatabase ) ),
          pSvCallRegExTable( std::make_shared<SvCallRegExTable>( pDatabase ) ),
          pSvCallTable( std::make_shared<SvCallTable>( pDatabase ) ),
          pSvCallSupportTable( std::make_shared<SvCallSupportTable>( pDatabase ) )
    {} // constructor

    SV_DB( std::string sName ) : SV_DB( sName, eCREATE_DB )
    {} // constructor

    SV_DB( std::string sName, std::string sMode ) : SV_DB( sName, sMode == "create" ? eCREATE_DB : eOPEN_DB )
    {} // constructor

    inline void createSequencerIndices( )
    {
        pSequencerTable->createIndices( );
        pReadTable->createIndices( );
    } // method

    inline void createCallerIndices( )
    {
        pSvCallerRunTable->createIndices( );
        pSvJumpTable->createIndices( );
    } // method

    inline void dropCallerIndices( )
    {
        pSvCallerRunTable->dropIndices( );
        pSvJumpTable->dropIndices( );
    } // method

    inline void setNumThreads( size_t uiN )
    {
        pDatabase->set_num_threads( (int)uiN );
    } // method

    inline int64_t getRunId( std::string& rS )
    {
        return pSvCallerRunTable->getRunId( rS );
    } // method

    inline int64_t getCallArea( int64_t iCallerRunId, double dMinScore )
    {
        return pSvCallTable->callArea( iCallerRunId, dMinScore );
    } // method

    inline uint32_t getNumOverlapsBetweenCalls( int64_t iCallerRunIdA, int64_t iCallerRunIdB, double dMinScore,
                                                int64_t iAllowedDist )
    {
        return pSvCallTable->numOverlaps( iCallerRunIdA, iCallerRunIdB, dMinScore, iAllowedDist );
    } // method

    inline uint32_t getNumCalls( int64_t iCallerRunId, double dMinScore )
    {
        return pSvCallTable->numCalls( iCallerRunId, dMinScore );
    } // method

    inline uint32_t getNumRuns( )
    {
        return pSvCallerRunTable->size( );
    } // method

    inline std::string getRunName( int64_t iId )
    {
        return pSvCallerRunTable->getRunName( iId );
    } // method

    inline std::string getRunDesc( int64_t iId )
    {
        return pSvCallerRunTable->getRunDesc( iId );
    } // method

    inline std::string getRunDate( int64_t iId )
    {
        return pSvCallerRunTable->getRunDate( iId );
    } // method

    inline bool runExists( int64_t iId )
    {
        return pSvCallerRunTable->runExists( iId );
    } // method

    inline bool nameExists( std::string sName )
    {
        return pSvCallerRunTable->nameExists( sName );
    } // method

    inline std::shared_ptr<Pack> reconstructSequencedGenome( std::shared_ptr<Pack> pRef, int64_t iCallerRun )
    {
        return pSvCallTable->reconstructSequencedGenome( pRef, iCallerRun );
    } // method

    class ReadInserter
    {
      private:
        // this is here so that it gets destructed after the transaction context
        std::shared_ptr<SV_DB> pDB;
        // must be after the DB so that it is deconstructed first
        CppSQLiteExtImmediateTransactionContext xTransactionContext;

        int64_t uiSequencerId;

      public:
        ReadInserter( std::shared_ptr<SV_DB> pDB, std::string sSequencerName )
            : pDB( pDB ),
              xTransactionContext( *pDB->pDatabase ),
              uiSequencerId( pDB->pSequencerTable->insertSequencer( sSequencerName ) )
        {} // constructor

        ReadInserter( const ReadInserter& rOther ) = delete; // delete copy constructor

        inline void insertRead( std::shared_ptr<NucSeq> pRead )
        {
            pDB->pReadTable->insertRead( uiSequencerId, pRead );
        } // method

        inline void insertPairedRead( std::shared_ptr<NucSeq> pReadA, std::shared_ptr<NucSeq> pReadB )
        {
            pDB->pPairedReadTable->insertRead( uiSequencerId, pReadA, pReadB );
        } // method
    }; // class

    class SvJumpInserter
    {
        // this is here so that it gets destructed after the transaction context
        std::shared_ptr<SV_DB> pDB;
        // must be after the DB so that it is deconstructed first
        CppSQLiteExtImmediateTransactionContext xTransactionContext;

      public:
        const int64_t iSvCallerRunId;

        class ReadContex
        {
          private:
            std::shared_ptr<SvJumpTable> pSvJumpTable;
            const int64_t iSvCallerRunId;
            const int64_t iReadId;

          public:
            ReadContex( std::shared_ptr<SvJumpTable> pSvJumpTable, const int64_t iSvCallerRunId, const int64_t iReadId )
                : pSvJumpTable( pSvJumpTable ), iSvCallerRunId( iSvCallerRunId ), iReadId( iReadId )
            {} // constructor

            inline void insertJump( SvJump& rJump )
            {
                // make shure the read id mathes the read context
                if( rJump.iReadId == -1 ) // if there is no read id given yet add it
                    rJump.iReadId = iReadId;
                else // otherwise assert it matches
                    assert( rJump.iReadId == iReadId );

                if( rJump.does_switch_strand( ) )
                    assert( rJump.from_start( ) > std::numeric_limits<int64_t>::max( ) / 2 );
                rJump.iId = pSvJumpTable->xInsertRow( iSvCallerRunId, rJump.iReadId, rJump.from_start( ),
                                                      rJump.from_end( ), (uint32_t)rJump.uiFrom, (uint32_t)rJump.uiTo,
                                                      (uint32_t)rJump.uiQueryFrom, (uint32_t)rJump.uiQueryTo,
                                                      rJump.bFromForward, rJump.bToForward, rJump.bFromSeedStart );
            } // method
        }; // class

        SvJumpInserter( std::shared_ptr<SV_DB> pDB,
                        const std::string& rsSvCallerName,
                        const std::string& rsSvCallerDesc )
            : pDB( pDB ),
              xTransactionContext( *pDB->pDatabase ),
              iSvCallerRunId( pDB->pSvCallerRunTable->insertRun( rsSvCallerName, rsSvCallerDesc ) )
        {} // constructor

        SvJumpInserter( const SvJumpInserter& ) = delete; // delete copy constructor

        inline ReadContex insertRead( std::shared_ptr<NucSeq> pRead )
        {
            return ReadContex( pDB->pSvJumpTable, iSvCallerRunId,
                               pDB->pReadTable->insertRead( iSvCallerRunId, pRead ) );
        } // method

        inline ReadContex readContext( int64_t iReadId )
        {
            return ReadContex( pDB->pSvJumpTable, iSvCallerRunId, iReadId );
        } // method

    }; // class

    class SvCallInserter
    {
        // this is here so that it gets destructed after the transaction context
        std::shared_ptr<SV_DB> pDB;
        // must be after the DB so that it is deconstructed first
        CppSQLiteExtImmediateTransactionContext xTransactionContext;

      public:
        const int64_t iSvCallerRunId;

        class CallContex
        {
          private:
            std::shared_ptr<SvCallSupportTable> pSvCallSupportTable;
            const int64_t iCallId;

          public:
            CallContex( std::shared_ptr<SvCallSupportTable> pSvCallSupportTable, const int64_t iCallId )
                : pSvCallSupportTable( pSvCallSupportTable ), iCallId( iCallId )
            {} // constructor

            inline void addSupport( SvJump& rJump )
            {
                pSvCallSupportTable->xInsertRow( iCallId, rJump.iId );
            } // method

            inline void addSupport( int64_t iId )
            {
                pSvCallSupportTable->xInsertRow( iCallId, iId );
            } // method
        }; // class

        SvCallInserter( std::shared_ptr<SV_DB> pDB, const int64_t iSvCallerRunId )
            : pDB( pDB ), xTransactionContext( *pDB->pDatabase ), iSvCallerRunId( iSvCallerRunId )
        {} // constructor

        SvCallInserter( const SvCallInserter& ) = delete; // delete copy constructor

        inline void insertCall( SvCall& rCall )
        {
            CallContex xContext( pDB->pSvCallSupportTable, pDB->pSvCallTable->insertCall( iSvCallerRunId, rCall ) );
            for( int64_t iId : rCall.vSupportingJumpIds )
                xContext.addSupport( iId );
        } // method

    }; // class

    inline void clearCallsTable( )
    {
        pSvCallSupportTable->clearTable( );
        pSvCallTable->clearTable( );
        pSvJumpTable->clearTable( );
        pSvCallerRunTable->clearTable( );
    } // method

    inline void clearCallsTableForCaller( std::string& rS )
    {
        pSvCallSupportTable->deleteRun( rS );
        pSvCallTable->deleteRun( rS );
        pSvJumpTable->deleteRun( rS );
        pSvCallerRunTable->deleteRun( rS );
    } // method

    inline uint32_t numJumps( )
    {
        return pSvJumpTable->numJumps( );
    } // method

    inline uint32_t numCalls( )
    {
        return pSvCallTable->numCalls( );
    } // method

}; // class

class SortedSvJumpFromSql
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t, int64_t, int64_t>
        xQueryStart;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t, int64_t, int64_t>
        xQueryEnd;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t, int64_t,
                               int64_t>::Iterator xTableIteratorStart;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t, int64_t,
                               int64_t>::Iterator xTableIteratorEnd;

  public:
    SortedSvJumpFromSql( std::shared_ptr<SV_DB> pDb, int64_t iSvCallerRunId )
        : pDb( pDb ),
          xQueryStart( *pDb->pDatabase,
                       "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                       "sort_pos_start, id, read_id "
                       "FROM sv_jump_table "
                       "WHERE sv_caller_run_id == ? "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd( *pDb->pDatabase,
                     "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                     "sort_pos_end, id, read_id "
                     "FROM sv_jump_table "
                     "WHERE sv_caller_run_id == ? "
                     "ORDER BY sort_pos_end" ),
          xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator( iSvCallerRunId ) ),
          xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId ) )
    {} // constructor

    bool hasNextStart( )
    {
        return !xTableIteratorStart.eof( );
    } // method

    bool hasNextEnd( )
    {
        return !xTableIteratorEnd.eof( );
    } // method

    bool nextStartIsSmaller( )
    {
        auto xStartTup = xTableIteratorStart.get( );
        auto xEndTup = xTableIteratorEnd.get( );
        return std::get<7>( xStartTup ) <= std::get<7>( xEndTup );
    } // method

    SvJump getNextStart( )
    {
        assert( hasNextStart( ) );

        auto xTup = xTableIteratorStart.get( );
        xTableIteratorStart.next( );
        return SvJump( std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                       std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ), std::get<8>( xTup ),
                       std::get<9>( xTup ) );
    } // method

    SvJump getNextEnd( )
    {
        assert( hasNextEnd( ) );

        auto xTup = xTableIteratorEnd.get( );
        xTableIteratorEnd.next( );
        return SvJump( std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                       std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ), std::get<8>( xTup ),
                       std::get<9>( xTup ) );
    } // method

}; // class

class AllNucSeqFromSql : public Module<NucSeq, true>
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t> xQuery;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t>::Iterator xTableIterator;

  public:
    AllNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT read_table.sequence, read_table.id "
                  "FROM read_table " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( ) )
    {
        if( xTableIterator.eof( ) )
            setFinished( );
    } // constructor

    std::shared_ptr<NucSeq> execute( )
    {
        if( xTableIterator.eof( ) )
            throw AnnotatedException( "No more NucSeq in NucSeqFromSql module" );

        auto xTup = xTableIterator.get( );
        // std::get<0>( xTup ).pNucSeq->sName = std::to_string( std::get<1>( xTup ) );
        std::get<0>( xTup ).pNucSeq->iId = std::get<1>( xTup );
        xTableIterator.next( );

        if( xTableIterator.eof( ) )
            setFinished( );
        return std::get<0>( xTup ).pNucSeq;
    } // method

    // override
    bool requiresLock( ) const
    {
        return true;
    } // method
}; // class

class NucSeqFromSql : public Module<NucSeq, true>
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t> xQuery;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t>::Iterator xTableIterator;

  public:
    NucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT read_table.sequence, read_table.id "
                  "FROM read_table "
                  "WHERE read_table.id NOT IN ( "
                  "   SELECT paired_read_table.first_read FROM paired_read_table "
                  "   UNION "
                  "   SELECT paired_read_table.second_read FROM paired_read_table "
                  ") " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( ) )
    {
        if( xTableIterator.eof( ) )
            setFinished( );
    } // constructor

    std::shared_ptr<NucSeq> execute( )
    {
        if( xTableIterator.eof( ) )
            throw AnnotatedException( "No more NucSeq in NucSeqFromSql module" );

        auto xTup = xTableIterator.get( );
        // std::get<0>( xTup ).pNucSeq->sName = std::to_string( std::get<1>( xTup ) );
        std::get<0>( xTup ).pNucSeq->iId = std::get<1>( xTup );
        xTableIterator.next( );

        if( xTableIterator.eof( ) )
            setFinished( );
        return std::get<0>( xTup ).pNucSeq;
    } // method

    // override
    bool requiresLock( ) const
    {
        return true;
    } // method
}; // class

class PairedNucSeqFromSql : public Module<ContainerVector<std::shared_ptr<NucSeq>>, true>
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, NucSeqSql, uint32_t, uint32_t> xQuery;
    CppSQLiteExtQueryStatement<NucSeqSql, NucSeqSql, uint32_t, uint32_t>::Iterator xTableIterator;
    const bool bRevCompMate;

  public:
    PairedNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT A.sequence, B.sequence, A.id, B.id "
                  "FROM read_table A, read_table B "
                  "INNER JOIN paired_read_table "
                  "ON paired_read_table.first_read == A.id "
                  "AND paired_read_table.second_read == B.id " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( ) ),
          bRevCompMate( rParameters.getSelected( )->xRevCompPairedReadMates->get( ) )
    {
        if( xTableIterator.eof( ) )
            setFinished( );
    } // constructor

    std::shared_ptr<ContainerVector<std::shared_ptr<NucSeq>>> execute( )
    {
        if( xTableIterator.eof( ) )
            throw AnnotatedException( "No more NucSeq in PairedNucSeqFromSql module" );

        auto pRet = std::make_shared<ContainerVector<std::shared_ptr<NucSeq>>>( );

        auto xTup = xTableIterator.get( );
        pRet->push_back( std::get<0>( xTup ).pNucSeq );
        pRet->back( )->iId = std::get<2>( xTup );
        pRet->push_back( std::get<1>( xTup ).pNucSeq );
        pRet->back( )->iId = std::get<3>( xTup );

        if( bRevCompMate )
        {
            pRet->back( )->vReverse( );
            pRet->back( )->vSwitchAllBasePairsToComplement( );
        } // if

        xTableIterator.next( );

        if( xTableIterator.eof( ) )
            setFinished( );
        return pRet;
    } // method

    // override
    bool requiresLock( ) const
    {
        return true;
    } // method
}; // class

class SvDbInserter : public Module<Container, false, ContainerVector<SvJump>, NucSeq>
{
    std::shared_ptr<SV_DB> pDb;
    std::mutex xMutex;

  public:
    SV_DB::SvJumpInserter xInserter;

    SvDbInserter( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, std::string sRunDesc )
        : pDb( pDb ), xInserter( pDb, "MA-SV", sRunDesc )
    {} // constructor

    std::shared_ptr<Container> execute( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        std::lock_guard<std::mutex> xGuard( xMutex );

        SV_DB::SvJumpInserter::ReadContex xReadContext = xInserter.readContext( pRead->iId );
        for( SvJump& rJump : *pJumps )
            xReadContext.insertJump( rJump ); // also updates the jump ids;
        return std::make_shared<Container>( );
        // end of score for xGuard
    } // method
}; // class


class SvCallerRunsFromDb
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<int64_t, std::string, std::string> xQuery;
    CppSQLiteExtQueryStatement<int64_t, std::string, std::string>::Iterator xTableIterator;

  public:
    SvCallerRunsFromDb( std::shared_ptr<SV_DB> pDb )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, name, desc "
                  "FROM sv_caller_run_table " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( ) )
    {} // constructor

    int64_t id( )
    {
        return std::get<0>( xTableIterator.get( ) );
    } // method

    std::string name( )
    {
        return std::get<1>( xTableIterator.get( ) );
    } // method

    std::string desc( )
    {
        return std::get<2>( xTableIterator.get( ) );
    } // method

    void next( )
    {
        xTableIterator.next( );
    } // method

    bool eof( )
    {
        return xTableIterator.eof( );
    } // method
}; // class

class SvCallsFromDb
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, double> xQuery;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t> xQuerySupport;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, double>::Iterator
        xTableIterator;

  public:
    SvCallsFromDb( std::shared_ptr<SV_DB> pDb, int64_t iSvCallerId )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, score "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id == ? " ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "sv_jump_table.id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerId ) )
    {} // constructor

    SvCall next( )
    {
        auto xTup = xTableIterator.get( );
        SvCall xRet( std::get<1>( xTup ), // uiFromStart
                     std::get<2>( xTup ), // uiToStart
                     std::get<3>( xTup ), // uiFromSize
                     std::get<4>( xTup ), // uiToSize
                     std::get<5>( xTup ), // bSwitchStrand
                     std::get<7>( xTup ) // dScore
        );
        xRet.pInsertedSequence = std::get<6>( xTup ).pNucSeq;
        xRet.iId = std::get<0>( xTup );
        auto xSupportIterator( xQuerySupport.vExecuteAndReturnIterator( std::get<0>( xTup ) ) );
        while( !xSupportIterator.eof( ) )
        {
            auto xTup = xSupportIterator.get( );
            xRet.vSupportingJumpIds.push_back( std::get<7>( xTup ) );
            xRet.vSupportingJumps.emplace_back( std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ),
                                                std::get<3>( xTup ), std::get<4>( xTup ), std::get<5>( xTup ),
                                                std::get<6>( xTup ), std::get<7>( xTup ) );
            xSupportIterator.next( );
        } // while
        xTableIterator.next( );
        return xRet;
    } // method

    bool hasNext( )
    {
        return !xTableIterator.eof( );
    } // method
}; // class

}; // namespace libMA

#ifdef WITH_PYTHON
#ifdef WITH_BOOST
void exportSoCDbWriter( );
#else
void exportSoCDbWriter( py::module& rxPyModuleId );
#endif
#endif
