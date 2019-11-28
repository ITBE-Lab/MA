/**
 * @file svDb.h
 * @details
 * The database interface for the structural variant caller
 */
#pragma once

#include "container/container.h"
#include "container/nucSeq.h"
#include "container/pack.h"
#include "container/soc.h"
#include "container/svJump.h"
#include "module/fileReader.h"
#include "module/module.h"
#include "util/exception.h"
#include "util/sqlite3.h"
#include "util/system.h"
#include <chrono>
#include <ctime>
#include <iomanip>
// include all table definitions
#include "container/sv_db/tables/contigCoverage.h"
#include "container/sv_db/tables/nameDesc.h"
#include "container/sv_db/tables/pairedRead.h"
#include "container/sv_db/tables/read.h"
#include "container/sv_db/tables/sequencer.h"
#include "container/sv_db/tables/svCall.h"
#include "container/sv_db/tables/svCallRegEx.h"
#include "container/sv_db/tables/svCallSupport.h"
#include "container/sv_db/tables/svCallerRun.h"
#include "container/sv_db/tables/svJump.h"

namespace libMA
{

class SV_DB : public Container
{
  public:
    const std::string sName;
    std::shared_ptr<std::mutex> pWriteLock;
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    std::shared_ptr<SequencerTable> pSequencerTable;
    std::shared_ptr<ContigCovTable> pContigCovTable;
    std::shared_ptr<ReadTable> pReadTable;
    std::shared_ptr<PairedReadTable> pPairedReadTable;
    std::shared_ptr<NameDescTable> pSvJumpRunTable;
    std::shared_ptr<SvJumpTable> pSvJumpTable;
    std::shared_ptr<SvCallerRunTable> pSvCallerRunTable;
    std::shared_ptr<SvCallRegExTable> pSvCallRegExTable;
    std::shared_ptr<SvCallTable> pSvCallTable;
    std::shared_ptr<SvCallSupportTable> pSvCallSupportTable;

    /**
     * @brief open a new database connection; with shared metadata (table pointers are shared)
     */
    SV_DB( SV_DB& rOther )
        : sName( rOther.sName ),
          pWriteLock( rOther.pWriteLock ),
          pDatabase( std::make_shared<CppSQLiteDBExtended>( "", rOther.sName, eOPEN_DB ) ),
          pSequencerTable( rOther.pSequencerTable ),
          pContigCovTable( rOther.pContigCovTable ),
          pReadTable( rOther.pReadTable ),
          pPairedReadTable( rOther.pPairedReadTable ),
          pSvJumpRunTable( rOther.pSvJumpRunTable ),
          pSvJumpTable( rOther.pSvJumpTable ),
          pSvCallerRunTable( rOther.pSvCallerRunTable ),
          pSvCallRegExTable( rOther.pSvCallRegExTable ),
          pSvCallTable( rOther.pSvCallTable ),
          pSvCallSupportTable( rOther.pSvCallSupportTable )
    {
        this->setNumThreads( 32 ); // @todo do this via a parameter
        pDatabase->execDML( "PRAGMA journal_mode=WAL;" ); // use write ahead mode
        pDatabase->execDML( "PRAGMA busy_timeout=0;" ); // do not throw sqlite busy errors
        // https://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
        // pDatabase->execDML( "PRAGMA synchronous = OFF;" ); // insert performance
        // pDatabase->execDML( "PRAGMA journal_mode = MEMORY;" ); // insert performance
    } // constructor

    SV_DB( std::string sName, enumSQLite3DBOpenMode xMode )
        : sName( sName ),
          pWriteLock( std::make_shared<std::mutex>( ) ),
          pDatabase( std::make_shared<CppSQLiteDBExtended>( "", sName, xMode ) ),
          pSequencerTable( std::make_shared<SequencerTable>( pDatabase ) ),
          pContigCovTable( std::make_shared<ContigCovTable>( pDatabase ) ),
          pReadTable( std::make_shared<ReadTable>( pDatabase ) ),
          pPairedReadTable( std::make_shared<PairedReadTable>( pDatabase, pReadTable ) ),
          pSvJumpRunTable( std::make_shared<NameDescTable>( pDatabase, "sv_jump_run_table" ) ),
          pSvJumpTable( std::make_shared<SvJumpTable>( pDatabase ) ),
          pSvCallerRunTable( std::make_shared<SvCallerRunTable>( pDatabase ) ),
          pSvCallRegExTable( std::make_shared<SvCallRegExTable>( pDatabase ) ),
          pSvCallTable( std::make_shared<SvCallTable>( pDatabase ) ),
          pSvCallSupportTable( std::make_shared<SvCallSupportTable>( pDatabase ) )
    {
        this->setNumThreads( 32 ); // @todo do this via a parameter
        pDatabase->execDML( "PRAGMA journal_mode=WAL;" ); // use write ahead mode
        pDatabase->execDML( "PRAGMA busy_timeout=0;" ); // do not throw sqlite busy errors
        if( xMode == eCREATE_DB )
        {
            // https://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
            pDatabase->execDML( "PRAGMA synchronous = OFF;" ); // insert performance
            pDatabase->execDML( "PRAGMA journal_mode = MEMORY;" ); // insert performance
        }
    } // constructor

    SV_DB( std::string sName ) : SV_DB( sName, eCREATE_DB )
    {} // constructor

    SV_DB( std::string sName, std::string sMode ) : SV_DB( sName, sMode == "create" ? eCREATE_DB : eOPEN_DB )
    {} // constructor

    inline void createJumpIndices( int64_t uiRun )
    {
        pSvJumpTable->createIndices( uiRun );
    } // method


    inline void addScoreIndex( int64_t iCallerRunId )
    {
        pSvCallTable->addScoreIndex( iCallerRunId );
    } // method

    inline void setNumThreads( size_t uiN )
    {
        pDatabase->set_num_threads( (int)uiN );
    } // method

    inline int64_t getRunId( std::string& rS )
    {
        return pSvCallerRunTable->getId( rS );
    } // method

    inline int64_t getCallArea( int64_t iCallerRunId, double dMinScore )
    {
        return pSvCallTable->callArea( iCallerRunId, dMinScore );
    } // method

    inline double getMaxScore( int64_t iCallerRunId )
    {
        return pSvCallTable->maxScore( iCallerRunId );
    } // method

    inline double getMinScore( int64_t iCallerRunId )
    {
        return pSvCallTable->minScore( iCallerRunId );
    } // method

    inline uint32_t getNumOverlapsBetweenCalls( int64_t iCallerRunIdA, int64_t iCallerRunIdB, double dMinScore,
                                                int64_t iAllowedDist )
    {
        return pSvCallTable->numOverlaps( iCallerRunIdA, iCallerRunIdB, dMinScore, iAllowedDist );
    } // method

    inline double getBlurOnOverlapsBetweenCalls( int64_t iCallerRunIdA, int64_t iCallerRunIdB, double dMinScore,
                                                 int64_t iAllowedDist )
    {
        return pSvCallTable->blurOnOverlaps( iCallerRunIdA, iCallerRunIdB, dMinScore, iAllowedDist );
    } // method

    inline uint32_t getNumInvalidCalls( int64_t iCallerRunIdA, double dMinScore, int64_t iAllowedDist )
    {
        return pSvCallTable->numInvalidCalls( iCallerRunIdA, dMinScore, iAllowedDist );
    } // method

    inline uint32_t getNumCalls( int64_t iCallerRunId, double dMinScore )
    {
        return pSvCallTable->numCalls( iCallerRunId, dMinScore );
    } // method

    inline uint32_t getNumRuns( )
    {
        return pSvCallerRunTable->size( );
    } // method

    inline std::vector<int64_t> getNumNts( int64_t iSequencerId )
    {
        return pContigCovTable->getNumNt( iSequencerId );
    } // method

    inline std::string getRunName( int64_t iId )
    {
        return pSvCallerRunTable->getName( iId );
    } // method

    inline std::string getRunDesc( int64_t iId )
    {
        return pSvCallerRunTable->getDesc( iId );
    } // method

    inline int64_t getRunJumpId( int64_t iId )
    {
        return pSvCallerRunTable->getSvJumpRunId( iId );
    } // method

    inline std::string getRunDate( int64_t iId )
    {
        return pSvCallerRunTable->getDate( iId );
    } // method

    inline bool runExists( int64_t iId )
    {
        return pSvCallerRunTable->exists( iId );
    } // method

    inline std::vector<int64_t> getNewestUniqueRuns( uint32_t uiNum, std::string sDesc )
    {
        return pSvCallerRunTable->getNewestUnique( uiNum, sDesc );
    } // method

    inline bool nameExists( std::string sName )
    {
        return pSvCallerRunTable->nameExists( sName );
    } // method

    inline int64_t insertSvCallerRun( std::string rsSvCallerName, std::string rsSvCallerDesc, int64_t uiJumpRunId )
    {
        return pSvCallerRunTable->insert( rsSvCallerName, rsSvCallerDesc, uiJumpRunId );
    }

    inline int64_t insertSvJumpRun( std::string rsSvCallerName, std::string rsSvCallerDesc )
    {
        return pSvJumpRunTable->insert( rsSvCallerName, rsSvCallerDesc );
    }

    inline std::shared_ptr<Pack> reconstructSequencedGenome( std::shared_ptr<Pack> pRef, int64_t iCallerRun )
    {
        return pSvCallTable->reconstructSequencedGenome( pRef, iCallerRun );
    } // method

    inline void updateCoverage( SvCall& rCall )
    {
        pSvCallTable->updateCoverage( rCall );
    } // method

    inline std::shared_ptr<NucSeq> getRead( int64_t iId )
    {
        return pReadTable->getRead( iId );
    } // method

    inline void deleteRun( int64_t uiCallerId )
    {
        CppSQLiteExtImmediateTransactionContext xTransactionContext( *pDatabase );
        // delete rows that can be found easily
        CppSQLiteExtStatement( *pDatabase, "DELETE FROM sv_caller_run_table WHERE id == ?" )
            .bindAndExecute( uiCallerId );
        CppSQLiteExtStatement( *pDatabase, "DELETE FROM sv_call_table WHERE sv_caller_run_id == ?" )
            .bindAndExecute( uiCallerId );
        CppSQLiteExtStatement( *pDatabase, "DELETE FROM sv_call_r_tree WHERE run_id_a == ? " )
            .bindAndExecute( uiCallerId );

        // vacuum up dangeling objects -> compile sqlite with outher switches so this is not needed anymore (@todo)
        CppSQLiteExtStatement( *pDatabase,
                               "DELETE FROM sv_call_support_table WHERE call_id NOT IN (SELECT call_id FROM "
                               "sv_call_table)" )
            .bindAndExecute( );
        CppSQLiteExtStatement(
            *pDatabase,
            "DELETE FROM sv_jump_run_table WHERE id NOT IN (SELECT sv_jump_run_id FROM sv_caller_run_table)" )
            .bindAndExecute( );
        CppSQLiteExtStatement(
            *pDatabase, "DELETE FROM sv_jump_table WHERE sv_jump_run_id NOT IN (SELECT id FROM sv_jump_run_table)" )
            .bindAndExecute( );
    } // method

    class ReadInserter
    {
      private:
        // this is here so that it gets destructed after the transaction context
        std::shared_ptr<SV_DB> pDB;
        // must be after the DB so that it is deconstructed first
        CppSQLiteExtImmediateTransactionContext xTransactionContext;

      public:
        int64_t uiSequencerId;

        ReadInserter( std::shared_ptr<SV_DB> pDB, std::string sSequencerName, std::shared_ptr<Pack> pPack )
            : pDB( pDB ),
              xTransactionContext( *pDB->pDatabase ),
              uiSequencerId( pDB->pSequencerTable->insertSequencer( sSequencerName ) )
        {
            pDB->pContigCovTable->insert( uiSequencerId, pPack );
        } // constructor

        ReadInserter( const ReadInserter& rOther ) = delete; // delete copy constructor

        inline void insertRead( std::shared_ptr<NucSeq> pRead )
        {
            pDB->pReadTable->insertRead( uiSequencerId, pRead );
        } // method

        inline void insertPairedRead( std::shared_ptr<NucSeq> pReadA, std::shared_ptr<NucSeq> pReadB )
        {
            pDB->pPairedReadTable->insertRead( uiSequencerId, pReadA, pReadB );
        } // method

        inline void insertFastaFiles( const ParameterSetManager& rParameters, const std::vector<fs::path>& vsFileNames )
        {
            FileListReader xReader( rParameters, vsFileNames );
            {
                ThreadPool xPool( 4 );
                std::mutex xReadLock, xWriteLock;

                xPool.enqueue(
                    []( size_t uiTid, FileListReader* pReader, ReadInserter* pInserter, std::mutex* pReadLock,
                        std::mutex* pWriteLock ) {
                        while( !pReader->isFinished( ) )
                        {
                            std::shared_ptr<NucSeq> pRead;
                            {
                                std::lock_guard<std::mutex> xGuard( *pReadLock );
                                pRead = pReader->execute( );
                            } // scope for xGuard
                            {
                                std::lock_guard<std::mutex> xGuard( *pWriteLock );
                                pInserter->insertRead( pRead );
                            } // scope for xGuard
                        } // while
                        return 0;
                    },
                    &xReader, this, &xReadLock, &xWriteLock );
            } // scope for thread pool
        } // method

        inline void insertPairedFastaFiles( const ParameterSetManager& rParameters,
                                            const std::vector<fs::path>& vsFileNames1,
                                            const std::vector<fs::path>& vsFileNames2 )
        {
            PairedFileReader xReader( rParameters, vsFileNames1, vsFileNames2 );
            {
                ThreadPool xPool( 4 );
                std::mutex xReadLock, xWriteLock;

                xPool.enqueue(
                    []( size_t uiTid, PairedFileReader* pReader, ReadInserter* pInserter, std::mutex* pReadLock,
                        std::mutex* pWriteLock ) {
                        while( !pReader->isFinished( ) )
                        {
                            std::shared_ptr<TP_PAIRED_READS> pvReads;
                            {
                                std::lock_guard<std::mutex> xGuard( *pReadLock );
                                pvReads = pReader->execute( );
                            } // scope for xGuard
                            {
                                std::lock_guard<std::mutex> xGuard( *pWriteLock );
                                pInserter->insertPairedRead( ( *pvReads )[ 0 ], ( *pvReads )[ 1 ] );
                            } // scope for xGuard
                        } // while
                        return 0;
                    },
                    &xReader, this, &xReadLock, &xWriteLock );
            } // scope for thread pool
        } // method
    }; // class

    class SvJumpInserter
    {
        // this is here so that it gets destructed after the transaction context
        std::shared_ptr<SV_DB> pDB;
        // must be after the DB so that it is deconstructed first
        CppSQLiteExtImmediateTransactionContext xTransactionContext;

      public:
        const int64_t iSvJumpRunId;

        class ReadContex
        {
          private:
            std::shared_ptr<SvJumpTable> pSvJumpTable;
            const int64_t iSvJumpRunId;
            const int64_t iReadId;

          public:
            ReadContex( std::shared_ptr<SvJumpTable> pSvJumpTable, const int64_t iSvJumpRunId, const int64_t iReadId )
                : pSvJumpTable( pSvJumpTable ), iSvJumpRunId( iSvJumpRunId ), iReadId( iReadId )
            {} // constructor

            inline void insertJump( SvJump& rJump )
            {
                // make sure the read id matches the read context
                if( rJump.iReadId == -1 ) // if there is no read id given yet add it
                    rJump.iReadId = iReadId;
                else // otherwise assert it matches
                    assert( rJump.iReadId == iReadId );

                if( rJump.does_switch_strand( ) )
                    assert( rJump.from_start( ) >= std::numeric_limits<int64_t>::max( ) / 2 );
                rJump.iId = pSvJumpTable->xInsertRow(
                    iSvJumpRunId, rJump.iReadId, rJump.from_start( ), rJump.from_end( ), (uint32_t)rJump.uiFrom,
                    (uint32_t)rJump.uiTo, (uint32_t)rJump.uiQueryFrom, (uint32_t)rJump.uiQueryTo,
                    (uint32_t)rJump.uiNumSupportingNt, rJump.bFromForward, rJump.bToForward, rJump.bFromSeedStart );
            } // method
        }; // class

        SvJumpInserter( std::shared_ptr<SV_DB> pDB, int64_t iSvJumpRunId )
            : pDB( pDB ), xTransactionContext( *pDB->pDatabase ), iSvJumpRunId( iSvJumpRunId )
        {} // constructor

        SvJumpInserter( std::shared_ptr<SV_DB> pDB,
                        const std::string& rsSvCallerName,
                        const std::string& rsSvCallerDesc )
            : pDB( pDB ),
              xTransactionContext( *pDB->pDatabase ),
              iSvJumpRunId( pDB->pSvJumpRunTable->insert( rsSvCallerName, rsSvCallerDesc ) )
        {} // constructor

        inline ReadContex readContext( int64_t iReadId )
        {
            return ReadContex( pDB->pSvJumpTable, iSvJumpRunId, iReadId );
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

            inline void remSupport( )
            {
                pSvCallSupportTable->deleteCall( iCallId );
            } // method
        }; // class

        SvCallInserter( std::shared_ptr<SV_DB> pDB, const int64_t iSvCallerRunId )
            : pDB( pDB ), xTransactionContext( *pDB->pDatabase ), iSvCallerRunId( iSvCallerRunId )
        {} // constructor

        SvCallInserter( std::shared_ptr<SV_DB> pDB,
                        const std::string& rsSvCallerName,
                        const std::string& rsSvCallerDesc,
                        const int64_t uiJumpRunId )
            : SvCallInserter( pDB, pDB->pSvCallerRunTable->insert( rsSvCallerName, rsSvCallerDesc, uiJumpRunId ) )
        {} // constructor

        SvCallInserter( const SvCallInserter& ) = delete; // delete copy constructor

        inline void insertCall( SvCall& rCall )
        {
            CallContex xContext( pDB->pSvCallSupportTable, pDB->pSvCallTable->insertCall( iSvCallerRunId, rCall ) );
            for( int64_t iId : rCall.vSupportingJumpIds )
                xContext.addSupport( iId );
        } // method

        inline void updateCall( SvCall& rCall )
        {
            CallContex xContext( pDB->pSvCallSupportTable, pDB->pSvCallTable->updateCall( iSvCallerRunId, rCall ) );
            // remove the link between jumps and this call
            xContext.remSupport( );
            // reinsert the link (no need to compare old and new set this way)
            for( int64_t iId : rCall.vSupportingJumpIds )
                xContext.addSupport( iId );
        } // method

    }; // class

    inline uint32_t numJumps( )
    {
        return pSvJumpTable->numJumps( );
    } // method

    inline uint32_t numCalls( )
    {
        return pSvCallTable->numCalls( );
    } // method

}; // class

// @todo does buffering in vector increase the speed here?
class SortedSvJumpFromSql
{
    const std::shared_ptr<Presetting> pSelectedSetting;
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t,
                               int64_t>
        xQueryStart;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t,
                               int64_t>
        xQueryEnd;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t,
                               int64_t>::Iterator xTableIteratorStart;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t,
                               int64_t>::Iterator xTableIteratorEnd;

  public:
    SortedSvJumpFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerRunId )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQueryStart( *pDb->pDatabase,
                       "SELECT sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                       "       from_seed_start, num_supporting_nt, id, read_id "
                       "FROM sv_jump_table "
                       "WHERE sv_jump_run_id == ? "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd( *pDb->pDatabase,
                     "SELECT sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                     "       from_seed_start, num_supporting_nt, id, read_id "
                     "FROM sv_jump_table "
                     "WHERE sv_jump_run_id == ? "
                     "ORDER BY sort_pos_end" ),
          xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator( iSvCallerRunId ) ),
          xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId ) )
    {} // constructor

    SortedSvJumpFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerRunId,
                         int64_t iX, int64_t iY, uint32_t uiW, uint32_t uiH )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQueryStart( *pDb->pDatabase,
                       "SELECT sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                       "       from_seed_start, num_supporting_nt, id, read_id "
                       "FROM sv_jump_table "
                       "WHERE sv_jump_run_id == ? "
                       "AND ( (from_pos >= ? AND from_pos <= ?) OR from_pos == ? ) "
                       "AND ( (to_pos >= ? AND to_pos <= ?) OR to_pos == ? ) "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd( *pDb->pDatabase,
                     "SELECT sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                     "       from_seed_start, num_supporting_nt, id, read_id "
                     "FROM sv_jump_table "
                     "WHERE sv_jump_run_id == ? "
                     "AND ( (from_pos >= ? AND from_pos <= ?) OR from_pos == ? ) "
                     "AND ( (to_pos >= ? AND to_pos <= ?) OR to_pos == ? ) "
                     "ORDER BY sort_pos_end" ),
          xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator(
              iSvCallerRunId, iX, iX + uiW, std::numeric_limits<uint32_t>::max( ), iY, iY + uiH,
              std::numeric_limits<uint32_t>::max( ) ) ),
          xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId, iX, iX + uiW,
                                                                  std::numeric_limits<uint32_t>::max( ), iY, iY + uiH,
                                                                  std::numeric_limits<uint32_t>::max( ) ) )
    {} // constructor

    SortedSvJumpFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerRunId,
                         int64_t iS, int64_t iE )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQueryStart( *pDb->pDatabase,
                       "SELECT sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                       "       from_seed_start, num_supporting_nt, id, read_id "
                       "FROM sv_jump_table "
                       "WHERE sv_jump_run_id == ? "
                       "AND sort_pos_start >= ? "
                       "AND sort_pos_start <= ? "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd( *pDb->pDatabase,
                     "SELECT sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward, to_forward, "
                     "      from_seed_start, num_supporting_nt, id, read_id "
                     "FROM sv_jump_table "
                     "WHERE sv_jump_run_id == ? "
                     "AND sort_pos_end >= ? "
                     "AND sort_pos_end <= ? "
                     "ORDER BY sort_pos_end" ),
          xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator( iSvCallerRunId, iS, iE ) ),
          xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId, iS, iE ) )
    {
        assert( iE >= iS );
#if DEBUG_LEVEL > 0
#if 0
        std::cout << "SortedSvJumpFromSql::xQueryStart" << std::endl;
        xQueryStart.bindAndExplain( iSvCallerRunId, iX, iX + iW );
        std::cout << "SortedSvJumpFromSql::xQueryEnd" << std::endl;
        xQueryEnd.bindAndExplain( iSvCallerRunId, iX, iX + iW );
#endif
#endif
    } // constructor

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
        if( !hasNextStart( ) )
            return false;
        if( !hasNextEnd( ) )
            return true;
        auto xStartTup = xTableIteratorStart.get( );
        auto xEndTup = xTableIteratorEnd.get( );
        return std::get<0>( xStartTup ) <= std::get<0>( xEndTup );
    } // method

    std::shared_ptr<SvJump> getNextStart( )
    {
        assert( hasNextStart( ) );

        auto xTup = xTableIteratorStart.get( );
        xTableIteratorStart.next( );
        return std::make_shared<SvJump>( pSelectedSetting, std::get<1>( xTup ), std::get<2>( xTup ),
                                         std::get<3>( xTup ), std::get<4>( xTup ), std::get<5>( xTup ),
                                         std::get<6>( xTup ), std::get<7>( xTup ), std::get<8>( xTup ),
                                         std::get<9>( xTup ), std::get<10>( xTup ) );
    } // method

    std::shared_ptr<SvJump> getNextEnd( )
    {
        assert( hasNextEnd( ) );

        auto xTup = xTableIteratorEnd.get( );
        xTableIteratorEnd.next( );
        return std::make_shared<SvJump>( pSelectedSetting, std::get<1>( xTup ), std::get<2>( xTup ),
                                         std::get<3>( xTup ), std::get<4>( xTup ), std::get<5>( xTup ),
                                         std::get<6>( xTup ), std::get<7>( xTup ), std::get<8>( xTup ),
                                         std::get<9>( xTup ), std::get<10>( xTup ) );
    } // method

}; // class

// @todo this does not make the process sufficiently parallel
class AllNucSeqFromSql : public Module<NucSeq, true>
{
    /// this wrapper is required so that the iterator is never copied
    class InteratorHolder
    {
      public:
        CppSQLiteExtQueryStatement<NucSeqSql, int64_t>::Iterator xIterator;

        InteratorHolder( CppSQLiteExtQueryStatement<NucSeqSql, int64_t>& xQuery, int64_t iSequencerId, uint32_t uiRes,
                         uint32_t uiModulo )
            : xIterator( xQuery.vExecuteAndReturnIterator( iSequencerId, uiModulo, uiRes ) )
        {} // constructor

        InteratorHolder( CppSQLiteExtQueryStatement<NucSeqSql, int64_t>& xQuery, int64_t iSequencerId )
            : xIterator( xQuery.vExecuteAndReturnIterator( iSequencerId ) )
        {} // constructor

        InteratorHolder( CppSQLiteExtQueryStatement<NucSeqSql, int64_t>& xQuery )
            : xIterator( xQuery.vExecuteAndReturnIterator( ) )
        {} // constructor
    }; // class
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, int64_t> xQuery;
    std::shared_ptr<InteratorHolder> pTableIterator;
    int64_t iSequencerId;
    uint32_t uiRes;
    uint32_t uiModulo;

  public:
    AllNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
                  "SELECT read_table.sequence, read_table.id "
                  "FROM read_table " ),
          iSequencerId( -1 ),
          uiRes( 0 ),
          uiModulo( 0 )
    {} // constructor

    AllNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSequencerId,
                      size_t uiRes, size_t uiModulo )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
                  ( uiModulo != 1 ? "SELECT read_table.sequence, read_table.id "
                                    "FROM read_table "
                                    "WHERE sequencer_id == ? "
                                    "AND read_table.id % ? == ? "
                                  : "SELECT read_table.sequence, read_table.id "
                                    "FROM read_table "
                                    "WHERE sequencer_id == ? " ) ),
          iSequencerId( iSequencerId ),
          uiRes( (uint32_t)uiRes ),
          uiModulo( (uint32_t)uiModulo )
    {
#if DEBUG_LEVEL > 0
#if 0
        std::cout << "AllNucSeqFromSql::xQuery" << std::endl;
        xQuery.bindAndExplain( iSequencerId, uiModulo, uiRes );
#endif
#endif
    } // constructor

    std::shared_ptr<NucSeq> execute( )
    {
        if( pTableIterator == nullptr && iSequencerId != -1 && uiModulo != 1 )
            pTableIterator = std::make_unique<InteratorHolder>( xQuery, iSequencerId, uiRes, uiModulo );
        else if( pTableIterator == nullptr && iSequencerId != -1 )
            pTableIterator = std::make_unique<InteratorHolder>( xQuery, iSequencerId );
        else if( pTableIterator == nullptr )
            pTableIterator = std::make_unique<InteratorHolder>( xQuery );

        if( pTableIterator->xIterator.eof( ) )
            throw AnnotatedException( "No more NucSeq in NucSeqFromSql module" );

        auto xTup = pTableIterator->xIterator.get( );
        // std::get<0>( xTup ).pNucSeq->sName = std::to_string( std::get<1>( xTup ) );
        std::get<0>( xTup ).pNucSeq->iId = std::get<1>( xTup );
        pTableIterator->xIterator.next( );

        if( pTableIterator->xIterator.eof( ) )
            setFinished( );
        return std::get<0>( xTup ).pNucSeq;
    } // method
}; // class

class NucSeqFromSql : public Module<NucSeq, true>
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t> xQuery;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t>::Iterator xTableIterator;

  public:
    NucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSequencerId )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
                  "SELECT read_table.sequence, read_table.id "
                  "FROM read_table "
                  "WHERE read_table.id NOT IN ( "
                  "   SELECT paired_read_table.first_read FROM paired_read_table "
                  "   UNION "
                  "   SELECT paired_read_table.second_read FROM paired_read_table "
                  ") "
                  "AND sequencer_id = ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSequencerId ) )
    {
        if( xTableIterator.eof( ) )
            setFinished( );
    } // constructor

    NucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
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
    PairedNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSequencerId )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
                  "SELECT A.sequence, B.sequence, A.id, B.id "
                  "FROM read_table A, read_table B "
                  "INNER JOIN paired_read_table "
                  "ON paired_read_table.first_read == A.id "
                  "AND paired_read_table.second_read == B.id "
                  "AND A.sequencer_id = ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSequencerId ) ),
          bRevCompMate( rParameters.getSelected( )->xRevCompPairedReadMates->get( ) )
    {
        if( xTableIterator.eof( ) )
            setFinished( );
    } // constructor

    PairedNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
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

  public:
    // this creates a transaction
    SV_DB::SvJumpInserter xInserter;

    SvDbInserter( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, std::string sRunDesc )
        : pDb( pDb ), xInserter( this->pDb, "MA-SV", sRunDesc )
    {} // constructor

    std::shared_ptr<Container> execute( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        std::lock_guard<std::mutex> xGuard( *pDb->pWriteLock );

        SV_DB::SvJumpInserter::ReadContex xReadContext = xInserter.readContext( pRead->iId );
        for( SvJump& rJump : *pJumps )
            xReadContext.insertJump( rJump ); // also updates the jump ids;

        return std::make_shared<Container>( );
        // end of score for xGuard
    } // method
}; // class

class BufferedSvDbInserter : public Module<Container, false, ContainerVector<SvJump>, NucSeq>
{
    std::shared_ptr<SV_DB> pDb;
    int64_t iSvJumpRunId;

  public:
    std::vector<std::pair<std::shared_ptr<ContainerVector<SvJump>>, int64_t>> vBuffer;
    // this creates a transaction

    BufferedSvDbInserter( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvJumpRunId )
        : pDb( pDb ), iSvJumpRunId( iSvJumpRunId )
    {} // constructor

    inline void commit( )
    {
        if( vBuffer.size( ) == 0 )
            return;
        SV_DB::SvJumpInserter xInserter( pDb, iSvJumpRunId );
        std::lock_guard<std::mutex> xGuard( *pDb->pWriteLock );
        for( auto xPair : vBuffer )
        {
            SV_DB::SvJumpInserter::ReadContex xReadContext = xInserter.readContext( xPair.second );
            for( SvJump& rJump : *xPair.first )
                xReadContext.insertJump( rJump ); // also updates the jump ids;
        } // for
        vBuffer.clear( );
        // end of scope for lock guard
    } // method

    ~BufferedSvDbInserter( )
    {
        commit( );
    } // destructor

    std::shared_ptr<Container> execute( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        vBuffer.emplace_back( pJumps, pRead->iId );
        return std::make_shared<Container>( );
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
    const std::shared_ptr<Presetting> pSelectedSetting;
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, uint32_t, uint32_t>
        xQuery;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, uint32_t, int64_t, int64_t>
        xQuerySupport;
    CppSQLiteExtQueryStatement<int64_t, uint32_t, uint32_t, uint32_t, uint32_t, bool, NucSeqSql, uint32_t,
                               uint32_t>::Iterator xTableIterator;

  public:
    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerId )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_nt, "
                  "       coverage "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id == ? " ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerId ) )
    {} // constructor

    // fetch overlapping or non overlapping calls:
    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerIdA,
                   int64_t iSvCallerIdB, bool bOverlapping, int64_t iAllowedDist )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQuery(
              *pDb->pDatabase,
              ( std::string(
                    "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_nt, "
                    "       coverage "
                    "FROM sv_call_table AS inner "
                    "WHERE sv_caller_run_id == ? "
                    "AND " ) +
                ( bOverlapping ? "" : "NOT " ) +
                // make sure that inner overlaps the outer:
                "EXISTS( "
                "     SELECT outer.id "
                "     FROM sv_call_table AS outer, sv_call_r_tree AS idx_outer "
                "     WHERE outer.id == idx_outer.id "
                "     AND idx_outer.run_id_b >= ? " // dim 1
                "     AND idx_outer.run_id_a <= ? " // dim 1
                "     AND idx_outer.maxX >= inner.from_pos - ? " // dim 2
                "     AND idx_outer.minX <= inner.from_pos + inner.from_size + ? " // dim 2
                "     AND idx_outer.maxY >= inner.to_pos - ? " // dim 3
                "     AND idx_outer.minY <= inner.to_pos + inner.to_size + ? " // dim 3
                "     AND outer.switch_strand == inner.switch_strand "
                ") "
                // make sure that inner does not overlap with any other call with higher score
                "AND NOT EXISTS( "
                "     SELECT inner2.id "
                "     FROM sv_call_table AS inner2, sv_call_r_tree AS idx_inner2 "
                "     WHERE inner2.id == idx_inner2.id "
                "     AND idx_inner2.id != inner.id "
                "     AND (inner2.supporting_nt*1.0)/inner2.coverage >= (inner.supporting_nt*1.0)/inner.coverage "
                "     AND idx_inner2.run_id_b >= inner.id " // dim 1
                "     AND idx_inner2.run_id_a <= inner.id " // dim 1
                "     AND idx_inner2.maxX >= inner.from_pos - ? " // dim 2
                "     AND idx_inner2.minX <= inner.from_pos + inner.from_size + ? " // dim 2
                "     AND idx_inner2.maxY >= inner.to_pos - ? " // dim 3
                "     AND idx_inner2.minY <= inner.to_pos + inner.to_size + ? " // dim 3
                "     AND inner2.switch_strand == inner.switch_strand "
                ") " )
                  .c_str( ) ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerIdA, iSvCallerIdB, iSvCallerIdB, iAllowedDist,
                                                            iAllowedDist, iAllowedDist, iAllowedDist, iAllowedDist,
                                                            iAllowedDist, iAllowedDist, iAllowedDist ) )
    {} // constructor

    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerId,
                   double dMinScore )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_nt, "
                  "       coverage "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id == ? "
                  "AND (supporting_nt*1.0)/coverage >= ? " ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerId, dMinScore ) )
    {} // constructor

    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerId,
                   uint32_t uiX, uint32_t uiY, uint32_t uiW, uint32_t uiH )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_nt, "
                  "       coverage "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id == ? "
                  "AND from_pos + from_size >= ? "
                  "AND to_pos + to_size >= ? "
                  "AND from_pos <= ? "
                  "AND to_pos <= ? " ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerId, uiX, uiY, uiX + uiW, uiY + uiH ) )
    {} // constructor

    SvCallsFromDb( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvCallerId, int64_t iX,
                   int64_t iY, int64_t iW, int64_t iH, double dMinScore )
        : pSelectedSetting( rParameters.getSelected( ) ),
          pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT id, from_pos, to_pos, from_size, to_size, switch_strand, inserted_sequence, supporting_nt, "
                  "       coverage "
                  "FROM sv_call_table "
                  "WHERE sv_caller_run_id == ? "
                  "AND from_pos + from_size >= ? "
                  "AND to_pos + to_size >= ? "
                  "AND from_pos <= ? "
                  "AND to_pos <= ? "
                  "AND (supporting_nt*1.0)/coverage >= ? " ),
          xQuerySupport( *pDb->pDatabase,
                         "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                         "num_supporting_nt, sv_jump_table.id, read_id "
                         "FROM sv_call_support_table "
                         "JOIN sv_jump_table ON sv_call_support_table.jump_id == sv_jump_table.id "
                         "WHERE sv_call_support_table.call_id == ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSvCallerId, iX, iY, iX + iW, iY + iH, dMinScore ) )
    {} // constructor

    SvCall next( )
    {
        auto xTup = xTableIterator.get( );
        SvCall xRet( std::get<1>( xTup ), // uiFromStart
                     std::get<2>( xTup ), // uiToStart
                     std::get<3>( xTup ), // uiFromSize
                     std::get<4>( xTup ), // uiToSize
                     std::get<5>( xTup ), // bSwitchStrand
                     std::get<7>( xTup ) // num_supporting_nt
        );
        xRet.uiCoverage = std::get<8>( xTup );
        xRet.pInsertedSequence = std::get<6>( xTup ).pNucSeq;
        xRet.iId = std::get<0>( xTup );
        auto xSupportIterator( xQuerySupport.vExecuteAndReturnIterator( std::get<0>( xTup ) ) );
        while( !xSupportIterator.eof( ) )
        {
            auto xTup = xSupportIterator.get( );
            xRet.vSupportingJumpIds.push_back( std::get<7>( xTup ) );
            xRet.vSupportingJumps.push_back( std::make_shared<SvJump>(
                pSelectedSetting, std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ), std::get<7>( xTup ), std::get<8>( xTup ),
                std::get<9>( xTup ) ) );
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
