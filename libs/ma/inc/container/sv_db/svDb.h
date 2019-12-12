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
#include "exception.h"
#include "db_sql.h"
#include "system.h"
#include <chrono>
#include <ctime>
#include <iomanip>
// include all table definitions
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

        // vacuum up dangeling objects -> compile sqlite with other switches so this is not needed anymore (@todo)
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

    inline uint32_t numJumps( )
    {
        return pSvJumpTable->numJumps( );
    } // method

    inline uint32_t numCalls( )
    {
        return pSvCallTable->numCalls( );
    } // method

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
void exportSoCDbWriter( py::module& rxPyModuleId );
#endif
