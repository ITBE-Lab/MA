/**
 * @file svDb.h
 * @details
 * The database interface for the structural variant caller
 */
#pragma once
// FIXME: These two should appear in later files ...
// The order of the following two includes is significant with MSVC.
// #include "db_sql.h" // SQLite connector
#include <MySQL_con.h> // MySQL connector

#include <chrono>
#include <ctime>
#include <iomanip>
#include <string>

#include "container/container.h"
#include "container/nucSeq.h"
#include "container/pack.h"
#include "container/soc.h"
#include "container/svJump.h"
#include "exception.h"
#include "module/fileReader.h"
#include "module/module.h"
#include "system.h"

// include all table definitions
#include "container/sv_db/tables/nameDesc.h"
#include "container/sv_db/tables/pairedRead.h"
#include "container/sv_db/tables/read.h"
#include "container/sv_db/tables/sequencer.h"
#include "container/sv_db/tables/svCall.h"
#include "container/sv_db/tables/svCallSupport.h"
#include "container/sv_db/tables/svCallerRun.h" // cleaned
#include "container/sv_db/tables/svJump.h"

namespace libMA
{

/** @brief: An instance of this object models a connection to a database object.
 */
template <typename DBCon> class SV_Schema : public Container
{
    using DBType = SV_Schema<DBCon>;

  public:
    const std::string sName; // Improvement: DB-Definition via Json
    const bool bInMemory;
    const bool bIsCopy;
    std::shared_ptr<std::mutex> pWriteLock; // For synchronization among different connections
    std::shared_ptr<DBCon> pDatabase; // Pointer to Database itself
    std::shared_ptr<SequencerTable<DBCon>> pSequencerTable;
    std::shared_ptr<ReadTable<DBCon>> pReadTable;
    std::shared_ptr<PairedReadTable<DBCon>> pPairedReadTable;
    std::shared_ptr<NameDescTable<DBCon>> pSvJumpRunTable;
    std::shared_ptr<SvJumpTable<DBCon>> pSvJumpTable;
    std::shared_ptr<SvCallerRunTable<DBCon>> pSvCallerRunTable;
    //- std::shared_ptr<SvCallRegExTable> pSvCallRegExTable;
    std::shared_ptr<SvCallTable<DBCon>> pSvCallTable;
    std::shared_ptr<SvCallSupportTable<DBCon>> pSvCallSupportTable;

    /// @brief create a new database connection
    SV_Schema( std::shared_ptr<DBCon> pDBPoolCon )
        : sName( "unknown" ),
          bInMemory( false ),
          bIsCopy( false ),
          pWriteLock( std::make_shared<std::mutex>( ) ),
          pDatabase( pDBPoolCon ),
          pSequencerTable( std::make_shared<SequencerTable<DBCon>>( pDatabase ) ),
          pReadTable( std::make_shared<ReadTable<DBCon>>( pDatabase ) ),
          pPairedReadTable( std::make_shared<PairedReadTable<DBCon>>( pDatabase, pReadTable ) ),
          pSvJumpRunTable( std::make_shared<NameDescTable<DBCon>>( pDatabase, "sv_jump_run_table" ) ),
          pSvJumpTable( std::make_shared<SvJumpTable<DBCon>>( pDatabase ) ),
          pSvCallerRunTable( std::make_shared<SvCallerRunTable<DBCon>>( pDatabase ) ),
          //- pSvCallRegExTable( std::make_shared<SvCallRegExTable>( pDatabase ) ),
          pSvCallTable( std::make_shared<SvCallTable<DBCon>>( pDatabase, pWriteLock, sName ) ),
          pSvCallSupportTable( std::make_shared<SvCallSupportTable<DBCon>>( pDatabase ) )
    {
        // DEBUG( std::cout << "Opened DB connection" << std::endl; )
        //--> Json SQLite : this->setNumThreads( 32 ); // @todo do this via a parameter
        //--> Json SQLite : pDatabase->execDML( "PRAGMA busy_timeout=0;" ); // do not throw sqlite busy errors
        //--> Json SQLite : //
        // https://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
        //--> Json SQLite : pDatabase->execDML( "PRAGMA synchronous = OFF;" ); // insert performance
        //--> Json SQLite : pDatabase->execDML( "PRAGMA journal_mode = WAL;" ); // insert performance -> read while
        // write
    } // constructor

    // Delete: SV_DB( std::string sName ) : SV_DB( sName, eCREATE_DB, false )
    // Delete: {} // constructor

    SV_Schema( std::string sName, std::string sMode, bool bInMemory ) : SV_Schema( std::make_shared<DBCon>( sName ) )
    {} // constructor

    SV_Schema( std::string sName, std::string sMode ) : SV_Schema( sName, sMode, false )
    {} // constructor

    SV_Schema( std::string sName ) : SV_Schema( sName, "open" )
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
        // FIXME: Unimplemented yet.
        // pDatabase->set_num_threads( (int)uiN );
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
        return pSvCallerRunTable->insert_( rsSvCallerName, rsSvCallerDesc, uiJumpRunId );
    }

    inline int64_t insertSvJumpRun( std::string rsSvCallerName, std::string rsSvCallerDesc )
    {
        // FIXME: CppSQLiteExtImmediateTransactionContext xTransactionContext( *pDatabase );
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
        //- TO DO: CppSQLiteExtImmediateTransactionContext xTransactionContext( *pDatabase );
        // delete rows that can be found easily
        SQLStatement<DBCon>( pDatabase, "DELETE FROM sv_caller_run_table WHERE id = ?" ).exec( uiCallerId );
        SQLStatement<DBCon>( pDatabase, "DELETE FROM sv_call_table WHERE sv_caller_run_id = ?" ).exec( uiCallerId );

        // vacuum up dangeling objects -> compile sqlite with other switches so this is not needed anymore (@todo)
        SQLStatement<DBCon>( pDatabase,
                             "DELETE FROM sv_call_support_table WHERE call_id NOT IN (SELECT call_id FROM "
                             "sv_call_table)" )
            .exec( );
        SQLStatement<DBCon>(
            pDatabase,
            "DELETE FROM sv_jump_run_table WHERE id NOT IN (SELECT sv_jump_run_id FROM sv_caller_run_table)" )
            .exec( );
        SQLStatement<DBCon>(
            pDatabase, "DELETE FROM sv_jump_table WHERE sv_jump_run_id NOT IN (SELECT id FROM sv_jump_run_table)" )
            .exec( );
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
