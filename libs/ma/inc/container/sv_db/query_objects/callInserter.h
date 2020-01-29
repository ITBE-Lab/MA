/**
 * @file callInserter.h
 * @brief Implements libMA::SvCallInserter; a transaction based structural variant call inserter
 * @author Markus Schmidt
 */
#include "db_config.h"
#ifndef USE_NEW_DB_API
#include "container/sv_db/svDb.h"
#else
#include "container/sv_db/_svDb.h" // NEW DATABASE INTERFACE
#endif

#pragma once

namespace libMA
{
/**
 * @brief A transaction based structural variant call inserter
 * @details
 * Objects of this class can be used to update or insert structural variant calls into a libMA::svDb.
 */
#ifdef USE_NEW_DB_API
template <typename DBCon>
#endif
class SvCallInserter
{
  public:
    // this is here so that it gets destructed after the transaction context
#ifndef USE_NEW_DB_API
    std::shared_ptr<SV_DB> pDB;
#else
    std::shared_ptr<_SV_DB<DBCon>> pDB;
#endif
  private:
    // must be after the DB so that it is deconstructed first
#ifndef USE_NEW_DB_API
    std::shared_ptr<CppSQLiteExtImmediateTransactionContext> pTransactionContext;
#else
// FIXME: Insert transaction-contexts in the MySQL interface
#endif

  public:
    /// @brief id of the caller run this inserter is bound to.
    const int64_t iSvCallerRunId;

    /**
     * @brief Context for one call
     * @details
     * Each call is supported by multiple libMA::SvJump.
     * These can be inserted by creating a cotext for the libMA::SvCall.
     */
    class CallContex
    {
      private:
#ifndef USE_NEW_DB_API
        std::shared_ptr<SvCallSupportTable> pSvCallSupportTable;
#else
        std::shared_ptr<SvCallSupportTable<DBCon>> pSvCallSupportTable;
#endif
        const int64_t iCallId;

      public:
        /**
         * @brief Create a context for the call with id iCallId.
         */
#ifndef USE_NEW_DB_API
        CallContex( std::shared_ptr<SvCallSupportTable> pSvCallSupportTable, const int64_t iCallId )
            : pSvCallSupportTable( pSvCallSupportTable ), iCallId( iCallId )
#else
        CallContex( std::shared_ptr<SvCallSupportTable<DBCon>> pSvCallSupportTable, const int64_t iCallId )
            : pSvCallSupportTable( pSvCallSupportTable ), iCallId( iCallId )
#endif
        {} // constructor

        /**
         * @brief links the given rJump with the call in context.
         * @details
         * This uses the iId of the jump and expects the jump to exist in the database.
         */
        inline void addSupport( SvJump& rJump )
        {
#ifndef USE_NEW_DB_API
            pSvCallSupportTable->xInsertRow( iCallId, rJump.iId ); 
#else
			pSvCallSupportTable->insert(iCallId, rJump.iId); // xInsertRow -> insert
#endif
        } // method

        /**
         * @brief links the given rJump (via it's DB ID) with the call in context.
         */
        inline void addSupport( int64_t iId )
        {
#ifndef USE_NEW_DB_API
            pSvCallSupportTable->xInsertRow( iCallId, iId ); 
#else
			pSvCallSupportTable->insert(iCallId, iId); // xInsertRow -> insert
#endif
        } // method

        /**
         * @brief removes ALL links between jumps and the call in context.
         */
        inline void remSupport( )
        {
            pSvCallSupportTable->deleteCall( iCallId );
        } // method
    }; // class

    /**
     * @brief create a SvCallInserter object for the given run id.
     * @details
     * Expects the run to exists in the DB.
     */
#ifndef USE_NEW_DB_API
    SvCallInserter( std::shared_ptr<SV_DB> pDB, const int64_t iSvCallerRunId )
        : pDB( pDB ),
          pTransactionContext( std::make_shared<CppSQLiteExtImmediateTransactionContext>( *pDB->pDatabase ) ),
          iSvCallerRunId( iSvCallerRunId )
    {} // constructor
#else
    SvCallInserter( std::shared_ptr<_SV_DB<DBCon>> pDB, const int64_t iSvCallerRunId )
        : pDB( pDB ),
          // FIXME: pTransactionContext( std::make_shared<CppSQLiteExtImmediateTransactionContext>( *pDB->pDatabase ) ),
          iSvCallerRunId( iSvCallerRunId )
    {} // constructor
#endif

    /**
     * @brief create a SvCallInserter object for a new run.
     * @details
     * This creates a new caller run with the given name and description.
     */
#ifndef USE_NEW_DB_API
    SvCallInserter( std::shared_ptr<SV_DB> pDB,
                    const std::string& rsSvCallerName,
                    const std::string& rsSvCallerDesc,
                    const int64_t uiJumpRunId )
        : SvCallInserter( pDB, pDB->pSvCallerRunTable->insert( rsSvCallerName, rsSvCallerDesc, uiJumpRunId ) )
    {} // constructor
#else
    SvCallInserter( std::shared_ptr<_SV_DB<DBCon>> pDB,
                    const std::string& rsSvCallerName,
                    const std::string& rsSvCallerDesc,
                    const int64_t uiJumpRunId )
        : SvCallInserter( pDB, pDB->pSvCallerRunTable->insert_( rsSvCallerName, rsSvCallerDesc, uiJumpRunId ) )
    {} // constructor
#endif
    /// @brief Object cannot be copied.
    SvCallInserter( const SvCallInserter& ) = delete; // delete copy constructor

    /**
     * @brief insert a new call and link to the supporting jumps.
     * @details
     * Makes use of the libMA::SvCallInserter::CallContex.
     * Expects that rCall does NOT exist in the DB yet.
     */
    inline void insertCall( SvCall& rCall )
    {
        CallContex xContext( pDB->pSvCallSupportTable, pDB->pSvCallTable->insertCall( iSvCallerRunId, rCall ) );
        for( int64_t iId : rCall.vSupportingJumpIds )
            xContext.addSupport( iId );
    } // method

    /**
     * @brief updates a call and all it's links to supporting jumps.
     * @details
     * makes use of the libMA::SvCallInserter::CallContex.
     * Expects that rCall does exist in the DB.
     */
    inline void updateCall( SvCall& rCall )
    {
        CallContex xContext( pDB->pSvCallSupportTable, pDB->pSvCallTable->updateCall( iSvCallerRunId, rCall ) );
        // remove the link between jumps and this call
        xContext.remSupport( );
        // reinsert the link (no need to compare old and new set this way)
        for( int64_t iId : rCall.vSupportingJumpIds )
            xContext.addSupport( iId );
    } // method

    /**
     * @brief terminates the current transaction
     */
    inline void endTransaction( )
    {
#ifndef USE_NEW_DB_API // FIXME: NEW DATABASE API
        pTransactionContext.reset( );
#endif
    }; // method

    /**
     * @brief terminates the current transaction and starts a new one
     */
    inline void reOpenTransaction( )
    {
#ifndef USE_NEW_DB_API // FIXME: NEW DATABASE API
        endTransaction( );
        pTransactionContext = std::make_shared<CppSQLiteExtImmediateTransactionContext>( *pDB->pDatabase );
#endif
    }; // method

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief used to expose libMA::SvCallInserter to python
 */
void exportSvCallInserter( py::module& rxPyModuleId );
#endif
