/**
 * @file callInserter.h
 * @brief Implements libMA::SvCallInserter; a transaction based structural variant call inserter
 * @author Markus Schmidt
 */
#include "container/sv_db/svSchema.h"

#pragma once

namespace libMA
{
/**
 * @brief A transaction based structural variant call inserter
 * @details
 * Objects of this class can be used to update or insert structural variant calls into a libMA::svDb.
 */
template <typename DBCon> class SvCallInserter
{
  public:
    // this is here so that it gets destructed after the transaction context
    std::shared_ptr<DBCon> pConnection;

    std::shared_ptr<SvCallTable<DBCon>> pSvCallTable;
    std::shared_ptr<SvCallSupportTable<DBCon>> pSvCallSupportTable;

  private:
    // must be after the DB so that it is deconstructed first

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
        std::shared_ptr<SvCallSupportTable<DBCon>> pSvCallSupportTable;
        const int64_t iCallId;

      public:
        /**
         * @brief Create a context for the call with id iCallId.
         */
        CallContex( std::shared_ptr<SvCallSupportTable<DBCon>> pSvCallSupportTable, const int64_t iCallId )
            : pSvCallSupportTable( pSvCallSupportTable ), iCallId( iCallId )
        {} // constructor

        /**
         * @brief links the given rJump with the call in context.
         * @details
         * This uses the iId of the jump and expects the jump to exist in the database.
         */
        inline void addSupport( SvJump& rJump )
        {
            pSvCallSupportTable->insert( iCallId, rJump.iId ); // xInsertRow -> insert
        } // method

        /**
         * @brief links the given rJump (via it's DB ID) with the call in context.
         */
        inline void addSupport( int64_t iId )
        {
            pSvCallSupportTable->insert( iCallId, iId ); // xInsertRow -> insert
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
    SvCallInserter( std::shared_ptr<DBCon> pConnection, const int64_t iSvCallerRunId )
        : pConnection( pConnection ),
          pSvCallTable( std::make_shared<SvCallTable<DBCon>>( pConnection ) ),
          pSvCallSupportTable( std::make_shared<SvCallSupportTable<DBCon>>( pConnection ) ),
          // @fixme: pTransactionContext( std::make_shared<CppSQLiteExtImmediateTransactionContext>( *pDB->pDatabase )
          // ),
          iSvCallerRunId( iSvCallerRunId )
    {} // constructor

    /**
     * @brief create a SvCallInserter object for a new run.
     * @details
     * This creates a new caller run with the given name and description.
     */
    SvCallInserter( std::shared_ptr<DBCon> pConnection,
                    const std::string& rsSvCallerName,
                    const std::string& rsSvCallerDesc,
                    const int64_t uiJumpRunId )
        : SvCallInserter(
              pConnection,
              SvCallerRunTable<DBCon>( pConnection ).insert_( rsSvCallerName, rsSvCallerDesc, uiJumpRunId ) )
    {} // constructor

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
        CallContex xContext( pSvCallSupportTable, pSvCallTable->insertCall( iSvCallerRunId, rCall ) );
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
        CallContex xContext( pSvCallSupportTable, pSvCallTable->updateCall( iSvCallerRunId, rCall ) );
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
        // @fixme: NEW DATABASE API
    } // method

    /**
     * @brief terminates the current transaction and starts a new one
     */
    inline void reOpenTransaction( )
    {
        // @fixme: NEW DATABASE API
    } // method

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief used to expose libMA::SvCallInserter to python
 */
void exportSvCallInserter( py::module& rxPyModuleId );
#endif
