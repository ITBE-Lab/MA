#include "container/sv_db/svDb.h"

#pragma once

namespace libMA
{

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

} // namespace libMA

#ifdef WITH_PYTHON
void exportSvCallInserter( py::module& rxPyModuleId );
#endif