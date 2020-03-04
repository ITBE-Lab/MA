/**
 * @file callInserter.h
 * @brief Implements libMA::SvCallInserter; a transaction based structural variant call inserter
 * @author Markus Schmidt
 */
#pragma once

#include "container/svJump.h"
#include "container/sv_db/tables/svCall.h"
#include "container/sv_db/tables/svCallSupport.h"
#include "container/sv_db/tables/svCallerRun.h"
#include "module/get_inserter_container_module.h"

namespace libMA
{


#if DEBUG_LEVEL == 0
// in release mode, we can use a BulkInserter since foreign_key_checks are
// disabled
#define BulkOrNot BulkInserterContainer
#else
// if we are in debug mode, calls need to be inserted immediately since
// the sv_call_support_table references the primary key of the call_table
#define BulkOrNot InserterContainer
#endif

/**
 * @brief A transaction based structural variant call inserter
 * @details
 * Objects of this class can be used to update or insert structural variant calls into a libMA::svDb.
 */
template <typename CallOrVector, typename DBCon>
class SvCallInserterContainerTmpl : public BulkOrNot<DBCon, AbstractInserterContainer, SvCallTable, CallOrVector>
{
  public:
    using ParentType = BulkOrNot<DBCon, AbstractInserterContainer, SvCallTable, CallOrVector>;

    std::shared_ptr<BulkInserterType<SvCallSupportTable<DBCon>>> pSupportInserter;


    SvCallInserterContainerTmpl( std::shared_ptr<PoolContainer<DBCon>> pPool,
                                 int64_t iId,
                                 std::shared_ptr<SharedInserterProfiler>
                                     pSharedProfiler )
        : ParentType::BulkOrNot( pPool, iId, pSharedProfiler ),
          pSupportInserter( pPool->xPool.run( ParentType::iConnectionId, []( auto pConnection ) {
              return SvCallSupportTable<DBCon>( pConnection )
                  .template getBulkInserter<SvCallSupportTable<DBCon>::uiBulkInsertSize::value>( );
          } ) )
    {} // constructor

    /**
     * @brief insert a new call and link to the supporting jumps.
     * @details
     * Makes use of the libMA::SvCallInserter::CallContex.
     * Expects that rCall does NOT exist in the DB yet.
     */
    inline size_t insertCall( SvCall& rCall )
    {
        auto xRectangle = WKBUint64Rectangle( rCall );
        int64_t iCallId = ParentType::pInserter->insert( ParentType::iId, //
                                                         (uint32_t)rCall.xXAxis.start( ), //
                                                         (uint32_t)rCall.xYAxis.start( ), //
                                                         (uint32_t)rCall.xXAxis.size( ), //
                                                         (uint32_t)rCall.xYAxis.size( ), //
                                                         rCall.bSwitchStrand, //
                                                         // can deal with nullpointers
                                                         makeSharedCompNucSeq( rCall.pInsertedSequence ), //
                                                         (uint32_t)rCall.uiNumSuppReads, //
                                                         (uint32_t)rCall.uiReferenceAmbiguity, //
                                                         -1, // regex_id
                                                         -1, // filter_id
                                                         xRectangle );
        rCall.iId = iCallId;
        for( int64_t iId : rCall.vSupportingJumpIds )
            pSupportInserter->insert( iCallId, iId );
        return 1 + rCall.vSupportingJumpIds.size( );
    } // method
  protected:
    virtual size_t EXPORTED insert_override( std::shared_ptr<SvCall> pCall )
    {
        return insertCall( *pCall );
    } // method

    virtual size_t EXPORTED insert_override( std::shared_ptr<CompleteBipartiteSubgraphClusterVector> pCalls )
    {
        size_t uiTotal = 0;
        for( auto pCall : pCalls->vContent )
            uiTotal += insertCall( *pCall );
        return uiTotal;
    } // method
  public:
    virtual void close( std::shared_ptr<PoolContainer<DBCon>> pPool )
    {
        pSupportInserter.reset( );
        ParentType::close( pPool );
    } // method
}; // class

template <typename DBCon> using SvCallInserterContainer = SvCallInserterContainerTmpl<SvCall, DBCon>;
template <typename DBCon>
using SvCallVectorInserterContainer = SvCallInserterContainerTmpl<CompleteBipartiteSubgraphClusterVector, DBCon>;

template <typename DBCon, typename DBConInit>
using GetCallInserterContainerModule =
    GetInserterContainerModule<SvCallInserterContainer, DBCon, DBConInit, SvCallerRunTable>;
template <typename DBCon, typename DBConInit>
using GetCallVectorInserterContainerModule =
    GetInserterContainerModule<SvCallVectorInserterContainer, DBCon, DBConInit, SvCallerRunTable>;

template <typename DBCon> using SvCallInserterModule = InserterModule<SvCallInserterContainer<DBCon>>;
template <typename DBCon> using SvCallVectorInserterModule = InserterModule<SvCallVectorInserterContainer<DBCon>>;

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief used to expose libMA::SvCallInserter to python
 */
void exportSvCallInserter( py::module& rxPyModuleId );
#endif
