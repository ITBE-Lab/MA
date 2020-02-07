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

/**
 * @brief A transaction based structural variant call inserter
 * @details
 * Objects of this class can be used to update or insert structural variant calls into a libMA::svDb.
 * @todo ask arne if we can somehow use a bulk inserter here
 */
template <typename CallOrVector, typename DBCon>
class SvCallInserterContainerTmpl : public InserterContainer<DBCon, SvCallTable, CallOrVector>
{
  public:
    using ParentType = InserterContainer<DBCon, SvCallTable, CallOrVector>;

    const static size_t BUFFER_SIZE = 500;
    using SupportInserterType = typename SvCallSupportTable<DBCon>::template SQLBulkInserterType<BUFFER_SIZE>;
    std::shared_ptr<SupportInserterType> pSupportInserter;


    SvCallInserterContainerTmpl( std::shared_ptr<ConnectionContainer<DBCon>> pConnection, int64_t iId )
        : ParentType::InserterContainer( pConnection, iId ),
          pSupportInserter( std::make_shared<SupportInserterType>( SvCallSupportTable<DBCon>( pConnection->pConnection ) ) )
    {} // constructor

    /**
     * @brief insert a new call and link to the supporting jumps.
     * @details
     * Makes use of the libMA::SvCallInserter::CallContex.
     * Expects that rCall does NOT exist in the DB yet.
     */
    inline void insertCall( SvCall& rCall )
    {
        auto xRectangle = WKBUint64Rectangle( rCall );
        int64_t iCallId = ParentType::pTable->insert( ParentType::iId, //
                                                      (uint32_t)rCall.xXAxis.start( ), //
                                                      (uint32_t)rCall.xYAxis.start( ), //
                                                      (uint32_t)rCall.xXAxis.size( ), //
                                                      (uint32_t)rCall.xYAxis.size( ), //
                                                      rCall.bSwitchStrand, //
                                                      // NucSeqSql can deal with nullpointers
                                                      NucSeqSql( rCall.pInsertedSequence ), //
                                                      (uint32_t)rCall.uiNumSuppReads, //
                                                      (uint32_t)rCall.uiReferenceAmbiguity, //
                                                      -1, //
                                                      xRectangle );
        rCall.iId = iCallId;
        for( int64_t iId : rCall.vSupportingJumpIds )
            pSupportInserter->insert( iCallId, iId );
    } // method

    virtual void EXPORTED insert( std::shared_ptr<SvCall> pCall )
    {
        insertCall( *pCall );
    } // method

    virtual void EXPORTED insert( std::shared_ptr<CompleteBipartiteSubgraphClusterVector> pCalls )
    {
        for( auto pCall : pCalls->vContent )
            insertCall( *pCall );
    } // method

    virtual void close( )
    {
        pSupportInserter.reset( );
        ParentType::close( );
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

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief used to expose libMA::SvCallInserter to python
 */
void exportSvCallInserter( py::module& rxPyModuleId );
#endif
