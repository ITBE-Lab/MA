/**
 * @file callInserter.h
 * @brief Implements libMA::SvCallInserter; a transaction based structural variant call inserter
 * @author Markus Schmidt
 */
#pragma once

#include "container/sv_db/svSchema.h"
#include "module/get_inserter_container_module.h"

namespace libMA
{
/**
 * @brief A transaction based structural variant call inserter
 * @details
 * Objects of this class can be used to update or insert structural variant calls into a libMA::svDb.
 */
template <typename DBCon> class SvCallInserterContainer : public InserterContainer<DBCon, SvCallTable>
{
    const static size_t BUFFER_SIZE = 500;
    using CALL_INSERTER_TYPE = decltype( SvCallTable<DBCon>::template getBulkInserter<BUFFER_SIZE> );
    CALL_INSERTER_TYPE pCallInserter;
    using SUPPORT_INSERTER_TYPE = decltype( SvCallSupportTable<DBCon>::template getBulkInserter<BUFFER_SIZE> );
    SUPPORT_INSERTER_TYPE pSupportInserter;

  public:
    /// @brief id of the caller run this inserter is bound to.
    const int64_t iSvCallerRunId;

    /**
     * @brief create a SvCallInserter object for the given run id.
     * @details
     * Expects the run to exists in the DB.
     */
    SvCallInserterContainer( std::shared_ptr<DBCon> pConnection, const int64_t iSvCallerRunId )
        : pCallInserter( SvCallTable<DBCon>( pConnection ).template getBulkInserter<BUFFER_SIZE>( ) ),
          pSupportInserter( SvCallSupportTable<DBCon>( pConnection ).template getBulkInserter<BUFFER_SIZE>( ) ),
          iSvCallerRunId( iSvCallerRunId )
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
        int64_t iCallId = pCallInserter->insert( nullptr, //
                                                 iSvCallerRunId, //
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
            pSupportInserter->insert( nullptr, iCallId, iId );
    } // method

    inline void insert( std::shared_ptr<SvCall> pCall )
    {
        insertCall( *pCall );
    } // method


    inline void close()
    {
        pInserter->reset();
        pSupportInserter->reset();
    } // method
}; // class

template <typename DBCon, typename DBConInit>
class GetCallInserterContainerModule
    : public GetInserterContainerModule<SvCallInserterContainer, DBCon, DBConInit, SvCallerRunTable>
{}; // class

template <typename DBCon> using CallInserterModule = InserterModule<SvCallInserterContainer<DBCon>, SvCall>;

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief used to expose libMA::SvCallInserter to python
 */
void exportSvCallInserter( py::module& rxPyModuleId );
#endif
