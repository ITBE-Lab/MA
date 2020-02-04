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
template <typename DBCon> class SvCallInserterContainer : public InserterContainer<DBCon, ReadTable, SvCall>
{
    const static size_t BUFFER_SIZE = InserterContainer<DBCon, ReadTable, SvCall>::BUFFER_SIZE;
    using SupportInserterType = typename SvCallSupportTable<DBCon>::template SQLBulkInserterType<BUFFER_SIZE>;
    std::shared_ptr<SupportInserterType> pSupportInserter;

    /**
     * @brief insert a new call and link to the supporting jumps.
     * @details
     * Makes use of the libMA::SvCallInserter::CallContex.
     * Expects that rCall does NOT exist in the DB yet.
     */
    inline void insertCall( SvCall& rCall )
    {
        auto xRectangle = WKBUint64Rectangle( rCall );
        int64_t iCallId = InserterContainer<DBCon, ReadTable, SvCall>::pInserter->insert( nullptr, //
                                             InserterContainer<DBCon, ReadTable, SvCall>::iId, //
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

    virtual void insert( std::shared_ptr<SvCall> pCall )
    {
        insertCall( *pCall );
    } // method

    virtual void close( )
    {
        pSupportInserter->reset( );
        InserterContainer<DBCon, ReadTable, SvCall>::close( );
    } // method
}; // class

template <typename DBCon, typename DBConInit>
class GetCallInserterContainerModule : public GetInserterContainerModule<SvCallInserterContainer,
                                                                         DBCon,
                                                                         DBConInit,
                                                                         SvCallerRunTable,
                                                                         std::string, // name
                                                                         std::string, // desc
                                                                         int64_t, // timestamp
                                                                         int64_t> // sv_jump_run_id
{}; // class

template <typename DBCon> using CallInserterModule = InserterModule<SvCallInserterContainer<DBCon>, SvCall>;

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief used to expose libMA::SvCallInserter to python
 */
void exportSvCallInserter( py::module& rxPyModuleId );
#endif
