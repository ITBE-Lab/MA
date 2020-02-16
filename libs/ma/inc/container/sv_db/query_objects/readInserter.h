/**
 * @file readInserter.h
 * @brief implements libMA::ReadInserter that inserts reads into the DB.
 * @author Markus Schmidt
 */
// The order of the following two includes is significant with MSVC.

#pragma once

#include "container/svJump.h"
#include "container/sv_db/tables/pairedRead.h"
#include "container/sv_db/tables/read.h"
#include "container/sv_db/tables/sequencer.h"
#include "module/get_inserter_container_module.h"

namespace libMA
{

/// @brief inserts reads into a DB
template <typename DBCon> class ReadInserterContainer : public BulkInserterContainer<DBCon, ReadTable, NucSeq>
{
  public:
    using ParentType = BulkInserterContainer<DBCon, ReadTable, NucSeq>;
    using ParentType::BulkInserterContainer;

    virtual void EXPORTED insert( std::shared_ptr<NucSeq> pRead )
    {
        ParentType::pInserter->insert( nullptr, ParentType::iId, pRead->sName, NucSeqSql( pRead ) );
    } // method
}; // class

template <typename DBCon, typename DBConInit>
using GetReadInserterContainerModule =
    GetInserterContainerModule<ReadInserterContainer, DBCon, DBConInit, SequencerTable>;
template <typename DBCon> using ReadInserterModule = InserterModule<ReadInserterContainer<DBCon>>;

/**
 * @brief inserts reads into a DB
 * @details
 * @todo ask arne if we can somehow use a bulk inserter here
 */
template <typename DBCon> class PairedReadInserterContainer : public InserterContainer<DBCon, ReadTable, NucSeq, NucSeq>
{
  public:
    using ParentType = InserterContainer<DBCon, ReadTable, NucSeq, NucSeq>;

    const static size_t BUFFER_SIZE = 500;
    std::shared_ptr<BulkInserterType<PairedReadTable<DBCon>>> pPairedReadInserter;

    PairedReadInserterContainer( std::shared_ptr<PoolContainer<DBCon>> pPool, int64_t iId )
        : ParentType::InserterContainer( pPool, iId ),
          pPairedReadInserter( pPool->xPool.run(
              ParentType::iConnectionId,
              []( auto pConnection ) //
              {
                  return std::make_shared<BulkInserterType<PairedReadTable<DBCon>>>( PairedReadTable<DBCon>(
                      pConnection, nullptr /* @todo nullptr is bad practice but save in this context */ ) );
              } ) )
    {} // constructor

    virtual void EXPORTED insert( std::shared_ptr<NucSeq> pReadA, std::shared_ptr<NucSeq> pReadB )
    {
        int64_t iReadIdA = ParentType::pInserter->insert( ParentType::iId, pReadA->sName, NucSeqSql( pReadA ) );
        int64_t iReadIdB = ParentType::pInserter->insert( ParentType::iId, pReadB->sName, NucSeqSql( pReadB ) );
        pPairedReadInserter->insert( iReadIdA, iReadIdB );
    } // method
}; // class

template <typename DBCon, typename DBConInit>
using GetPairedReadInserterContainerModule =
    GetInserterContainerModule<PairedReadInserterContainer, DBCon, DBConInit, SequencerTable>;
template <typename DBCon> using PairedReadInserterModule = InserterModule<PairedReadInserterContainer<DBCon>>;

} // namespace libMA


#ifdef WITH_PYTHON
/// @brief expose libMA::ReadInsert to python
void exportReadInserter( py::module& rxPyModuleId );
#endif
