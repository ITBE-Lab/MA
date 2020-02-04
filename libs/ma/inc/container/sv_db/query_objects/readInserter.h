/**
 * @file readInserter.h
 * @brief implements libMA::ReadInserter that inserts reads into the DB.
 * @author Markus Schmidt
 */
// The order of the following two includes is significant with MSVC.

#pragma once

#include "container/sv_db/svSchema.h"
#include "module/get_inserter_container_module.h"

namespace libMA
{

/// @brief inserts reads into a DB
template <typename DBCon> class ReadInserterContainer : public InserterContainer<DBCon, ReadTable, NucSeq>
{
  public:
    using InserterContainer<DBCon, ReadTable, NucSeq>::InserterContainer;

    virtual void insert( std::shared_ptr<NucSeq> pRead )
    {
        pInserter->insert( nullptr, iId, pRead->sName, NucSeqSql( pRead ) );
    } // method
}; // class

/// @brief inserts reads into a DB
template <typename DBCon> class PairedReadInserterContainer : public InserterContainer<DBCon, ReadTable, NucSeq, NucSeq>
{
  public:
    virtual void insert( std::shared_ptr<NucSeq> pReadA, std::shared_ptr<NucSeq> pReadB )
    {
        pInserter->insert( nullptr, iId, pReadA->sName, NucSeqSql( pReadA ) );
        pInserter->insert( nullptr, iId, pReadB->sName, NucSeqSql( pReadB ) );
        // @todo paired connection is lost here!
    } // method
}; // class

template <typename DBCon, typename DBConInit>
using GetReadInserterContainerModule =
    GetInserterContainerModule<ReadInserterContainer, DBCon, DBConInit, SequencerTable, std::string>;
template <typename DBCon> using ReadInserterModule = InserterModule<ReadInserterContainer<DBCon>, NucSeq>;

template <typename DBCon, typename DBConInit>
using GetPairedReadInserterContainerModule =
    GetInserterContainerModule<PairedReadInserterContainer, DBCon, DBConInit, SequencerTable, std::string>;
template <typename DBCon>
using PairedReadInserterModule = InserterModule<PairedReadInserterContainer<DBCon>, NucSeq, NucSeq>;

} // namespace libMA


#ifdef WITH_PYTHON
/// @brief expose libMA::ReadInsert to python
void exportReadInserter( py::module& rxPyModuleId );
#endif
