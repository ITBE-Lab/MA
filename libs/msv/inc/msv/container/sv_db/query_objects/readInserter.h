/**
 * @file readInserter.h
 * @brief implements libMSV::ReadInserter that inserts reads into the DB.
 * @author Markus Schmidt
 */
// The order of the following two includes is significant with MSVC.

#pragma once

#include "ms/module/get_inserter_container_module.h"
#include "msv/container/svJump.h"
#include "msv/container/sv_db/tables/pairedRead.h"
#include "msv/container/sv_db/tables/read.h"
#include "msv/container/sv_db/tables/sequencer.h"


using namespace libMA;
using namespace libMS;

namespace libMSV
{

/// @brief inserts reads into a DB
template <typename DBCon>
class ReadInserterContainer : public BulkInserterContainer<DBCon, AbstractInserterContainer, ReadTable, NucSeq>
{
  public:
    using ParentType = BulkInserterContainer<DBCon, libMS::AbstractInserterContainer, ReadTable, NucSeq>;

// Expose constructor of base class
#if defined( __clang__ )
    ReadInserterContainer( std::shared_ptr<PoolContainer<DBCon>> pPool, int64_t iId,
                           std::shared_ptr<SharedInserterProfiler> pSharedProfiler )
        : ParentType::BulkInserterContainer( pPool, iId, pSharedProfiler )
    {} // constructur
#else
    using ParentType::BulkInserterContainer;
#endif

  protected:
    virtual size_t insert_override( std::shared_ptr<NucSeq> pRead )
    {
        ParentType::pInserter->insert( ParentType::iId, pRead->sName, makeSharedCompNucSeq( *pRead ) );
        return 1;
    } // method
}; // class

template <typename DBCon, typename DBConInit>
using GetReadInserterContainerModule =
    GetInserterContainerModule<ReadInserterContainer, DBCon, DBConInit, SequencerTable>;
template <typename DBCon> using ReadInserterModule = InserterModule<ReadInserterContainer<DBCon>>;

/**
 * @brief inserts reads into a DB
 * @details
 */
template <typename DBCon>
class PairedReadInserterContainer
    : public BulkInserterContainer<DBCon, AbstractInserterContainer, ReadTable, NucSeq, NucSeq>
{
  public:
    using ParentType = BulkInserterContainer<DBCon, libMS::AbstractInserterContainer, ReadTable, NucSeq, NucSeq>;

    std::shared_ptr<BulkInserterType<PairedReadTable<DBCon>>> pPairedReadInserter;

    PairedReadInserterContainer( std::shared_ptr<PoolContainer<DBCon>> pPool,
                                 int64_t iId,
                                 std::shared_ptr<SharedInserterProfiler>
                                     pSharedProfiler )
        : ParentType::BulkInserterContainer( pPool, iId, pSharedProfiler ),
          pPairedReadInserter( pPool->xPool.run(
              ParentType::iConnectionId,
              []( auto pConnection ) //
              {
                  return PairedReadTable<DBCon>( pConnection,
                                                 /* nullptr is save in this context since PairedReadTable only uses it's
                                                    second constructor argument if insert is called. */
                                                 nullptr )
                      .template getBulkInserter<PairedReadTable<DBCon>::uiBulkInsertSize::value>( );
              } ) )
    {} // constructor
  protected:
    virtual size_t insert_override( std::shared_ptr<NucSeq> pReadA, std::shared_ptr<NucSeq> pReadB )
    {
        auto iReadIdA =
            ParentType::pInserter->insert( ParentType::iId, pReadA->sName, makeSharedCompNucSeq( *pReadA ) );
        auto iReadIdB =
            ParentType::pInserter->insert( ParentType::iId, pReadB->sName, makeSharedCompNucSeq( *pReadB ) );
        pPairedReadInserter->insert( iReadIdA, iReadIdB );
        return 3;
    } // method
}; // class

template <typename DBCon, typename DBConInit>
using GetPairedReadInserterContainerModule =
    GetInserterContainerModule<PairedReadInserterContainer, DBCon, DBConInit, SequencerTable>;
template <typename DBCon> using PairedReadInserterModule = InserterModule<PairedReadInserterContainer<DBCon>>;

} // namespace libMSV


#ifdef WITH_PYTHON
/// @brief expose libMSV::ReadInsert to python
void exportReadInserter( libMS::SubmoduleOrganizer& xOrganizer );
#endif
