/**
 * @file kMerInserter.h
 * @brief implements libMSV::KMerInserter that inserts k-mers objects into the DB.
 * @author Markus Schmidt
 */
#pragma once

#include "ms/module/get_inserter_container_module.h"
#include "msv/container/sv_db/tables/kMerFilter.h"
#include "msv/container/sv_db/tables/sequencer.h"
#include "msv/module/count_k_mers.h"

using namespace libMA;
using namespace libMS;

namespace libMSV
{
/**
 * @brief Insertion of libMSV::SvJump into the database.
 */
template <typename DBCon>
class KMerInserterContainer : public BulkInserterContainer<DBCon, AbstractInserterContainer, KMerFilterTable, NucSeq>
{
  public:
    using ParentType = BulkInserterContainer<DBCon, libMS::AbstractInserterContainer, KMerFilterTable, NucSeq>;
    using ParentType::BulkInserterContainer;

  protected:
    virtual size_t insert_override( std::shared_ptr<NucSeq> pRead )
    {
        size_t uiCnt = 0;
        KMerCounter::toKMers( pRead, 0, pRead->length( ), 1, [&]( const NucSeq& xNucSeq ) {
            uiCnt++;
            auto pInsert = std::make_shared<NucSeq>( xNucSeq );
            ParentType::pInserter->insert( ParentType::iId, NucSeqSql(pInsert), (uint32_t)1 );
            return true;
        } );
        return uiCnt;
    } // method
}; // class

template <typename DBCon, typename DBConInit>
using GetKMerInserterContainerModule =
    GetInserterContainerModule<KMerInserterContainer, DBCon, DBConInit, SequencerTable>;


template <typename DBCon> using KMerInserterModule = InserterModule<KMerInserterContainer<DBCon>>;

} // namespace libMSV

#ifdef WITH_PYTHON
/// @brief used to expose libMSV::KMerInserter to python
void exportKMerInserter( libMS::SubmoduleOrganizer& xOrganizer );
#endif
