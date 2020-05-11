/**
 * @file jumpInserter.h
 * @brief implements libMSV::SvJumpInserter that inserts libMSV::SvJump objects into the DB.
 * @author Markus Schmidt
 */
#pragma once

#include "ms/module/get_inserter_container_module.h"
#include "msv/container/svJump.h"
#include "msv/container/sv_db/tables/nameDesc.h"
#include "msv/container/sv_db/tables/svJump.h"
#include "util/geom.h"

using namespace libMA;
using namespace libMS;

namespace libMSV
{
/**
 * @brief Insertion of libMSV::SvJump into the database.
 */
template <typename DBCon>
class JumpInserterContainer
    : public BulkInserterContainer<DBCon, AbstractInserterContainer, SvJumpTable, ContainerVector<SvJump>, NucSeq>
{
  public:
    using ParentType =
        BulkInserterContainer<DBCon, libMS::AbstractInserterContainer, SvJumpTable, ContainerVector<SvJump>, NucSeq>;
    using ParentType::BulkInserterContainer;

  protected:
    virtual size_t insert_override( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        const int64_t iReadId = pRead->iId;
        for( SvJump& rJump : *pJumps )
        {
            // make sure the read id matches the read context
            if( rJump.iReadId == -1 ) // if there is no read id given yet add it
                rJump.iReadId = iReadId;
            else // otherwise assert it matches
                assert( rJump.iReadId == iReadId );

            assert( rJump.iReadId != -1 );

            auto xRectangle = WKBUint64Rectangle( geom::Rectangle<nucSeqIndex>(
                rJump.from_start_same_strand( ), rJump.to_start( ), rJump.from_size( ), rJump.to_size( ) ) );
            rJump.iId = ParentType::pInserter->insert(
                ParentType::iId /* <- id of the sv jump run */, rJump.iReadId, rJump.from_start( ), rJump.from_end( ),
                (uint32_t)rJump.uiFrom, (uint32_t)rJump.uiTo, (uint32_t)rJump.uiQueryFrom, (uint32_t)rJump.uiQueryTo,
                (uint32_t)rJump.uiNumSupportingNt, rJump.bFromForward, rJump.bToForward, rJump.bFromSeedStart,
                xRectangle );
        } // for
        return pJumps->size( );
    } // method
}; // class

template <typename DBCon, typename DBConInit>
using GetJumpInserterContainerModule =
    GetInserterContainerModule<JumpInserterContainer, DBCon, DBConInit, SvJumpRunTable>;


template <typename DBCon> using JumpInserterModule = InserterModule<JumpInserterContainer<DBCon>>;


} // namespace libMSV

#ifdef WITH_PYTHON
/// @brief used to expose libMSV::SvJumpInserter to python
void exportSvJumpInserter( libMS::SubmoduleOrganizer& xOrganizer );
#endif
