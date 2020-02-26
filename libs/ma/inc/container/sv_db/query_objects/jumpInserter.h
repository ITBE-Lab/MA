/**
 * @file jumpInserter.h
 * @brief implements libMA::SvJumpInserter that inserts libMA::SvJump objects into the DB.
 * @author Markus Schmidt
 */
#pragma once

#include "container/svJump.h"
#include "container/sv_db/tables/nameDesc.h"
#include "container/sv_db/tables/svJump.h"
#include "module/get_inserter_container_module.h"

namespace libMA
{
/**
 * @brief Insertion of libMA::SvJump into the database.
 */
template <typename DBCon>
class JumpInserterContainer
    : public BulkInserterContainer<DBCon, AbstractInserterContainer, SvJumpTable, ContainerVector<SvJump>, NucSeq>
{
  public:
    using ParentType =
        BulkInserterContainer<DBCon, AbstractInserterContainer, SvJumpTable, ContainerVector<SvJump>, NucSeq>;
    using ParentType::BulkInserterContainer;

  protected:
    virtual size_t EXPORTED
    insert_override( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
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

            if( rJump.does_switch_strand( ) )
                assert( rJump.from_start( ) >= std::numeric_limits<int64_t>::max( ) / 2 );
            rJump.iId = ParentType::pInserter->insert(
                ParentType::iId /* <- id of the sv jump run */, rJump.iReadId, rJump.from_start( ), rJump.from_end( ),
                (uint32_t)rJump.uiFrom, (uint32_t)rJump.uiTo, (uint32_t)rJump.uiQueryFrom, (uint32_t)rJump.uiQueryTo,
                (uint32_t)rJump.uiNumSupportingNt, rJump.bFromForward, rJump.bToForward, rJump.bFromSeedStart );
        } // for
        return pJumps->size( );
    } // method
}; // class

template <typename DBCon, typename DBConInit>
using GetJumpInserterContainerModule =
    GetInserterContainerModule<JumpInserterContainer, DBCon, DBConInit, SvJumpRunTable>;


template <typename DBCon> using JumpInserterModule = InserterModule<JumpInserterContainer<DBCon>>;


} // namespace libMA

#ifdef WITH_PYTHON
/// @brief used to expose libMA::SvJumpInserter to python
void exportSvJumpInserter( py::module& rxPyModuleId );
#endif
