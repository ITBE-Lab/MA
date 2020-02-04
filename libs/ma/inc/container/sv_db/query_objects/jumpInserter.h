/**
 * @file jumpInserter.h
 * @brief implements libMA::SvJumpInserter that inserts libMA::SvJump objects into the DB.
 * @author Markus Schmidt
 */
#pragma once

#include "container/sv_db/svSchema.h"
#include "container/sv_db/tables/nameDesc.h"
#include "module/get_inserter_container_module.h"

namespace libMA
{
/**
 * @brief Insertion of libMA::SvJump into the database.
 */
template <typename DBCon>
class JumpInserterContainer : public InserterContainer<DBCon, ReadTable, ContainerVector<SvJump>, NucSeq>
{
    /**
     * @brief Inner helper class. Memorizes the context for a single read.
     * @details Context: Read ID, ID of current caller run, holds the table via shared pointer
     */
    class ReadContex
    {
      private:
        INSERTER_TYPE pInserter;
        const int64_t iSvJumpRunId;
        const int64_t iReadId;

      public:
        /// @brief create the context for the jump run iSvJumpRunId and read iReadId.
        ReadContex( INSERTER_TYPE pInserter, const int64_t iSvJumpRunId, const int64_t iReadId )
            : pInserter( pInserter ), iSvJumpRunId( iSvJumpRunId ), iReadId( iReadId )
        {} // constructor

        /** @brief Insert the jump rJump into the DB.
         * @details
         * expects rJump not to be in the DB.
         * rJump.iId will be assigned the new id of rJump
         * rJump.iReadId must must match the read this contig was created for.
         * if rJump.iReadId == -1 the jump's read id will be overritten with the correct read id.
         */
        inline void insertJump( SvJump& rJump )
        {
            // make sure the read id matches the read context
            if( rJump.iReadId == -1 ) // if there is no read id given yet add it
                rJump.iReadId = iReadId;
            else // otherwise assert it matches
                assert( rJump.iReadId == iReadId );

            if( rJump.does_switch_strand( ) )
                assert( rJump.from_start( ) >= std::numeric_limits<int64_t>::max( ) / 2 );
            rJump.iId = pInserter->insert( nullptr, iSvJumpRunId, rJump.iReadId, rJump.from_start( ), rJump.from_end( ),
                                           (uint32_t)rJump.uiFrom, (uint32_t)rJump.uiTo, (uint32_t)rJump.uiQueryFrom,
                                           (uint32_t)rJump.uiQueryTo, (uint32_t)rJump.uiNumSupportingNt,
                                           rJump.bFromForward, rJump.bToForward, rJump.bFromSeedStart );
        } // method
    }; // class
  public:
    /**
     * @brief Open a context for the read with id = iReadId, which can be later used for inserting jumps.
     */
    inline ReadContex readContext( int64_t iReadId )
    {
        return ReadContex( pInserter, iId, iReadId );
    } // method

    virtual void insert( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        auto xReadContext = readContext( pRead->iId );
        for( SvJump& rJump : *pJumps )
            xReadContext.insertJump( rJump ); // also updates the jump ids;
    } // method
}; // class

template <typename DBCon> class JumpInserterContainer : public Container
{
    const static size_t BUFFER_SIZE = 500;
    using INSERTER_TYPE = decltype( SvJumpTable<DBCon>::template getBulkInserter<BUFFER_SIZE> );
    INSERTER_TYPE pInserter;

  public:
    /// @brief the id of the run this inserter is attached to.
    const int64_t iSvJumpRunId;


    /**
     * @brief creates a jump inserter for the run with id = iSvJumpRunId
     * @param pDB the sv database
     * @param iSvJumpRunId caller run id
     */
    JumpInserterContainer( std::shared_ptr<DBCon> pConnection, int64_t iSvJumpRunId )
        : pInserter( SvJumpTable<DBCon>( pConnection ).template getBulkInserter<BUFFER_SIZE>( ) ),
          iSvJumpRunId( iSvJumpRunId )
    {} // constructor

    /**
     * @brief Open a context for the read with id = iReadId, which can be later used for inserting jumps.
     */
    inline ReadContex readContext( int64_t iReadId )
    {
        return ReadContex( pInserter, iSvJumpRunId, iReadId );
    } // method

    inline void insert( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        auto xReadContext = readContext( pRead->iId );
        for( SvJump& rJump : *pJumps )
            xReadContext.insertJump( rJump ); // also updates the jump ids;
    } // method

    /// @brief closes the inserter (they cannot be reopened again)
    /// @details intended for python as replacement for the deconstructor
    inline void close( )
    {
        pInserter->reset( );
    } // method
}; // class

template <typename DBCon, typename DBConInit>
using GetJumpInserterContainerModule =
    GetInserterContainerModule<JumpInserterContainer, DBCon, DBConInit, NameDescTable<"sv_jump_run_table">>;


template <typename DBCon>
using JumpInserterModule = InserterModule<JumpInserterContainer<DBCon>>;


} // namespace libMA

#ifdef WITH_PYTHON
/// @brief used to expose libMA::SvJumpInserter to python
void exportSvJumpInserter( py::module& rxPyModuleId );
#endif
