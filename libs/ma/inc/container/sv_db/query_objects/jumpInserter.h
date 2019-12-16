/**
 * @file jumpInserter.h
 * @brief implements libMA::SvJumpInserter that inserts libMA::SvJump objects into the DB.
 * @author Markus Schmidt
 */
#include "container/sv_db/svDb.h"

#pragma once

namespace libMA
{

/**
 * @brief Idea: Insertion of libMA::SvJump into the database.
 */
class SvJumpInserter
{
    // this is here so that it gets destructed after the transaction context
    std::shared_ptr<SV_DB> pDB;
    // must be after the DB so that it is deconstructed first
    std::shared_ptr<CppSQLiteExtImmediateTransactionContext> pTransactionContext;

  public:
    /// @brief the id of the run this inserter is attached to.
    const int64_t iSvJumpRunId;

    /**
     * @brief Inner helper class. Memorizes the context for a single read.
     * @details Context: Read ID, ID of current caller run, holds the table via shared pointer
     */
    class ReadContex
    {
      private:
        std::shared_ptr<SvJumpTable> pSvJumpTable;
        const int64_t iSvJumpRunId;
        const int64_t iReadId;

      public:
        /// @brief create the context for the jump run iSvJumpRunId and read iReadId.
        ReadContex( std::shared_ptr<SvJumpTable> pSvJumpTable, const int64_t iSvJumpRunId, const int64_t iReadId )
            : pSvJumpTable( pSvJumpTable ), iSvJumpRunId( iSvJumpRunId ), iReadId( iReadId )
        {} // constructor

        /** @brief insert the jump rJump into the DB.
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
            rJump.iId = pSvJumpTable->xInsertRow(
                iSvJumpRunId, rJump.iReadId, rJump.from_start( ), rJump.from_end( ), (uint32_t)rJump.uiFrom,
                (uint32_t)rJump.uiTo, (uint32_t)rJump.uiQueryFrom, (uint32_t)rJump.uiQueryTo,
                (uint32_t)rJump.uiNumSupportingNt, rJump.bFromForward, rJump.bToForward, rJump.bFromSeedStart );
        } // method
    }; // class

    /**
     * @brief creates a jump inserter for the run with id = iSvJumpRunId
     * @param pDB the sv database
     * @param iSvJumpRunId caller run id
     * @param bTransactionLess do not create a transaction
     */
    SvJumpInserter( std::shared_ptr<SV_DB> pDB, int64_t iSvJumpRunId, bool bTransactionLess )
        : pDB( pDB ),
          pTransactionContext( bTransactionLess
                                   ? nullptr
                                   : std::make_shared<CppSQLiteExtImmediateTransactionContext>( *pDB->pDatabase ) ),
          iSvJumpRunId( iSvJumpRunId )
    {} // constructor

    /**
     * @brief creates a jump inserter for the run with id = iSvJumpRunId
     * @param pDB the sv database
     * @param iSvJumpRunId caller run id
     */
    SvJumpInserter( std::shared_ptr<SV_DB> pDB, int64_t iSvJumpRunId ) : SvJumpInserter( pDB, iSvJumpRunId, false )
    {} // constructor

    /**
     * @brief creates a jump inserter for a new jump-run.
     * @details
     * The new run gets a name and description
     */
    SvJumpInserter( std::shared_ptr<SV_DB> pDB, const std::string& rsSvCallerName, const std::string& rsSvCallerDesc )
        : pDB( pDB ),
          pTransactionContext( std::make_shared<CppSQLiteExtImmediateTransactionContext>( *pDB->pDatabase ) ),
          iSvJumpRunId( pDB->pSvJumpRunTable->insert( rsSvCallerName, rsSvCallerDesc ) )
    {} // constructor

    /**
     * @brief Open a context for the read with id = iReadId, which can be later used for inserting jumps.
     */
    inline ReadContex readContext( int64_t iReadId )
    {
        return ReadContex( pDB->pSvJumpTable, iSvJumpRunId, iReadId );
    } // method

}; // class


/**
 * @brief Wraps a jump inserter, so that it can become part of a computational graph.
 */
class SvDbInserter : public Module<Container, false, ContainerVector<SvJump>, NucSeq>
{
    std::shared_ptr<SV_DB> pDb;

  public:
    /// @brief the jump inserter; This creates a transaction
    SvJumpInserter xInserter;

    ///@brief creates a new jump-run with the name MA-SV.
    SvDbInserter( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, std::string sRunDesc )
        : pDb( pDb ), xInserter( this->pDb, "MA-SV", sRunDesc )
    {} // constructor

    /// @brief insert all jumps in pJumps for the read pRead
    std::shared_ptr<Container> execute( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        std::lock_guard<std::mutex> xGuard( *pDb->pWriteLock );

        SvJumpInserter::ReadContex xReadContext = xInserter.readContext( pRead->iId );
        for( SvJump& rJump : *pJumps )
            xReadContext.insertJump( rJump ); // also updates the jump ids;

        return std::make_shared<Container>( );
        // end of score for xGuard
    } // method
}; // class

/**
 * @brief Wraps a jump inserter, so that it can become part of a computational graph.
 * @details
 * buffers all jumps in a vector and then bulk inserts them once commit is called.
 * The destructor calls commit automatically.
 * If the buffer holds 10.000 elements a bulk insert is triggered as well.
 */
class BufferedSvDbInserter : public Module<Container, false, ContainerVector<SvJump>, NucSeq>
{
    std::shared_ptr<SV_DB> pDb;
    int64_t iSvJumpRunId;

  public:
    /// @brief the buffer containing all jumps that are not yet commited
    std::vector<std::pair<std::shared_ptr<ContainerVector<SvJump>>, int64_t>> vBuffer;


    /**
     * @brief create the inserter for the run with id iSvJumpRunId
     * @details
     * the run with id iSvJumpRunId must exist.
     */
    BufferedSvDbInserter( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvJumpRunId )
        : pDb( pDb ), iSvJumpRunId( iSvJumpRunId )
    {} // constructor

    /// @brief bulk insert all jumps in the buffer; then clear the buffer
    inline void commit( bool bTransactionLess = true, bool bForce = false )
    {
        if( !bForce && vBuffer.size( ) < 10000 )
            return;
        if( vBuffer.size( ) == 0 )
            return;
        SvJumpInserter xInserter( pDb, iSvJumpRunId );
        std::lock_guard<std::mutex> xGuard( *pDb->pWriteLock );
        for( auto xPair : vBuffer )
        {
            SvJumpInserter::ReadContex xReadContext = xInserter.readContext( xPair.second );
            for( SvJump& rJump : *xPair.first )
                xReadContext.insertJump( rJump ); // also updates the jump ids;
        } // for
        vBuffer.clear( );
        // end of scope for lock guard
    } // method

    /// @brief triggers commit
    ~BufferedSvDbInserter( )
    {
        commit( false, true );
    } // destructor

    /// @brief buffer all jumps in pJumps for the read pRead
    std::shared_ptr<Container> execute( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        vBuffer.emplace_back( pJumps, pRead->iId );
        commit();
        return std::make_shared<Container>( );
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/// @brief used to expose libMA::SvJumpInserter to python
void exportSvJumpInserter( py::module& rxPyModuleId );
#endif