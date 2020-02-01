/**
 * @file jumpInserter.h
 * @brief implements libMA::SvJumpInserter that inserts libMA::SvJump objects into the DB.
 * @author Markus Schmidt
 */
#pragma once

#include "container/sv_db/svDb.h"
#include "db_config.h"

namespace libMA
{
/**
 * @brief Idea: Insertion of libMA::SvJump into the database.
 */
template <typename DBCon> class SvJumpInserter
{
  public:
    // this is here so that it gets destructed after the transaction context
    std::shared_ptr<_SV_DB<DBCon>> pDB;

  private:
    // must be after the DB so that it is deconstructed first


    // REPLACED: std::shared_ptr<CppSQLiteExtImmediateTransactionContext> pTransactionContext;
    typename DBCon::sharedGuardedTrxnType pGuardedTrxn; // technically a shared pointer

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
        std::shared_ptr<SvJumpTable<DBCon>> pSvJumpTable;
        const int64_t iSvJumpRunId;
        const int64_t iReadId;

      public:
        /// @brief create the context for the jump run iSvJumpRunId and read iReadId.
        ReadContex( std::shared_ptr<SvJumpTable<DBCon>> pSvJumpTable, const int64_t iSvJumpRunId,
                    const int64_t iReadId )
            : pSvJumpTable( pSvJumpTable ), iSvJumpRunId( iSvJumpRunId ), iReadId( iReadId )
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
            // DEL: rJump.iId = pSvJumpTable->xInsertRow(
            rJump.iId = pSvJumpTable->insert( iSvJumpRunId, rJump.iReadId, rJump.from_start( ), rJump.from_end( ),
                                              (uint32_t)rJump.uiFrom, (uint32_t)rJump.uiTo, (uint32_t)rJump.uiQueryFrom,
                                              (uint32_t)rJump.uiQueryTo, (uint32_t)rJump.uiNumSupportingNt,
                                              rJump.bFromForward, rJump.bToForward, rJump.bFromSeedStart );
        } // method
    }; // class

    /**
     * @brief creates a jump inserter for the run with id = iSvJumpRunId
     * @param pDB the sv database
     * @param iSvJumpRunId caller run id
     * @param bTransactionLess do not create a transaction
     */
    SvJumpInserter( std::shared_ptr<_SV_DB<DBCon>> pDB, int64_t iSvJumpRunId, bool bTransactionLess )
        : pDB( std::make_shared<_SV_DB<DBCon>>(
              *pDB ) ), // create a copy of the connection to the db
                        // REPLACED: pTransactionContext(bTransactionLess
                        // REPLACED: 	? nullptr
                        // REPLACED: 	: std::make_shared<CppSQLiteExtImmediateTransactionContext>(*pDB->pDatabase)),
          pGuardedTrxn( bTransactionLess ? nullptr : pDB->pDatabase->sharedGuardedTrxn( ) ),

          iSvJumpRunId( iSvJumpRunId )
    {} // constructor

    /**
     * @brief creates a jump inserter for the run with id = iSvJumpRunId
     * @param pDB the sv database
     * @param iSvJumpRunId caller run id
     */
    SvJumpInserter( std::shared_ptr<_SV_DB<DBCon>> pDB, int64_t iSvJumpRunId )
        : SvJumpInserter( pDB, iSvJumpRunId, false )
    {} // constructor

    /**
     * @brief creates a jump inserter for a new jump-run.
     * @details
     * The new run gets a name and description
     */
    SvJumpInserter( std::shared_ptr<_SV_DB<DBCon>> pDB, const std::string& rsSvCallerName,
                    const std::string& rsSvCallerDesc )
        : pDB( pDB ),
          // REPLACED: pTransactionContext(std::make_shared<CppSQLiteExtImmediateTransactionContext>(*pDB->pDatabase)),
          pGuardedTrxn( pDB->pDatabase->sharedGuardedTrxn( ) ),
          iSvJumpRunId( pDB->pSvJumpRunTable->insert( rsSvCallerName, rsSvCallerDesc ) )
    {} // constructor

    /**
     * @brief Open a context for the read with id = iReadId, which can be later used for inserting jumps.
     */
    inline ReadContex readContext( int64_t iReadId )
    {
        return ReadContex( pDB->pSvJumpTable, iSvJumpRunId, iReadId );
    } // method

    /**
     * @brief terminates the current transaction
     */
    inline void endTransaction( )
    {
        pGuardedTrxn.reset( );
    } // method

    /**
     * @brief terminates the current transaction and starts a new one
     */
    inline void reOpenTransaction( )
    {
        endTransaction( );
        pGuardedTrxn = pDB->pDatabase->sharedGuardedTrxn( );
    }; // method
}; // class


/**
 * @brief Wraps a jump inserter, so that it can become part of a computational graph.
 */
template <typename DBCon> class SvDbInserter : public Module<Container, false, ContainerVector<SvJump>, NucSeq>
{
    std::shared_ptr<_SV_DB<DBCon>> pDb;

  public:
    /// @brief the jump inserter; This creates a transaction
    SvJumpInserter<DBCon> xInserter;

    ///@brief creates a new jump-run with the name MA-SV.
    SvDbInserter( const ParameterSetManager& rParameters, std::shared_ptr<_SV_DB<DBCon>> pDb, std::string sRunDesc )
        : pDb( pDb ), xInserter( this->pDb, "MA-SV", sRunDesc )
    {} // constructor

    /// @brief insert all jumps in pJumps for the read pRead
    std::shared_ptr<Container> execute( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        std::lock_guard<std::mutex> xGuard( *pDb->pWriteLock );

        typename SvJumpInserter<DBCon>::ReadContex xReadContext = xInserter.readContext( pRead->iId );
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
template <typename DBCon> class BufferedSvDbInserter : public Module<Container, false, ContainerVector<SvJump>, NucSeq>
{
    std::shared_ptr<SvJumpInserter<DBCon>> pInserter;

  public:
    /// @brief the buffer containing all jumps that are not yet commited
    std::vector<std::pair<std::shared_ptr<ContainerVector<SvJump>>, int64_t>> vBuffer;


    /**
     * @brief create the inserter for the run with id iSvJumpRunId
     * @details
     */
    BufferedSvDbInserter( const ParameterSetManager& rParameters, std::shared_ptr<SvJumpInserter<DBCon>> pInserter )
        : pInserter( pInserter )
    {} // constructor

    /// @brief bulk insert all jumps in the buffer; then clear the buffer
    // Optimization for SQLite only ...
    inline void commit( bool bForce = false )
    {
        if( bForce == false && vBuffer.size( ) < 10000 )
            return;
        if( vBuffer.size( ) == 0 )
            return;

        {
            std::lock_guard<std::mutex> xGuard( *pInserter->pDB->pWriteLock );
            for( auto xPair : vBuffer )
            {
                typename SvJumpInserter<DBCon>::ReadContex xReadContext = pInserter->readContext( xPair.second );
                for( SvJump& rJump : *xPair.first )
                    xReadContext.insertJump( rJump ); // also updates the jump ids;
                pInserter->reOpenTransaction( );
            } // for
        } // scope for xGuard
        vBuffer.clear( );

        // end of scope for lock guard
    } // method

    /// @brief triggers commit
    ~BufferedSvDbInserter( )
    {
        commit( true );
    } // destructor

    /// @brief buffer all jumps in pJumps for the read pRead
    std::shared_ptr<Container> execute( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        vBuffer.emplace_back( pJumps, pRead->iId );
        commit( );
        return std::make_shared<Container>( );
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/// @brief used to expose libMA::SvJumpInserter to python
void exportSvJumpInserter( py::module& rxPyModuleId );
#endif
