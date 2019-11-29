#include "container/sv_db/svDb.h"

#pragma once

namespace libMA
{

/**@brief Idea: Insertion of Jumps (...) into the database.
 */
class SvJumpInserter
{
    // this is here so that it gets destructed after the transaction context
    std::shared_ptr<SV_DB> pDB;
    // must be after the DB so that it is deconstructed first
    CppSQLiteExtImmediateTransactionContext xTransactionContext;

  public:
    const int64_t iSvJumpRunId;

    /**@brief Inner helper class. Memories the coontext for a single read.
     * @details COntext: Read ID, ID of current caller run, holds the table via shared pointer 
     */
    class ReadContex
    {
      private:
        std::shared_ptr<SvJumpTable> pSvJumpTable;
        const int64_t iSvJumpRunId;
        const int64_t iReadId;

      public:
        ReadContex( std::shared_ptr<SvJumpTable> pSvJumpTable, const int64_t iSvJumpRunId, const int64_t iReadId )
            : pSvJumpTable( pSvJumpTable ), iSvJumpRunId( iSvJumpRunId ), iReadId( iReadId )
        {} // constructor

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
     * @brief Constructor 
     * @param pDB database
     * @param iSvJumpRunId caller run id
     */
    SvJumpInserter( std::shared_ptr<SV_DB> pDB, int64_t iSvJumpRunId )
        : pDB( pDB ), xTransactionContext( *pDB->pDatabase ), iSvJumpRunId( iSvJumpRunId )
    {} // constructor

    /**@brief Constructor 
     * @param database
     * @detail Extracts the caller run id from databse
     */
    SvJumpInserter( std::shared_ptr<SV_DB> pDB, const std::string& rsSvCallerName, const std::string& rsSvCallerDesc )
        : pDB( pDB ),
          xTransactionContext( *pDB->pDatabase ),
          iSvJumpRunId( pDB->pSvJumpRunTable->insert( rsSvCallerName, rsSvCallerDesc ) )
    {} // constructor

    
    /**@brief OPen ths context for the read, which can be later used for inserting jumps.
     */
    inline ReadContex readContext( int64_t iReadId )
    {
        return ReadContex( pDB->pSvJumpTable, iSvJumpRunId, iReadId );
    } // method

}; // class


/**@brief Wraps a jump inserter, so that it can become part of a computational graph.
 */
class SvDbInserter : public Module<Container, false, ContainerVector<SvJump>, NucSeq>
{
    std::shared_ptr<SV_DB> pDb;

  public:
    // this creates a transaction
    SvJumpInserter xInserter;

    SvDbInserter( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, std::string sRunDesc )
        : pDb( pDb ), xInserter( this->pDb, "MA-SV", sRunDesc )
    {} // constructor

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

class BufferedSvDbInserter : public Module<Container, false, ContainerVector<SvJump>, NucSeq>
{
    std::shared_ptr<SV_DB> pDb;
    int64_t iSvJumpRunId;

  public:
    std::vector<std::pair<std::shared_ptr<ContainerVector<SvJump>>, int64_t>> vBuffer;
    // this creates a transaction

    BufferedSvDbInserter( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSvJumpRunId )
        : pDb( pDb ), iSvJumpRunId( iSvJumpRunId )
    {} // constructor

    inline void commit( )
    {
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

    ~BufferedSvDbInserter( )
    {
        commit( );
    } // destructor

    std::shared_ptr<Container> execute( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        vBuffer.emplace_back( pJumps, pRead->iId );
        return std::make_shared<Container>( );
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
void exportSvJumpInserter( py::module& rxPyModuleId );
#endif