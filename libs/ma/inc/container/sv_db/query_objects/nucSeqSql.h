/**
 * @file nucSeqSql.h
 * @brief implements libMA::AllNucSeqFromSql that fetches reads from the DB.
 * @author Markus Schmidt
 */
#pragma once

#include "container/nucSeq.h"
#include "container/sv_db/pool_container.h"
#include "container/sv_db/tables/pairedRead.h"
#include "container/sv_db/tables/read.h"
#include "module/module.h"
#include "sql_api.h"

namespace libMA
{

template <typename DBCon> class NucSeqQueryContainer : public Container
{
  public:
    int iConnectionId;
    int64_t iSequencerId;
    uint32_t uiRes;
    uint32_t uiModulo;
    PriKeyDefaultType uiMinId;
    SQLQuery<DBCon, std::shared_ptr<CompressedNucSeq>, PriKeyDefaultType> xQuery;

    /** @brief actually execute the query
     * @details
     * eventhough pConnection is not needed in the function, the connection must be held, so that concurrent
     * accesses on one connection are avoided.
     */
    inline void exec( std::shared_ptr<DBCon> pConnection )
    {
        if( iSequencerId != -1 && uiModulo != 1 )
            xQuery.execAndFetch( uiMinId, iSequencerId, uiModulo, uiRes );
        else if( iSequencerId != -1 )
            xQuery.execAndFetch( uiMinId, iSequencerId );
        else if( uiModulo != 1 )
            xQuery.execAndFetch( uiMinId, uiModulo, uiRes );
        else
            xQuery.execAndFetch( uiMinId );
    } // method

    NucSeqQueryContainer( int iConnectionId, std::shared_ptr<DBCon> pConnection, int64_t iSequencerId, uint32_t uiRes,
                          uint32_t uiModulo, bool bPaired, bool bUnpaired )
        : iConnectionId( iConnectionId ),
          iSequencerId( iSequencerId ),
          uiRes( uiRes ),
          uiModulo( uiModulo ),
          uiMinId( 0 ),
          xQuery( pConnection,
                  std::string( "SELECT sequence, id "
                               "FROM read_table "
                               "WHERE id >= ? " )
                      .append( iSequencerId != -1 ? "AND sequencer_id = ? " : "" )
                      .append( uiModulo != 1 ? "AND id % ? = ? " : "" )
                      .append( bPaired ? ""
                                       : "AND id NOT IN ( "
                                         "   SELECT first_read FROM paired_read_table "
                                         "   UNION "
                                         "   SELECT second_read FROM paired_read_table "
                                         ") " )
                      .append( bUnpaired ? ""
                                         : "AND id IN ( "
                                           "   SELECT first_read FROM paired_read_table "
                                           "   UNION "
                                           "   SELECT second_read FROM paired_read_table "
                                           ") " )
                      .append( "ORDER BY id ASC "
                               "LIMIT 1000 " ) )
    {
        exec( pConnection );
    } // constructor

    inline void next( std::shared_ptr<PoolContainer<DBCon>> pPool )
    {
        if( !xQuery.next( ) )
            pPool->xPool.run( iConnectionId, [&]( auto pConnection ) { this->exec( pConnection ); } );
    } // method

    inline bool eof( )
    {
        return xQuery.eof( );
    } // method

    inline std::tuple<std::shared_ptr<CompressedNucSeq>, PriKeyDefaultType> get( )
    {
        auto xRet = xQuery.get( );
        uiMinId = std::get<1>( xRet ) + 1; // inc min id
        return xRet;
    } // method
}; // class

template <typename DBCon>
class GetNucSeqFromSqlQuery : public Module<NucSeqQueryContainer<DBCon>, false, PoolContainer<DBCon>>
{
    int64_t iSequencerId;
    uint32_t uiRes;
    uint32_t uiModulo;
    bool bPaired;
    bool bUnpaired;

  public:
    GetNucSeqFromSqlQuery( const ParameterSetManager& rParameters, int64_t iSequencerId, size_t uiRes, size_t uiModulo,
                           bool bPaired, bool bUnpaired )
        : iSequencerId( iSequencerId ),
          uiRes( (uint32_t)uiRes ),
          uiModulo( (uint32_t)uiModulo ),
          bPaired( bPaired ),
          bUnpaired( bUnpaired )
    {
        if( !bPaired && !bUnpaired )
            throw std::runtime_error( "NucSeqFromSqlQuery is not fetching anything" );
    } // constructor

    /// @brief returns a query that can fetch NucSeqs.
    virtual std::shared_ptr<NucSeqQueryContainer<DBCon>> EXPORTED execute( std::shared_ptr<PoolContainer<DBCon>> pPool )
    {
        int iConnectionId = pPool->xPool.getDedicatedConId( );
        return pPool->xPool.run( iConnectionId, [&]( auto pConnection ) {
            return std::make_shared<NucSeqQueryContainer<DBCon>>( iConnectionId, pConnection, iSequencerId, uiRes,
                                                                  uiModulo, bPaired, bUnpaired );
        } );
    } // method
}; // class


/**
 * @brief fetches reads from a database
 */
template <typename DBCon>
class NucSeqFetcher : public Module<NucSeq, true, PoolContainer<DBCon>, NucSeqQueryContainer<DBCon>>
{
  public:
    NucSeqFetcher( const ParameterSetManager& rParameters )
    {}

    /// @brief returns one read at a time until isFinished returns true.
    virtual std::shared_ptr<NucSeq> EXPORTED execute( std::shared_ptr<PoolContainer<DBCon>> pPool,
                                                      std::shared_ptr<NucSeqQueryContainer<DBCon>> pQuery )
    {
        if( pQuery->eof( ) )
            return nullptr;

        auto xTup = pQuery->get( );
        auto pRet = std::get<0>( xTup )->pUncomNucSeq;
        pRet->iId = (int64_t)std::get<1>( xTup );
        assert( pRet->iId != -1 );
        pQuery->next( pPool ); // increment the iterator
        return pRet;
    } // method
}; // class

} // namespace libMA


#ifdef WITH_PYTHON
/// @brief expose the NucSeqFromSql classes to python
void exportNucSeqSql( py::module& rxPyModuleId );
#endif
