/**
 * @file nucSeqSql.h
 * @brief implements libMA::AllNucSeqFromSql that fetches reads from the DB.
 * @author Markus Schmidt
 */
#pragma once

#include "db_base.h"
#include "sql_api.h"
#include "container/nucSeq.h"
#include "container/sv_db/pool_container.h"
#include "container/sv_db/tables/pairedRead.h"
#include "container/sv_db/tables/read.h"
#include "module/module.h"

namespace libMA
{

template <typename DBCon>
class NucSeqQueryContainer : public SQLQuery<DBCon, std::shared_ptr<CompressedNucSeq>, PriKeyDefaultType>,
                             public Container
{
  public:
    int iConnectionId;

    NucSeqQueryContainer( int iConnectionId, std::shared_ptr<DBCon> pConnection, bool bDoSequencerId, bool bDoModulo,
                          bool bUnpairedOnly )
        : SQLQuery<DBCon, std::shared_ptr<CompressedNucSeq>, PriKeyDefaultType>(
              pConnection,
              std::string( "SELECT read_table.sequence, read_table.id "
                           "FROM read_table " )
                  .append( bDoSequencerId ? "WHERE sequencer_id = ? " : "" )
                  .append( bDoModulo ? "AND read_table.id % ? = ? " : "" )
                  .append( bUnpairedOnly && bDoSequencerId ? "AND " : bUnpairedOnly ? "WHERE " : "" )
                  .append( bUnpairedOnly ? "read_table.id NOT IN ( "
                                           "   SELECT paired_read_table.first_read FROM paired_read_table "
                                           "   UNION "
                                           "   SELECT paired_read_table.second_read FROM paired_read_table "
                                           ") "
                                         : "" ) ),
          iConnectionId( iConnectionId )
    {
        // we do not allow bDoModulo=true while bDoSequencerId=false
        assert( !bDoModulo || bDoSequencerId );
    } // constructor
}; // class

template <typename DBCon>
class GetNucSeqFromSqlQuery : public Module<NucSeqQueryContainer<DBCon>, false, PoolContainer<DBCon>>
{
    int64_t iSequencerId;
    uint32_t uiRes;
    uint32_t uiModulo;
    bool bAll;

  public:
    GetNucSeqFromSqlQuery( const ParameterSetManager& rParameters, int64_t iSequencerId, size_t uiRes, size_t uiModulo,
                           bool bPaired, bool bUnpaired )
        : iSequencerId( iSequencerId ),
          uiRes( (uint32_t)uiRes ),
          uiModulo( (uint32_t)uiModulo ),
          bAll( bPaired && bUnpaired )
    {
        if( !bPaired && !bUnpaired )
            throw std::runtime_error( "NucSeqFromSqlQuery is not fetching anything" );
        if( !bAll && bPaired )
            throw std::runtime_error( "NucSeqFromSqlQuery: fetching paired NucSeqs only is unimplemented" );
    } // constructor

    /// @brief returns a query that can fetch NucSeqs.
    virtual std::shared_ptr<NucSeqQueryContainer<DBCon>> EXPORTED execute( std::shared_ptr<PoolContainer<DBCon>> pPool )
    {
        int iConnectionId = pPool->xPool.getDedicatedConId( );
        return pPool->xPool.run(
            iConnectionId,
            []( auto pConnection,
                int iConnectionId,
                int64_t iSequencerId,
                uint32_t uiRes,
                uint32_t uiModulo,
                bool bAll ) //
            {
                auto pQuery = std::make_shared<NucSeqQueryContainer<DBCon>>(
                    iConnectionId, pConnection, iSequencerId != -1, iSequencerId != -1 && uiModulo != 1, !bAll );

                if( iSequencerId != -1 && uiModulo != 1 )
                    pQuery->execAndFetch( iSequencerId, uiModulo, uiRes );
                else if( iSequencerId != -1 )
                    pQuery->execAndFetch( iSequencerId );
                else
                    pQuery->execAndFetch( );

                /* @note: this expects there to be at least as may reads as there are threads.
                 * This is necessary due to the comp. graph, where a volatile module cannot be
                 * dry on initialization.
                 */
                if( pQuery->eof( ) )
                    throw std::runtime_error( std::string( "No NucSeqs in database for iSequencerId=" )
                                                  .append( std::to_string( iSequencerId ) )
                                                  .append( " uiModulo=" )
                                                  .append( std::to_string( uiModulo ) )
                                                  .append( " uiRes=" )
                                                  .append( std::to_string( uiRes ) )
                                                  .append( " bAll=" )
                                                  .append( bAll ? "true" : "false" ) );

                return pQuery;
            },
            iConnectionId, iSequencerId, uiRes, uiModulo, bAll );
    } // method
}; // class


/**
 * @brief fetches reads from a database
 */
template <typename DBCon> class NucSeqFetcher : public Module<NucSeq, true, NucSeqQueryContainer<DBCon>>
{
  public:
    NucSeqFetcher( const ParameterSetManager& rParameters )
    {}

    /// @brief returns one read at a time until isFinished returns true.
    virtual std::shared_ptr<NucSeq> EXPORTED execute( std::shared_ptr<NucSeqQueryContainer<DBCon>> pQuery )
    {
        if( pQuery->eof( ) )
            throw std::runtime_error( "No more NucSeqs" );

        auto xTup = pQuery->get( );
        auto pRet = std::get<0>( xTup )->pUncomNucSeq;
        pRet->iId = (int64_t)std::get<1>( xTup );
        if( !pQuery->next( ) )
            return nullptr;
        assert( pRet->iId != -1 );
        return pRet;
    } // method
}; // class

} // namespace libMA


#ifdef WITH_PYTHON
/// @brief expose the NucSeqFromSql classes to python
void exportNucSeqSql( py::module& rxPyModuleId );
#endif
