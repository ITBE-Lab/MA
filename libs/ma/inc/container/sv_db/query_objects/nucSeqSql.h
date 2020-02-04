/**
 * @file nucSeqSql.h
 * @brief implements libMA::AllNucSeqFromSql that fetches reads from the DB.
 * @author Markus Schmidt
 */
#pragma once
#include "container/sv_db/svSchema.h"
#include "container/sv_db/connection_container.h"

namespace libMA
{

template <typename DBCon> class NucSeqQueryContainer : public SQLQuery<DBCon, NucSeqSql, int64_t>, public Container
{
    using SQLQuery::SQLQuery;
}; // class

template <typename DBCon> class NucSeqFromSqlQuery : public Module<NucSeqQueryContainer<DBCon>, false, DBCon>
{
    int64_t iSequencerId;
    uint32_t uiRes;
    uint32_t uiModulo;
    bool bAll;
    bool bUnpaired;

    NucSeqFromSqlQuery( const ParameterSetManager& rParameters, int64_t iSequencerId, size_t uiRes, size_t uiModulo,
                        bool bPaired, bool bUnpaired )
        : iSequencerId( iSequencerId ),
          uiRes( (uint32_t)uiRes ),
          uiModulo( (uint32_t)uiModulo ),
          bAll( bPaired && bUnpaired ),
          bUnpaired( bUnpaired )
    {
        if( !bPaired && !bUnpaired )
            throw std::runtime_error( "NucSeqFromSqlQuery is not fetching anything" );
        if( !bAll && !bUnpaired )
            throw std::runtime_error( "NucSeqFromSqlQuery: fetching unpaired NucSeqs only is unimplemented" );
    } // constructor

    /// @brief returns a query that can fetch NucSeqs.
    std::shared_ptr<NucSeqQueryContainer<DBCon>> execute( std::shared_ptr<ConnectionContainer<DBCon>> pConnection )
    {
        auto pQuery = std::make_shared<NucSeqQueryContainer<DBCon>>(
            pConnection,
            std::string( "SELECT read_table.sequence, read_table.id "
                         "FROM read_table " )
                .append( iSequencerId != -1 && uiModulo != 1 ? "WHERE sequencer_id = ? "
                                                               "AND read_table.id % ? = ? "
                                                             : "" )
                .append( bAll ? "" : ( iSequencerId != -1 && uiModulo != 1 ? "AND" : "WHERE" ) )
                .append( bAll ? ""
                              : " read_table.id NOT IN ( "
                                "   SELECT paired_read_table.first_read FROM paired_read_table "
                                "   UNION "
                                "   SELECT paired_read_table.second_read FROM paired_read_table "
                                ") " ) );

        if( iSequencerId != -1 && uiModulo != 1 )
            pQuery->execAndFetch( iSequencerId, uiRes, uiModulo );
        else if( iSequencerId != -1 )
            pQuery->execAndFetch( iSequencerId );
        else
            pQuery->execAndFetch( );

        if( pQuery->eof( ) )
            throw AnnotatedException( "No NucSeqs in database" );

        return pQuery;
    } // method
}; // class


/**
 * @brief fetches reads from a database
 */
template <typename DBCon> class NucSeqFetcher : public Module<NucSeq, true, NucSeqQueryContainer<DBCon>>
{
  public:
    NucSeqFetcher( )
    {}

    /// @brief returns one read at a time until isFinished returns true.
    std::shared_ptr<NucSeq> execute( std::shared_ptr<NucSeqQueryContainer<DBCon>> pQuery )
    {
        if( pQuery->eof( ) )
            throw AnnotatedException( "No more NucSeqs" );

        auto xTup = pQuery->get( );
        std::get<0>( xTup ).pNucSeq->iId = std::get<1>( xTup );
        pQuery->next( );
        if( pQuery->eof( ) )
            setFinished( );
        return std::get<0>( xTup ).pNucSeq;
    } // method
}; // class

} // namespace libMA


#ifdef WITH_PYTHON
/// @brief expose the NucSeqFromSql classes to python
void exportNucSeqSql( py::module& rxPyModuleId );
#endif
