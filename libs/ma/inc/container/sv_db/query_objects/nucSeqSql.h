#include "container/sv_db/svDb.h"

#pragma once

namespace libMA
{

// @todo this does not make the process sufficiently parallel
class AllNucSeqFromSql : public Module<NucSeq, true>
{
    /// this wrapper is required so that the iterator is never copied
    class InteratorHolder
    {
      public:
        CppSQLiteExtQueryStatement<NucSeqSql, int64_t>::Iterator xIterator;

        InteratorHolder( CppSQLiteExtQueryStatement<NucSeqSql, int64_t>& xQuery, int64_t iSequencerId, uint32_t uiRes,
                         uint32_t uiModulo )
            : xIterator( xQuery.vExecuteAndReturnIterator( iSequencerId, uiModulo, uiRes ) )
        {} // constructor

        InteratorHolder( CppSQLiteExtQueryStatement<NucSeqSql, int64_t>& xQuery, int64_t iSequencerId )
            : xIterator( xQuery.vExecuteAndReturnIterator( iSequencerId ) )
        {} // constructor

        InteratorHolder( CppSQLiteExtQueryStatement<NucSeqSql, int64_t>& xQuery )
            : xIterator( xQuery.vExecuteAndReturnIterator( ) )
        {} // constructor
    }; // class
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, int64_t> xQuery;
    std::shared_ptr<InteratorHolder> pTableIterator;
    int64_t iSequencerId;
    uint32_t uiRes;
    uint32_t uiModulo;

  public:
    AllNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
                  "SELECT read_table.sequence, read_table.id "
                  "FROM read_table " ),
          iSequencerId( -1 ),
          uiRes( 0 ),
          uiModulo( 0 )
    {} // constructor

    AllNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSequencerId,
                      size_t uiRes, size_t uiModulo )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
                  ( uiModulo != 1 ? "SELECT read_table.sequence, read_table.id "
                                    "FROM read_table "
                                    "WHERE sequencer_id == ? "
                                    "AND read_table.id % ? == ? "
                                  : "SELECT read_table.sequence, read_table.id "
                                    "FROM read_table "
                                    "WHERE sequencer_id == ? " ) ),
          iSequencerId( iSequencerId ),
          uiRes( (uint32_t)uiRes ),
          uiModulo( (uint32_t)uiModulo )
    {
#if DEBUG_LEVEL > 0
#if 0
        std::cout << "AllNucSeqFromSql::xQuery" << std::endl;
        xQuery.bindAndExplain( iSequencerId, uiModulo, uiRes );
#endif
#endif
    } // constructor

    std::shared_ptr<NucSeq> execute( )
    {
        if( pTableIterator == nullptr && iSequencerId != -1 && uiModulo != 1 )
            pTableIterator = std::make_unique<InteratorHolder>( xQuery, iSequencerId, uiRes, uiModulo );
        else if( pTableIterator == nullptr && iSequencerId != -1 )
            pTableIterator = std::make_unique<InteratorHolder>( xQuery, iSequencerId );
        else if( pTableIterator == nullptr )
            pTableIterator = std::make_unique<InteratorHolder>( xQuery );

        if( pTableIterator->xIterator.eof( ) )
            throw AnnotatedException( "No more NucSeq in NucSeqFromSql module" );

        auto xTup = pTableIterator->xIterator.get( );
        // std::get<0>( xTup ).pNucSeq->sName = std::to_string( std::get<1>( xTup ) );
        std::get<0>( xTup ).pNucSeq->iId = std::get<1>( xTup );
        pTableIterator->xIterator.next( );

        if( pTableIterator->xIterator.eof( ) )
            setFinished( );
        return std::get<0>( xTup ).pNucSeq;
    } // method
}; // class

class NucSeqFromSql : public Module<NucSeq, true>
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t> xQuery;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t>::Iterator xTableIterator;

  public:
    NucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSequencerId )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
                  "SELECT read_table.sequence, read_table.id "
                  "FROM read_table "
                  "WHERE read_table.id NOT IN ( "
                  "   SELECT paired_read_table.first_read FROM paired_read_table "
                  "   UNION "
                  "   SELECT paired_read_table.second_read FROM paired_read_table "
                  ") "
                  "AND sequencer_id = ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSequencerId ) )
    {
        if( xTableIterator.eof( ) )
            setFinished( );
    } // constructor

    NucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
                  "SELECT read_table.sequence, read_table.id "
                  "FROM read_table "
                  "WHERE read_table.id NOT IN ( "
                  "   SELECT paired_read_table.first_read FROM paired_read_table "
                  "   UNION "
                  "   SELECT paired_read_table.second_read FROM paired_read_table "
                  ") " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( ) )
    {
        if( xTableIterator.eof( ) )
            setFinished( );
    } // constructor

    std::shared_ptr<NucSeq> execute( )
    {
        if( xTableIterator.eof( ) )
            throw AnnotatedException( "No more NucSeq in NucSeqFromSql module" );

        auto xTup = xTableIterator.get( );
        // std::get<0>( xTup ).pNucSeq->sName = std::to_string( std::get<1>( xTup ) );
        std::get<0>( xTup ).pNucSeq->iId = std::get<1>( xTup );
        xTableIterator.next( );

        if( xTableIterator.eof( ) )
            setFinished( );
        return std::get<0>( xTup ).pNucSeq;
    } // method

    // override
    bool requiresLock( ) const
    {
        return true;
    } // method
}; // class

class PairedNucSeqFromSql : public Module<ContainerVector<std::shared_ptr<NucSeq>>, true>
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, NucSeqSql, uint32_t, uint32_t> xQuery;
    CppSQLiteExtQueryStatement<NucSeqSql, NucSeqSql, uint32_t, uint32_t>::Iterator xTableIterator;
    const bool bRevCompMate;

  public:
    PairedNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, int64_t iSequencerId )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
                  "SELECT A.sequence, B.sequence, A.id, B.id "
                  "FROM read_table A, read_table B "
                  "INNER JOIN paired_read_table "
                  "ON paired_read_table.first_read == A.id "
                  "AND paired_read_table.second_read == B.id "
                  "AND A.sequencer_id = ? " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( iSequencerId ) ),
          bRevCompMate( rParameters.getSelected( )->xRevCompPairedReadMates->get( ) )
    {
        if( xTableIterator.eof( ) )
            setFinished( );
    } // constructor

    PairedNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( std::make_shared<SV_DB>( *pDb ) ),
          xQuery( *this->pDb->pDatabase,
                  "SELECT A.sequence, B.sequence, A.id, B.id "
                  "FROM read_table A, read_table B "
                  "INNER JOIN paired_read_table "
                  "ON paired_read_table.first_read == A.id "
                  "AND paired_read_table.second_read == B.id " ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( ) ),
          bRevCompMate( rParameters.getSelected( )->xRevCompPairedReadMates->get( ) )
    {
        if( xTableIterator.eof( ) )
            setFinished( );
    } // constructor

    std::shared_ptr<ContainerVector<std::shared_ptr<NucSeq>>> execute( )
    {
        if( xTableIterator.eof( ) )
            throw AnnotatedException( "No more NucSeq in PairedNucSeqFromSql module" );

        auto pRet = std::make_shared<ContainerVector<std::shared_ptr<NucSeq>>>( );

        auto xTup = xTableIterator.get( );
        pRet->push_back( std::get<0>( xTup ).pNucSeq );
        pRet->back( )->iId = std::get<2>( xTup );
        pRet->push_back( std::get<1>( xTup ).pNucSeq );
        pRet->back( )->iId = std::get<3>( xTup );

        if( bRevCompMate )
        {
            pRet->back( )->vReverse( );
            pRet->back( )->vSwitchAllBasePairsToComplement( );
        } // if

        xTableIterator.next( );

        if( xTableIterator.eof( ) )
            setFinished( );
        return pRet;
    } // method

    // override
    bool requiresLock( ) const
    {
        return true;
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
void exportNucSeqSql( py::module& rxPyModuleId );
#endif