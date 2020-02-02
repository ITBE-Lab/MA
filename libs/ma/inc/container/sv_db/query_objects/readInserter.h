/**
 * @file readInserter.h
 * @brief implements libMA::ReadInserter that inserts reads into the DB.
 * @author Markus Schmidt
 */
// The order of the following two includes is significant with MSVC.

#pragma once

#include "container/sv_db/svSchema.h"

namespace libMA
{
// #define WITH_TRANSACTION_CONTEXT

/// @brief inserts reads into the DB via a transaction
template <typename DBCon> class _ReadInserter
{
  private:
    // this is here so that it gets destructed after the transaction context
    std::shared_ptr<SV_Schema<DBCon>> pDB;
// must be after the DB so that it is deconstructed first
#ifdef WITH_TRANSACTION_CONTEXT
    // One transaction throughout the existence of the object.
    CppSQLiteExtImmediateTransactionContext xTransactionContext;
#endif

  public:
    /// @brief the sequencer id this inserter is attached to
    int64_t uiSequencerId;

    /// @brief creates a new sequencer entry with the name sSequencerName and transaction
    _ReadInserter( std::shared_ptr<SV_Schema<DBCon>> pDB, std::string sSequencerName, std::shared_ptr<Pack> pPack )
        : pDB( pDB ),
#ifdef WITH_TRANSACTION_CONTEXT
          xTransactionContext( *pDB->pDatabase ),
#endif
          uiSequencerId( pDB->pSequencerTable->insertSequencer( sSequencerName ) )
    {} // constructor

    /// @brief cannot be copied
    _ReadInserter( const _ReadInserter& rOther ) = delete; // delete copy constructor

    /// @brief insert a unpaired read
    inline void insertRead( std::shared_ptr<NucSeq> pRead )
    {
        std::lock_guard<std::mutex> xGuard( *pDB->pWriteLock );
        pDB->pReadTable->insertRead( uiSequencerId, pRead );
    } // method

    /// @brief insert paired reads
    inline void insertPairedRead( std::shared_ptr<NucSeq> pReadA, std::shared_ptr<NucSeq> pReadB )
    {
        std::lock_guard<std::mutex> xGuard( *pDB->pWriteLock );
        pDB->pPairedReadTable->insertRead( uiSequencerId, pReadA, pReadB );
    } // method

    /// @brief insert a file of unpaired reads
    inline void insertFastaFiles( const ParameterSetManager& rParameters, const std::vector<fs::path>& vsFileNames )
    {
        ThreadPool xPool( rParameters.getNumThreads( ) );

        for( size_t uiT = 0; uiT < vsFileNames.size( ); uiT++ )
            xPool.enqueue(
                []( size_t, const std::vector<fs::path>* pvsFileNames, _ReadInserter* pInserter, size_t _uiT,
                    const ParameterSetManager* pParameters ) {
                    FileReader xReader( *pParameters, ( *pvsFileNames )[ _uiT ] );
                    while( !xReader.isFinished( ) )
                        pInserter->insertRead( xReader.execute( ) );
                    return 0;
                },
                &vsFileNames, this, uiT, &rParameters );
    } // method

    /// @brief insert files of paired reads
    inline void insertPairedFastaFiles( const ParameterSetManager& rParameters,
                                        const std::vector<fs::path>& vsFileNames1,
                                        const std::vector<fs::path>& vsFileNames2 )
    {
        ThreadPool xPool( rParameters.getNumThreads( ) );

        assert( vsFileNames1.size( ) == vsFileNames2.size( ) );

        for( size_t uiT = 0; uiT < vsFileNames1.size( ); uiT++ )
            xPool.enqueue(
                []( size_t, const std::vector<fs::path>* pvsFileNames1, const std::vector<fs::path>* pvsFileNames2,
                    _ReadInserter* pInserter, size_t _uiT, const ParameterSetManager* pParameters ) {
                    FileReader xReader1( *pParameters, ( *pvsFileNames1 )[ _uiT ] );
                    FileReader xReader2( *pParameters, ( *pvsFileNames2 )[ _uiT ] );
                    while( !xReader1.isFinished( ) && !xReader2.isFinished( ) )
                        pInserter->insertPairedRead( xReader1.execute( ), xReader2.execute( ) );
                    assert( xReader1.isFinished( ) == xReader2.isFinished( ) );
                    return 0;
                },
                &vsFileNames1, &vsFileNames2, this, uiT, &rParameters );
    } // method
}; // class

} // namespace libMA


#ifdef WITH_PYTHON
/// @brief expose libMA::ReadInsert to python
void exportReadInserter( py::module& rxPyModuleId );
#endif
