/**
 * @file readInserter.h
 * @brief implements libMA::ReadInserter that inserts reads into the DB.
 * @author Markus Schmidt
 */
#include "container/sv_db/svDb.h"

#pragma once

namespace libMA
{

/// @brief inserts reads into the DB via a transaction
class ReadInserter
{
  private:
    // this is here so that it gets destructed after the transaction context
    std::shared_ptr<SV_DB> pDB;
    // must be after the DB so that it is deconstructed first
    CppSQLiteExtImmediateTransactionContext xTransactionContext;

  public:
    /// @brief the sequencer id this inserter is attached to
    int64_t uiSequencerId;

    /// @brief creates a new sequencer entry with the name sSequencerName and transaction
    ReadInserter( std::shared_ptr<SV_DB> pDB, std::string sSequencerName, std::shared_ptr<Pack> pPack )
        : pDB( pDB ),
          xTransactionContext( *pDB->pDatabase ),
          uiSequencerId( pDB->pSequencerTable->insertSequencer( sSequencerName ) )
    {} // constructor

    /// @brief cannot be copied
    ReadInserter( const ReadInserter& rOther ) = delete; // delete copy constructor

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
                []( size_t, const std::vector<fs::path>* pvsFileNames, ReadInserter* pInserter, size_t uiT,
                    const ParameterSetManager* pParameters ) {
                    FileReader xReader( *pParameters, (*pvsFileNames)[ uiT ] );
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
                    ReadInserter* pInserter, size_t uiT, const ParameterSetManager* pParameters ) {
                    FileReader xReader1( *pParameters, (*pvsFileNames1)[ uiT ] );
                    FileReader xReader2( *pParameters, (*pvsFileNames2)[ uiT ] );
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