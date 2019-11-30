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
    {
        pDB->pContigCovTable->insert( uiSequencerId, pPack );
    } // constructor

    /// @brief cannot be copied
    ReadInserter( const ReadInserter& rOther ) = delete; // delete copy constructor

    /// @brief insert a unpaired read
    inline void insertRead( std::shared_ptr<NucSeq> pRead )
    {
        pDB->pReadTable->insertRead( uiSequencerId, pRead );
    } // method

    /// @brief insert paired reads
    inline void insertPairedRead( std::shared_ptr<NucSeq> pReadA, std::shared_ptr<NucSeq> pReadB )
    {
        pDB->pPairedReadTable->insertRead( uiSequencerId, pReadA, pReadB );
    } // method

    /// @brief insert a file of unpaired reads
    inline void insertFastaFiles( const ParameterSetManager& rParameters, const std::vector<fs::path>& vsFileNames )
    {
        FileListReader xReader( rParameters, vsFileNames );
        {
            ThreadPool xPool( 4 );
            std::mutex xReadLock, xWriteLock;

            xPool.enqueue(
                []( size_t uiTid, FileListReader* pReader, ReadInserter* pInserter, std::mutex* pReadLock,
                    std::mutex* pWriteLock ) {
                    while( !pReader->isFinished( ) )
                    {
                        std::shared_ptr<NucSeq> pRead;
                        {
                            std::lock_guard<std::mutex> xGuard( *pReadLock );
                            pRead = pReader->execute( );
                        } // scope for xGuard
                        {
                            std::lock_guard<std::mutex> xGuard( *pWriteLock );
                            pInserter->insertRead( pRead );
                        } // scope for xGuard
                    } // while
                    return 0;
                },
                &xReader, this, &xReadLock, &xWriteLock );
        } // scope for thread pool
    } // method

    /// @brief insert files of paired reads
    inline void insertPairedFastaFiles( const ParameterSetManager& rParameters,
                                        const std::vector<fs::path>& vsFileNames1,
                                        const std::vector<fs::path>& vsFileNames2 )
    {
        PairedFileReader xReader( rParameters, vsFileNames1, vsFileNames2 );
        {
            ThreadPool xPool( 4 );
            std::mutex xReadLock, xWriteLock;

            xPool.enqueue(
                []( size_t uiTid, PairedFileReader* pReader, ReadInserter* pInserter, std::mutex* pReadLock,
                    std::mutex* pWriteLock ) {
                    while( !pReader->isFinished( ) )
                    {
                        std::shared_ptr<TP_PAIRED_READS> pvReads;
                        {
                            std::lock_guard<std::mutex> xGuard( *pReadLock );
                            pvReads = pReader->execute( );
                        } // scope for xGuard
                        {
                            std::lock_guard<std::mutex> xGuard( *pWriteLock );
                            pInserter->insertPairedRead( ( *pvReads )[ 0 ], ( *pvReads )[ 1 ] );
                        } // scope for xGuard
                    } // while
                    return 0;
                },
                &xReader, this, &xReadLock, &xWriteLock );
        } // scope for thread pool
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/// @brief expose libMA::ReadInsert to python
void exportReadInserter( py::module& rxPyModuleId );
#endif