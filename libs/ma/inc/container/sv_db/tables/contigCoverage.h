/**
 * @file contigCoverage.h
 * @details
 * Database interface for the structural variant caller.
 * One table of the database.
 */
#pragma once

#include "container/pack.h"
#include "util/sqlite3.h"

namespace libMA
{
class SV_DB;


typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int64_t, // sequencer_id
                                                 int64_t, // contig_nr
                                                 int64_t // num_generated_nt
                                                 >
    TP_CONTIG_COV_TABLE;
class ContigCovTable : public TP_CONTIG_COV_TABLE
{
    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    CppSQLiteExtStatement xIncNt, xResetNt;
    CppSQLiteExtQueryStatement<int64_t> xGetNumNt;
    bool bPrintedCovList = false;

  public:
    ContigCovTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
        : TP_CONTIG_COV_TABLE( *pDatabase, // the database where the table resides
                               "contig_cov_table", // name of the table in the database
                               // column definitions of the table
                               std::vector<std::string>{"sequencer_id", "contig_nr", "num_generated_nt"},
                               // constraints for table
                               std::vector<std::string>{"UNIQUE (sequencer_id, contig_nr)",
                                                        "FOREIGN KEY (sequencer_id) REFERENCES sequencer_table(id)"} ),
          pDatabase( pDatabase ),
          xIncNt( *pDatabase,
                  "UPDATE contig_cov_table "
                  "SET num_generated_nt = num_generated_nt + ? "
                  "WHERE sequencer_id == ? "
                  "AND contig_nr == ? " ),
          xResetNt( *pDatabase,
                    "UPDATE contig_cov_table "
                    "SET num_generated_nt = 0 "
                    "WHERE sequencer_id == ? " ),
          xGetNumNt( *pDatabase,
                     "SELECT num_generated_nt "
                     "FROM contig_cov_table "
                     "WHERE sequencer_id == ? "
                     "ORDER BY contig_nr " )
    {
        // pDatabase->execDML( "CREATE INDEX IF NOT EXISTS sequencer_id_index ON sequencer_table (id)" );
    } // default constructor

    inline int64_t insert( int64_t iSequencerId, int64_t iContigId )
    {
        return xInsertRow( iSequencerId, iContigId, 0 );
    } // method

    inline void insert( int64_t iSequencerId, std::shared_ptr<Pack> pPack )
    {
        for( size_t uiI = 0; uiI < pPack->uiNumContigs( ); uiI++ )
            insert( iSequencerId, uiI );
    } // method

    inline void incrementNt( int64_t iSequencerId, int64_t iContigId, int64_t iAmount )
    {
        xIncNt.bindAndExecute( iAmount, iSequencerId, iContigId );
    } // method

    inline void resetCount( int64_t iSequencerId )
    {
        xResetNt.bindAndExecute( iSequencerId );
    } // method

    inline std::vector<int64_t> getNumNt( int64_t iSequencerId )
    {
        return xGetNumNt.executeAndStoreInVector<0>( iSequencerId );
    } // method

    inline std::vector<double> getEstimatedCoverageList( int64_t iSequencerId, std::shared_ptr<Pack> pPack )
    {
        std::vector<double> vRet;
        auto vNumNt = this->getNumNt( iSequencerId );
        assert( vNumNt.size( ) == pPack->uiNumContigs( ) );
        vRet.reserve( vNumNt.size( ) );
        if( !bPrintedCovList )
        {
            int64_t iTotal = 0;
            for( size_t uiI = 0; uiI < vNumNt.size( ); uiI++ )
                iTotal += vNumNt[ uiI ];
            std::cout << "estimated coverage per contig (showing >= 3x):" << std::endl;
            std::cout << "contig_id\tcoverage\tnum_nt\t%" << std::endl;
            for( size_t uiI = 0; uiI < vNumNt.size( ); uiI++ )
                if( vNumNt[ uiI ] / (double)pPack->xVectorOfSequenceDescriptors[ uiI ].uiLengthUnpacked >= 3 &&
                    100 * vNumNt[ uiI ] >= iTotal )
                    std::cout << uiI << "\t"
                              << ( (int)10 * vNumNt[ uiI ] /
                                   (double)pPack->xVectorOfSequenceDescriptors[ uiI ].uiLengthUnpacked ) /
                                     10.0
                              << "x\t" << vNumNt[ uiI ] << "\t" << ( 100 * vNumNt[ uiI ] ) / iTotal << "%" << std::endl;
            std::cout << std::endl;
            bPrintedCovList = true;
        } // if
        for( size_t uiI = 0; uiI < vNumNt.size( ); uiI++ )
            vRet.push_back( vNumNt[ uiI ] / ( (double)pPack->xVectorOfSequenceDescriptors[ uiI ].uiLengthUnpacked ) );
        return vRet;
    } // method

    class CovInserter
    {
      public:
        int64_t iSequencerId;
        std::shared_ptr<SV_DB> pDb;
        std::shared_ptr<Pack> pPack;
        std::vector<int64_t> vNumNts;
        std::mutex xLock;

        CovInserter( int64_t iSequencerId, std::shared_ptr<Pack> pPack, std::shared_ptr<SV_DB> pDb );

        void commit( );

        ~CovInserter( )
        {
            commit( );
        } // deconstructor

        /// @brief seeds need to be sorted by query pos
        inline void insert( Seeds& rSeeds, nucSeqIndex uiQlen )
        {
            if( iSequencerId == -1 )
                return;
            for( size_t uiI = 0; uiI < rSeeds.size( ); uiI++ )
            {
                // size of this seed
                int64_t iSize = rSeeds[ uiI ].size( );

#if 0 // 1 -> add gap in between seeds to estimation / 0 -> don't
      // add gap to previous seed or start of query
                    if( uiI == 0 )
                        iSize += rSeeds[ uiI ].start( );
                    else if( rSeeds[ uiI ].start( ) < rSeeds[ uiI - 1 ].end( ) )
                        iSize += ( rSeeds[ uiI - 1 ].end( ) - rSeeds[ uiI ].start( ) ) / 2;

                    // add gap to next seed or end of query
                    if( uiI + 1 == rSeeds.size( ) )
                        iSize += uiQlen - rSeeds[ uiI ].end( );
                    else if( rSeeds[ uiI ].end( ) < rSeeds[ uiI + 1 ].start( ) )
                        iSize += ( rSeeds[ uiI + 1 ].start( ) - rSeeds[ uiI ].end( ) ) / 2;
#endif

                // increase the count
                size_t uiIdx = pPack->uiSequenceIdForPosition( rSeeds[ uiI ].start_ref( ) );
                {
                    std::lock_guard<std::mutex> xGuard( xLock );
                    vNumNts[ uiIdx ] += iSize;
                } // scope for xGuard
            } // for
        } // method
    }; // class
}; // class

} // namespace libMA