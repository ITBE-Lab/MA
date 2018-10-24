/**
 * @file svDb.h
 * @details
 * The database interface for the structural variant caller
 */

#include "container/container.h"
#include "container/nucSeq.h"
#include "container/soc.h"
#include "module/module.h"
#include "util/exception.h"
#include "util/sqlite3.h"

namespace libMA
{

class SV_DB : public CppSQLite3DB, public Container
{
  private:
    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string, // sequencer name
                                                     std::string // parameter setting (for the moment merely the mode)
                                                     >
        TP_SEQUENCER_TABLE;
    class SequencerTable : public TP_SEQUENCER_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        CppSQLiteExtQueryStatement<int32_t> xGetSequencerId;

      public:
        // @todo sequencer_table.name should be UNIQUE (this requires that CppSQLiteExtTable can deal with constrints)
        SequencerTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SEQUENCER_TABLE( *pDatabase, // the database where the table resides
                                  "sequencer_table", // name of the table in the database
                                  // column definitions of the table
                                  std::vector<std::string>{"name", "mode"},
                                  std::vector<std::string>{"unique_sequencer_name_constraint UNIQUE (name)"} ),
              pDatabase( pDatabase ),
              xGetSequencerId( *pDatabase, "SELECT id FROM sequencer_table WHERE name = ?" )
        {} // default constructor

        ~SequencerTable( )
        {
            if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
                pDatabase->execDML( "CREATE INDEX sequencer_id_index ON sequencer_table (id)" );
        } // deconstructor

        inline size_t insertSequencer( std::string& sSequencerName )
        {
            xInsertRow( sSequencerName, "not supported yet" );
            return xGetSequencerId.scalar( sSequencerName );
        } // method
    }; // class


    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int32_t, // sequencer id
                                                     std::string, // read name
                                                     NucSeqSql, // read sequence
                                                     int32_t // paired read id
                                                     >
        TP_READ_TABLE;
    class ReadTable : public TP_READ_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        CppSQLiteExtQueryStatement<int32_t> xGetReadId;
        CppSQLiteExtStatement xUpdateReadId;

      public:
        // @todo (read_table.name, read_table.sequencer_id) should be UNIQUE
        ReadTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_READ_TABLE(
                  *pDatabase, // the database where the table resides
                  "read_table", // name of the table in the database
                  // column definitions of the table
                  std::vector<std::string>{"sequencer_id", "name", "sequence", "paired_read_id"},
                  std::vector<std::string>{"read_name_sequencer_id_constraint UNIQUE (sequencer_id, name)"} ),
              pDatabase( pDatabase ),
              xGetReadId( *pDatabase, "SELECT id FROM read_table WHERE sequencer_id = ? AND name = ?" ),
              xUpdateReadId( *pDatabase, "UPDATE read_table SET paired_read_id = ? WHERE id = ?" )
        {} // default constructor

        ~ReadTable( )
        {
            if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
                pDatabase->execDML( "CREATE INDEX read_id_index ON read_table (id)" );
        } // deconstructor

        inline int32_t insertRead( int32_t uiSequencerId, std::shared_ptr<NucSeq> pRead )
        {
            xInsertRow( uiSequencerId, pRead->sName, NucSeqSql( pRead ), -1 );
            return xGetReadId.scalar( uiSequencerId, pRead->sName );
        } // method

        inline std::pair<int32_t, int32_t> insertPairedRead( int32_t uiSequencerId, //
                                                             std::shared_ptr<NucSeq>
                                                                 pReadA,
                                                             std::shared_ptr<NucSeq>
                                                                 pReadB )
        {
            // insert both reads
            xInsertRow( uiSequencerId, pReadA->sName, NucSeqSql( pReadA ), -1 );
            size_t uiReadAid = xGetReadId.scalar( uiSequencerId, pReadA->sName );
            xInsertRow( uiSequencerId, pReadB->sName, NucSeqSql( pReadB ), -1 );
            size_t uiReadBid = xGetReadId.scalar( uiSequencerId, pReadB->sName );
            // update the paired read ids
            xUpdateReadId.bindAndExecute( uiReadAid, uiReadBid );
            xUpdateReadId.bindAndExecute( uiReadBid, uiReadAid );
            // return the id's
            return std::make_pair( uiReadAid, uiReadBid );
        } // method
    }; // class


    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int32_t, // read id
                                                     uint32_t, // soc start
                                                     uint32_t, // soc end
                                                     uint32_t, // soc score
                                                     double // strand ratio
                                                     >
        TP_SOC_TABLE;
    class SoCTable : public TP_SOC_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;

      public:
        SoCTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SOC_TABLE( *pDatabase, // the database where the table resides
                            "soc_table", // name of the table in the database
                            // column definitions of the table
                            std::vector<std::string>{"read_id", "soc_start", "soc_end", "soc_score", "strand_ratio"} ),
              pDatabase( pDatabase )
        {} // default constructor

        ~SoCTable( )
        {
            if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
            {
                pDatabase->execDML( "CREATE INDEX soc_start_index ON soc_table (soc_start)" );
                pDatabase->execDML( "CREATE INDEX soc_end_index ON soc_table (soc_end)" );
            } // if
        } // deconstructor
    }; // class

    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    std::shared_ptr<SequencerTable> pSequencerTable;
    std::shared_ptr<ReadTable> pReadTable;
    std::shared_ptr<SoCTable> pSocTable;

    friend class NucSeqFromSql;

  public:
    SV_DB( std::string sName, enumSQLite3DBOpenMode xMode )
        : pDatabase( new CppSQLiteDBExtended( "./", sName, xMode ) ),
          pSequencerTable( new SequencerTable( pDatabase ) ),
          pReadTable( new ReadTable( pDatabase ) ),
          pSocTable( new SoCTable( pDatabase ) )
    {} // constructor

    SV_DB( std::string sName ) : SV_DB( sName, eCREATE_DB )
    {} // constructor

    SV_DB( std::string sName, std::string sMode ) : SV_DB( sName, sMode == "create" ? eCREATE_DB : eOPEN_DB )
    {} // constructor

    SV_DB( const SV_DB& rOther ) = delete; // delete copy constructor

    class SoCInserter
    {
      private:
        // this is here so that it gets destructed after the transaction context
        std::shared_ptr<SV_DB> pDB;

        size_t uiSequencerId;
        // must be first so that it is deconstructed first
        CppSQLiteExtImmediateTransactionContext xTransactionContext;

        class ReadContex
        {
          private:
            size_t uiReadId;
            std::shared_ptr<SoCTable> pSocTable;

          public:
            ReadContex( size_t uiReadId, std::shared_ptr<SoCTable> pSocTable )
                : uiReadId( uiReadId ), pSocTable( pSocTable )
            {} // constructor

            inline void operator( )( uint32_t uiStart, uint32_t uiEnd, uint32_t uiScore, double dStrandRatio )
            {
                pSocTable->xInsertRow( uiReadId, uiStart, uiEnd, uiScore, dStrandRatio );
            } // method
        }; // class

      public:
        SoCInserter( std::shared_ptr<SV_DB> pDB, std::string sSequencerName )
            : pDB( pDB ),
              uiSequencerId( pDB->pSequencerTable->insertSequencer( sSequencerName ) ),
              xTransactionContext( *pDB->pDatabase )
        {} // constructor

        SoCInserter( const SoCInserter& rOther ) = delete;

        inline ReadContex getReadContext( std::shared_ptr<NucSeq> pRead )
        {
            return ReadContex( pDB->pReadTable->insertRead( uiSequencerId, pRead ), pDB->pSocTable );
        } // method

        inline std::pair<ReadContex, ReadContex> getPairedReadContext( //
            std::shared_ptr<NucSeq> pReadA, std::shared_ptr<NucSeq> pReadB )
        {
            auto xIdPair = pDB->pReadTable->insertPairedRead( uiSequencerId, pReadA, pReadB );
            return std::make_pair( ReadContex( xIdPair.first, pDB->pSocTable ),
                                   ReadContex( xIdPair.second, pDB->pSocTable ) );
        } // method
    }; // class

}; // class

class SoCDbWriter : public Module<Container, false, NucSeq, SoCPriorityQueue>
{
  private:
    std::shared_ptr<SV_DB::SoCInserter> pInserter;
    std::shared_ptr<std::mutex> pMutex;
    size_t uiNumSoCsToRecord = 1; // @todo expose

  public:
    SoCDbWriter( std::shared_ptr<SV_DB::SoCInserter> pInserter ) : pInserter( pInserter ), pMutex( new std::mutex )
    {} // constructor

    std::shared_ptr<Container> execute( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<SoCPriorityQueue> pSoCQueue )
    {
        std::lock_guard<std::mutex> xGuard( *pMutex );
        auto xInsertContext = pInserter->getReadContext( pQuery );
        for( size_t uiI = 0; uiI < uiNumSoCsToRecord; uiI++ )
        {
            if( pSoCQueue->empty( ) )
                break;

            auto xSoCTup = pSoCQueue->pop_info( );
            xInsertContext( std::get<0>( xSoCTup ), std::get<1>( xSoCTup ), std::get<2>( xSoCTup ),
                            std::get<3>( xSoCTup ) );
        } // for

        return std::make_shared<Container>( );
    } // method
}; // class

class NucSeqFromSql : public Module<NucSeq, true>
{
    std::vector<NucSeqSql> vSequences;

  public:
    NucSeqFromSql( std::shared_ptr<SV_DB> pDb, std::string sSql )
    {
        CppSQLiteExtQueryStatement<NucSeqSql> xSql( *pDb->pDatabase,
                                                    ( "SELECT sequence FROM read_table " + sSql ).c_str( ) );
        vSequences = xSql.executeAndStoreInVector<0>( );
        if( vSequences.empty( ) )
            setFinished( );
    } // constructor

    std::shared_ptr<NucSeq> execute( )
    {
        if( vSequences.empty( ) )
            throw AnnotatedException( "No more NucSeq in NucSeqFromSql module" );
        auto pRet = vSequences.back( ).pNucSeq;
        vSequences.pop_back( );
        if( vSequences.empty( ) )
            setFinished( );
        return pRet;
    } // method

    // override
    bool requiresLock( ) const
    {
        return true;
    } // method
}; // class

}; // namespace libMA

#ifdef WITH_PYTHON
void exportSoCDbWriter( );
#endif
