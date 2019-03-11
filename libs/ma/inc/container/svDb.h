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

        inline uint32_t insertSequencer( std::string& sSequencerName )
        {
            xInsertRow( sSequencerName, "not supported yet" );
            return xGetSequencerId.scalar( sSequencerName );
        } // method
    }; // class


    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int32_t, // sequencer id
                                                     std::string, // read name
                                                     NucSeqSql // read sequence
                                                     >
        TP_READ_TABLE;
    class ReadTable : public TP_READ_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        CppSQLiteExtQueryStatement<int32_t> xGetReadId;
        bool bDoDuplicateWarning = true;

      public:
        // @todo (read_table.name, read_table.sequencer_id) should be UNIQUE
        ReadTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_READ_TABLE(
                  *pDatabase, // the database where the table resides
                  "read_table", // name of the table in the database
                  // column definitions of the table
                  std::vector<std::string>{"sequencer_id", "name", "sequence"},
                  std::vector<std::string>{"read_name_sequencer_id_constraint UNIQUE (sequencer_id, name)"} ),
              pDatabase( pDatabase ),
              xGetReadId( *pDatabase, "SELECT id FROM read_table WHERE sequencer_id = ? AND name = ?" )
        {} // default constructor

        ~ReadTable( )
        {
            if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
                pDatabase->execDML( "CREATE INDEX read_id_index ON read_table (id)" );
        } // deconstructor

        inline int32_t insertRead( int32_t uiSequencerId, std::shared_ptr<NucSeq> pRead )
        {
            for( size_t i = 0; i < 5; i++ ) // at most do 5 tries
            {
                try
                {
                    xInsertRow( uiSequencerId, pRead->sName, NucSeqSql( pRead ) );
                    return xGetReadId.scalar( uiSequencerId, pRead->sName );
                } // try
                catch( CppSQLite3Exception& xException )
                {
                    if( bDoDuplicateWarning )
                    {
                        std::cerr << "WARNING: " << xException.errorMessage( ) << std::endl;
                        std::cerr << "Does your data contain duplicate reads? Current read name: " << pRead->sName
                                  << std::endl;
                        std::cerr << "Changing read name to: " << pRead->sName + "_2" << std::endl;
                        std::cerr << "This warning is only displayed once" << std::endl;
                        bDoDuplicateWarning = false;
                    } // if
                    pRead->sName += "_2";
                } // catch
            } // for
            throw AnnotatedException( "Could not insert read after 5 tries" );
        } // method

    }; // class


    typedef CppSQLiteExtTable<int32_t, // read id
                              int32_t, // sequencer id
                              std::string, // read name
                              NucSeqSql // read sequence
                              >
        TP_PAIRED_READ_TABLE;
    class PairedReadTable : public TP_PAIRED_READ_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        bool bDoDuplicateWarning = true;

      public:
        // @todo (read_table.name, read_table.sequencer_id) should be UNIQUE
        PairedReadTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_PAIRED_READ_TABLE(
                  *pDatabase, // the database where the table resides
                  "paired_read_table", // name of the table in the database
                  // column definitions of the table
                  std::vector<std::string>{"id", "sequencer_id", "name", "sequence"},
                  false, // do not generate automatic primary key...
                  std::vector<std::string>{"paired_read_name_sequencer_id_constraint UNIQUE (sequencer_id, name)",
                                           "paired_read_table_key PRIMARY KEY (id)"} ),
              pDatabase( pDatabase )
        {} // default constructor

        ~PairedReadTable( )
        {
            if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
                pDatabase->execDML( "CREATE INDEX paired_read_id_index ON paired_read_table (id)" );
        } // deconstructor

        inline void insertRead( int32_t uiId, int32_t uiSequencerId, std::shared_ptr<NucSeq> pRead )
        {
            for( size_t i = 0; i < 5; i++ ) // at most do 5 tries
            {
                try
                {
                    xInsertRow( uiId, uiSequencerId, pRead->sName, NucSeqSql( pRead ) );
                    return; // if we reach here we have succesfully inserted the read; no need to repeat the loop
                } // try
                catch( CppSQLite3Exception& xException )
                {
                    if( bDoDuplicateWarning )
                    {
                        this->vDump(std::cout);
                        std::cerr << "WARNING: " << xException.errorMessage( ) << std::endl;
                        std::cerr << "Does your data contain duplicate reads? Current read name: " << pRead->sName
                                  << std::endl;
                        std::cerr << "Changing read name to: " << pRead->sName + "_2" << std::endl;
                        std::cerr << "This warning is only displayed once" << std::endl;
                        bDoDuplicateWarning = false;
                    } // if
                    pRead->sName += "_2";
                } // catch
            } // for
            throw AnnotatedException( "Could not insert read after 5 tries" );
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
                pDatabase->execDML( "CREATE INDEX soc_index ON soc_table (soc_start, soc_end)" );
            } // if
        } // deconstructor
    }; // class


    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<uint32_t, // pos on reference
                                                     uint32_t, // size
                                                     int32_t, // associated SoC id
                                                     std::string // state?
                                                     >
        TP_SV_TABLE;
    class DeletionTable : public TP_SV_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;

      public:
        DeletionTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SV_TABLE( *pDatabase, // the database where the table resides
                           "deletion_table", // name of the table in the database
                           // column definitions of the table
                           std::vector<std::string>{"start_on_ref", "size", "soc_id", "state"} ),
              pDatabase( pDatabase )
        {} // default constructor

        ~DeletionTable( )
        {
            // if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
            // {
            //     pDatabase->execDML( "CREATE INDEX soc_start_index ON soc_table (soc_start)" );
            //     pDatabase->execDML( "CREATE INDEX soc_end_index ON soc_table (soc_end)" );
            // } // if
        } // deconstructor
    }; // class

    class InsertionTable : public TP_SV_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;

      public:
        InsertionTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SV_TABLE( *pDatabase, // the database where the table resides
                           "insertion_table", // name of the table in the database
                           // column definitions of the table
                           std::vector<std::string>{"pos_on_ref", "size", "soc_id", "state"} ),
              pDatabase( pDatabase )
        {} // default constructor

        ~InsertionTable( )
        {
            // if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
            // {
            //     pDatabase->execDML( "CREATE INDEX soc_start_index ON soc_table (soc_start)" );
            //     pDatabase->execDML( "CREATE INDEX soc_end_index ON soc_table (soc_end)" );
            // } // if
        } // deconstructor
    }; // class

    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<uint32_t, // pos on reference
                                                     int32_t, // associated SoC id
                                                     std::string // state?
                                                     >
        TP_MUTATION_TABLE;
    class MutationTable : public TP_MUTATION_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;

      public:
        MutationTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_MUTATION_TABLE( *pDatabase, // the database where the table resides
                                 "mutation_table", // name of the table in the database
                                 // column definitions of the table
                                 std::vector<std::string>{"pos_on_ref", "soc_id", "state"} ),
              pDatabase( pDatabase )
        {} // default constructor

        ~MutationTable( )
        {
            // if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
            // {
            //     pDatabase->execDML( "CREATE INDEX soc_start_index ON soc_table (soc_start)" );
            //     pDatabase->execDML( "CREATE INDEX soc_end_index ON soc_table (soc_end)" );
            // } // if
        } // deconstructor
    }; // class

    // // probably i wont use this
    // typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int32_t, // associated SoC id
    //                                                  std::string // bucket type
    //                                                  >
    //     TP_SV_BUCKET_TABLE;
    // class SvBucketTable : public TP_SV_BUCKET_TABLE
    // {
    //     std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    //
    //   public:
    //     SvBucketTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
    //         : TP_SV_BUCKET_TABLE( *pDatabase, // the database where the table resides
    //                        "sv_bucket_table", // name of the table in the database
    //                        // column definitions of the table
    //                        std::vector<std::string>{"soc_id", "bucket_type"} ),
    //           pDatabase( pDatabase )
    //     {} // default constructor
    //
    //     ~SvBucketTable( )
    //     {
    //         // if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
    //         // {
    //         //     pDatabase->execDML( "CREATE INDEX soc_start_index ON soc_table (soc_start)" );
    //         //     pDatabase->execDML( "CREATE INDEX soc_end_index ON soc_table (soc_end)" );
    //         // } // if
    //     } // deconstructor
    // }; // class

    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    std::shared_ptr<SequencerTable> pSequencerTable;
    std::shared_ptr<ReadTable> pReadTable;
    std::shared_ptr<PairedReadTable> pPairedReadTable;
    std::shared_ptr<SoCTable> pSocTable;

    friend class NucSeqFromSql;
    friend class PairedNucSeqFromSql;

  public:
    SV_DB( std::string sName, enumSQLite3DBOpenMode xMode )
        : pDatabase( new CppSQLiteDBExtended( "", sName, xMode ) ),
          pSequencerTable( new SequencerTable( pDatabase ) ),
          pReadTable( new ReadTable( pDatabase ) ),
          pPairedReadTable( new PairedReadTable( pDatabase ) ),
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

        uint32_t uiSequencerId;
        // must be first so that it is deconstructed first
        CppSQLiteExtImmediateTransactionContext xTransactionContext;

        class ReadContex
        {
          private:
            uint32_t uiReadId;
            std::shared_ptr<SoCTable> pSocTable;

          public:
            ReadContex( uint32_t uiReadId, std::shared_ptr<SoCTable> pSocTable )
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
            throw std::runtime_error( "unimplemented..." );
            // auto xIdPair = pDB->pReadTable->insertPairedRead( uiSequencerId, pReadA, pReadB );
            // return std::make_pair( ReadContex( xIdPair.first, pDB->pSocTable ),
            //                       ReadContex( xIdPair.second, pDB->pSocTable ) );
        } // method

        inline void insertRead( std::shared_ptr<NucSeq> pRead )
        {
            pDB->pReadTable->insertRead( uiSequencerId, pRead );
        } // method

        inline void insertPairedRead( std::shared_ptr<NucSeq> pReadA, std::shared_ptr<NucSeq> pReadB )
        {
            pDB->pPairedReadTable->insertRead( pDB->pReadTable->insertRead( uiSequencerId, pReadA ), uiSequencerId,
                                               pReadB );
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
    SoCDbWriter( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB::SoCInserter> pInserter )
        : pInserter( pInserter ), pMutex( new std::mutex )
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
            xInsertContext( (uint32_t)std::get<0>( xSoCTup ), (uint32_t)std::get<1>( xSoCTup ), std::get<2>( xSoCTup ),
                            std::get<3>( xSoCTup ) );
        } // for

        return std::make_shared<Container>( );
    } // method
}; // class

class NucSeqFromSql : public Module<NucSeq, true>
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t> xQuery;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t>::Iterator xTableIterator;

  public:
    NucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, std::string sSql )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  ( "SELECT read_table.sequence, read_table.id FROM read_table WHERE read_table.id NOT IN (SELECT "
                    "paired_read_table.id FROM paired_read_table) " +
                    sSql )
                      .c_str( ) ),
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
        std::get<0>( xTup ).pNucSeq->sName = std::to_string( std::get<1>( xTup ) );
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
    CppSQLiteExtQueryStatement<NucSeqSql, NucSeqSql, uint32_t> xQuery;
    CppSQLiteExtQueryStatement<NucSeqSql, NucSeqSql, uint32_t>::Iterator xTableIterator;

  public:
    PairedNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, std::string sSql )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  ( "SELECT read_table.sequence, paired_read_table.sequence, read_table.id FROM read_table INNER JOIN "
                    "paired_read_table ON "
                    "read_table.id = paired_read_table.id " +
                    sSql )
                      .c_str( ) ),
          xTableIterator( xQuery.vExecuteAndReturnIterator( ) )
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
        pRet->back( )->sName = std::to_string( std::get<2>( xTup ) );
        pRet->push_back( std::get<1>( xTup ).pNucSeq );
        pRet->back( )->sName = std::to_string( std::get<2>( xTup ) );
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

}; // namespace libMA

#ifdef WITH_PYTHON
#ifdef WITH_BOOST
void exportSoCDbWriter( );
#else
void exportSoCDbWriter( py::module& rxPyModuleId );
#endif
#endif
