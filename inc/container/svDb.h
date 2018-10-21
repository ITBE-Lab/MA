/**
 * @file svDb.h
 * @details
 * The database interface for the structural variant caller
 */

#include "container/container.h"
#include "container/nucSeq.h"
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
                                  std::vector<std::string>{"name UNIQUE", "mode"} // column definitions of the table
                                  ),
              pDatabase( pDatabase ),
              xGetSequencerId( *pDatabase, "SELECT id FROM sequencer_table WHERE name == ?" )
        {
            pDatabase->execDML( "CREATE INDEX sequencer_id_index ON sequencer_table (id)" );
        } // default constructor

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
              xGetReadId( *pDatabase, "SELECT id FROM read_table WHERE sequencer_id == ? AND name == ?" ),
              xUpdateReadId( *pDatabase, "UPDATE read_table SET paired_read_id = ? WHERE is == ?" )
        {
            pDatabase->execDML( "CREATE INDEX read_id_index ON read_table (id)" );
        } // default constructor

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
                                                     uint32_t // soc score
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
                            std::vector<std::string>{"read_id", "soc_start", "soc_end", "soc_score"} ),
              pDatabase( pDatabase )
        {
            pDatabase->execDML( "CREATE INDEX soc_start_index ON soc_table (soc_start)" );
        } // default constructor
    }; // class

    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    std::shared_ptr<SequencerTable> pSequencerTable;
    std::shared_ptr<ReadTable> pReadTable;
    std::shared_ptr<SoCTable> pSocTable;

  public:
    SV_DB( std::string sName )
        : pDatabase( new CppSQLiteDBExtended( "./", sName, eCREATE_DB ) ),
          pSequencerTable( new SequencerTable( pDatabase ) ),
          pReadTable( new ReadTable( pDatabase ) ),
          pSocTable( new SoCTable( pDatabase ) )
    {} // constructor

    SV_DB( const SV_DB& rOther ) = delete; // delete copy constructor

    class SoCInserter
    {
      private:
        std::shared_ptr<ReadTable> pReadTable;
        std::shared_ptr<SoCTable> pSocTable;
        size_t uiSequencerId;
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

            inline void operator( )( uint32_t uiStart, uint32_t uiEnd, uint32_t uiScore )
            {
                pSocTable->xInsertRow( uiReadId, uiStart, uiEnd, uiScore );
            } // method
        }; // class

      public:
        SoCInserter( SV_DB& rDB, size_t uiSequencerId )
            : pReadTable( rDB.pReadTable ),
              pSocTable( rDB.pSocTable ),
              uiSequencerId( uiSequencerId ),
              xTransactionContext( *rDB.pDatabase )
        {} // constructor

        inline ReadContex getReadContext( std::shared_ptr<NucSeq> pRead )
        {
            return ReadContex( pReadTable->insertRead( uiSequencerId, pRead ), pSocTable );
        } // method
    }; // class

    inline SoCInserter addSequencerType( std::string& sSequencerName )
    {
        return SoCInserter( *this, pSequencerTable->insertSequencer( sSequencerName ) );
    } // method

}; // class

}; // namespace libMA