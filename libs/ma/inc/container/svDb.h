/**
 * @file svDb.h
 * @details
 * The database interface for the structural variant caller
 */

#include "container/container.h"
#include "container/nucSeq.h"
#include "container/soc.h"
#include "container/svJump.h"
#include "module/module.h"
#include "util/exception.h"
#include "util/sqlite3.h"

namespace libMA
{

class SV_DB : public CppSQLite3DB, public Container
{
  private:
    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string // sequencer name
                                                     >
        TP_SEQUENCER_TABLE;
    class SequencerTable : public TP_SEQUENCER_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;

      public:
        SequencerTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SEQUENCER_TABLE( *pDatabase, // the database where the table resides
                                  "sequencer_table", // name of the table in the database
                                  // column definitions of the table
                                  std::vector<std::string>{"name"},
                                  // constraints for table
                                  std::vector<std::string>{"UNIQUE (name)"} ),
              pDatabase( pDatabase )
        {} // default constructor

        ~SequencerTable( )
        {
            if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
                pDatabase->execDML( "CREATE INDEX sequencer_id_index ON sequencer_table (id)" );
        } // deconstructor

        inline int64_t insertSequencer( std::string& sSequencerName )
        {
            return xInsertRow( sSequencerName );
        } // method
    }; // class

    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int64_t, // sequencer id (foreign key)
                                                     std::string, // read name
                                                     NucSeqSql // read sequence
                                                     >
        TP_READ_TABLE;
    class ReadTable : public TP_READ_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        bool bDoDuplicateWarning = true;

      public:
        CppSQLiteExtQueryStatement<int32_t> xGetReadId;

        ReadTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_READ_TABLE( *pDatabase, // the database where the table resides
                             "read_table", // name of the table in the database
                             // column definitions of the table
                             std::vector<std::string>{"sequencer_id", "name", "sequence"},
                             // constraints for table
                             std::vector<std::string>{"UNIQUE (sequencer_id, name)", //
                                                      "FOREIGN KEY (sequencer_id) REFERENCES sequencer_table(id)"} ),
              pDatabase( pDatabase ),
              xGetReadId( *pDatabase, "SELECT id FROM read_table WHERE sequencer_id = ? AND name = ?" )
        {} // default constructor

        ~ReadTable( )
        {
            if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
                pDatabase->execDML( "CREATE INDEX read_id_index ON read_table (id)" );
        } // deconstructor

        inline int64_t insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pRead )
        {
            for( size_t i = 0; i < 5; i++ ) // at most do 5 tries
            {
                try
                {
                    return xInsertRow( uiSequencerId, pRead->sName, NucSeqSql( pRead ) );
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

    typedef CppSQLiteExtTable<int64_t, // first read (foreign key)
                              int64_t // second read (foreign key)
                              >
        TP_PAIRED_READ_TABLE;
    class PairedReadTable : public TP_PAIRED_READ_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        std::shared_ptr<ReadTable> pReadTable;

      public:
        PairedReadTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase, std::shared_ptr<ReadTable> pReadTable );

        ~PairedReadTable( )
        {} // deconstructor

        inline std::pair<int64_t, int64_t> insertRead( int64_t uiSequencerId, std::shared_ptr<NucSeq> pReadA,
                                                       std::shared_ptr<NucSeq> pReadB )
        {
            int64_t uiReadIdA = pReadTable->insertRead( uiSequencerId, pReadA );
            int64_t uiReadIdB = pReadTable->insertRead( uiSequencerId, pReadB );
            xInsertRow( uiReadIdA, uiReadIdB );
            return std::make_pair( uiReadIdA, uiReadIdB );
        } // method
    }; // class

    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string, // name
                                                     std::string // desc
                                                     >
        TP_SV_CALLER_RUN_TABLE;
    class SvCallerRunTable : public TP_SV_CALLER_RUN_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;

      public:
        SvCallerRunTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SV_CALLER_RUN_TABLE( *pDatabase, // the database where the table resides
                                      "sv_caller_run_table", // name of the table in the database
                                      // column definitions of the table
                                      std::vector<std::string>{"name", "desc"} ),
              pDatabase( pDatabase )
        {} // default constructor

        ~SvCallerRunTable( )
        {
            if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
                pDatabase->execDML( "CREATE INDEX sv_caller_run_table_id_index ON sv_caller_run_table (id)" );
        } // deconstructor
    }; // class

    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int64_t, // sv_caller_run_id (foreign key)
                                                     int64_t, // read_id (foreign key)
                                                     int64_t, // sort_pos_start
                                                     int64_t, // sort_pos_end
                                                     uint32_t, // from_pos
                                                     uint32_t, // to_pos
                                                     uint32_t, // query_from
                                                     uint32_t, // query_to
                                                     bool, // from_forward @todo save space by compressing booleans?
                                                     bool, // to_forward
                                                     bool // from_seed_start
                                                     >
        TP_SV_JUMP_TABLE;
    class SvJumpTable : public TP_SV_JUMP_TABLE
    {
        std::shared_ptr<CppSQLiteDBExtended> pDatabase;
        CppSQLiteExtQueryStatement<uint32_t> xQuerySize;

      public:
        SvJumpTable( std::shared_ptr<CppSQLiteDBExtended> pDatabase )
            : TP_SV_JUMP_TABLE(
                  *pDatabase, // the database where the table resides
                  "sv_jump_table", // name of the table in the database
                  // column definitions of the table
                  std::vector<std::string>{"sv_caller_run_id", "read_id", "sort_pos_start", "sort_pos_end", "from_pos",
                                           "to_pos", "query_from", "query_to", "from_forward", "to_forward",
                                           "from_seed_start"},
                  // constraints for table
                  std::vector<std::string>{"FOREIGN KEY (sv_caller_run_id) REFERENCES sv_caller_run_table(id)",
                                           "FOREIGN KEY (read_id) REFERENCES read_table(id)"} ),
              pDatabase( pDatabase ),
              xQuerySize( *pDatabase, "SELECT COUNT(*) FROM sv_jump_table" )
        {} // default constructor

        ~SvJumpTable( )
        {
            if( pDatabase->eDatabaseOpeningMode == eCREATE_DB )
            {
                // https://www.sqlite.org/queryplanner.html -> 3.2. Searching And Sorting With A Covering Index
                // index intended for the sweep over the start of all sv-rectangles
                pDatabase->execDML(
                    "CREATE INDEX sv_jump_table_sort_index_start ON sv_jump_table"
                    "(sv_caller_run_id, sort_pos_start, from_pos, to_pos, query_from, query_to, from_forward,"
                    " to_forward, from_seed_start)" );
                // index intended for the sweep over the end of all sv-rectangles
                pDatabase->execDML(
                    "CREATE INDEX sv_jump_table_sort_index_end ON sv_jump_table"
                    "(sv_caller_run_id, sort_pos_end, from_pos, to_pos, query_from, query_to, from_forward,"
                    " to_forward, from_seed_start)" );
            } // if
        } // deconstructor

        inline uint32_t numJumps( )
        {
            return xQuerySize.scalar( );
        } // method
    }; // class


    std::shared_ptr<CppSQLiteDBExtended> pDatabase;
    std::shared_ptr<SequencerTable> pSequencerTable;
    std::shared_ptr<ReadTable> pReadTable;
    std::shared_ptr<PairedReadTable> pPairedReadTable;
    std::shared_ptr<SvCallerRunTable> pSvCallerRunTable;
    std::shared_ptr<SvJumpTable> pSvJumpTable;

    friend class NucSeqFromSql;
    friend class AllNucSeqFromSql;
    friend class PairedNucSeqFromSql;
    friend class PairedReadTable;
    friend class SortedSvJumpFromSql;

  public:
    SV_DB( std::string sName, enumSQLite3DBOpenMode xMode )
        : pDatabase( std::make_shared<CppSQLiteDBExtended>( "", sName, xMode ) ),
          pSequencerTable( std::make_shared<SequencerTable>( pDatabase ) ),
          pReadTable( std::make_shared<ReadTable>( pDatabase ) ),
          pPairedReadTable( std::make_shared<PairedReadTable>( pDatabase, pReadTable ) ),
          pSvCallerRunTable( std::make_shared<SvCallerRunTable>( pDatabase ) ),
          pSvJumpTable( std::make_shared<SvJumpTable>( pDatabase ) )
    {} // constructor

    SV_DB( std::string sName ) : SV_DB( sName, eCREATE_DB )
    {} // constructor

    SV_DB( std::string sName, std::string sMode ) : SV_DB( sName, sMode == "create" ? eCREATE_DB : eOPEN_DB )
    {} // constructor

    class ReadInserter
    {
      private:
        // this is here so that it gets destructed after the transaction context
        std::shared_ptr<SV_DB> pDB;
        // must be after the DB so that it is deconstructed first
        CppSQLiteExtImmediateTransactionContext xTransactionContext;

        int64_t uiSequencerId;

      public:
        ReadInserter( std::shared_ptr<SV_DB> pDB, std::string sSequencerName )
            : pDB( pDB ),
              xTransactionContext( *pDB->pDatabase ),
              uiSequencerId( pDB->pSequencerTable->insertSequencer( sSequencerName ) )
        {} // constructor

        ReadInserter( const ReadInserter& rOther ) = delete; // delete copy constructor

        inline void insertRead( std::shared_ptr<NucSeq> pRead )
        {
            pDB->pReadTable->insertRead( uiSequencerId, pRead );
        } // method

        inline void insertPairedRead( std::shared_ptr<NucSeq> pReadA, std::shared_ptr<NucSeq> pReadB )
        {
            pDB->pPairedReadTable->insertRead( uiSequencerId, pReadA, pReadB );
        } // method
    }; // class

    class SvJumpInserter
    {
        // this is here so that it gets destructed after the transaction context
        std::shared_ptr<SV_DB> pDB;
        // must be after the DB so that it is deconstructed first
        CppSQLiteExtImmediateTransactionContext xTransactionContext;

      public:
        const int64_t iSvCallerRunId;

        class ReadContex
        {
          private:
            std::shared_ptr<SvJumpTable> pSvJumpTable;
            const int64_t iSvCallerRunId;
            const int64_t iReadId;

          public:
            ReadContex( std::shared_ptr<SvJumpTable> pSvJumpTable, const int64_t iSvCallerRunId, const int64_t iReadId )
                : pSvJumpTable( pSvJumpTable ), iSvCallerRunId( iSvCallerRunId ), iReadId( iReadId )
            {} // constructor

            inline void insertJump( const SvJump& rJump )
            {
                if( rJump.does_switch_strand( ) )
                    assert( rJump.from_start( ) > std::numeric_limits<int64_t>::max( ) / 2 );
                pSvJumpTable->xInsertRow( iSvCallerRunId, iReadId, rJump.from_start( ), rJump.from_end( ),
                                          (uint32_t)rJump.uiFrom, (uint32_t)rJump.uiTo, (uint32_t)rJump.uiQueryFrom,
                                          (uint32_t)rJump.uiQueryTo, rJump.bFromForward, rJump.bToForward,
                                          rJump.bFromSeedStart );
            } // method
        }; // class

        SvJumpInserter( std::shared_ptr<SV_DB> pDB,
                        const std::string& rsSvCallerName,
                        const std::string& rsSvCallerDesc )
            : pDB( pDB ),
              xTransactionContext( *pDB->pDatabase ),
              iSvCallerRunId( pDB->pSvCallerRunTable->xInsertRow( rsSvCallerName, rsSvCallerDesc ) )
        {} // constructor

        SvJumpInserter( const SvJumpInserter& ) = delete; // delete copy constructor

        inline ReadContex insertRead( std::shared_ptr<NucSeq> pRead )
        {
            return ReadContex( pDB->pSvJumpTable, iSvCallerRunId,
                               pDB->pReadTable->insertRead( iSvCallerRunId, pRead ) );
        } // method

        inline ReadContex readContext( int64_t iReadId )
        {
            return ReadContex( pDB->pSvJumpTable, iSvCallerRunId, iReadId );
        } // method

    }; // class

    inline void clearCallsTable( )
    {
        pSvJumpTable->clearTable( );
        pSvCallerRunTable->clearTable( );
    } // method

    inline uint32_t numJumps( )
    {
        return pSvJumpTable->numJumps( );
    } // method

}; // class

class SortedSvJumpFromSql
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t> xQueryStart;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t> xQueryEnd;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t>::Iterator
        xTableIteratorStart;
    CppSQLiteExtQueryStatement<uint32_t, uint32_t, uint32_t, uint32_t, bool, bool, bool, int64_t>::Iterator
        xTableIteratorEnd;

  public:
    SortedSvJumpFromSql( std::shared_ptr<SV_DB> pDb, int64_t iSvCallerRunId )
        : pDb( pDb ),
          xQueryStart( *pDb->pDatabase,
                       "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, "
                       "sort_pos_start "
                       "FROM sv_jump_table "
                       "WHERE sv_caller_run_id == ? "
                       "ORDER BY sort_pos_start" ),
          xQueryEnd(
              *pDb->pDatabase,
              "SELECT from_pos, to_pos, query_from, query_to, from_forward, to_forward, from_seed_start, sort_pos_end "
              "FROM sv_jump_table "
              "WHERE sv_caller_run_id == ? "
              "ORDER BY sort_pos_end" ),
          xTableIteratorStart( xQueryStart.vExecuteAndReturnIterator( iSvCallerRunId ) ),
          xTableIteratorEnd( xQueryEnd.vExecuteAndReturnIterator( iSvCallerRunId ) )
    {} // constructor

    bool hasNextStart( )
    {
        return !xTableIteratorStart.eof( );
    } // method

    bool hasNextEnd( )
    {
        return !xTableIteratorEnd.eof( );
    } // method

    bool nextStartIsSmaller( )
    {
        auto xStartTup = xTableIteratorStart.get( );
        auto xEndTup = xTableIteratorEnd.get( );
        return std::get<7>( xStartTup ) <= std::get<7>( xEndTup );
    } // method

    SvJump getNextStart( )
    {
        assert( hasNextStart( ) );

        auto xTup = xTableIteratorStart.get( );
        xTableIteratorStart.next( );
        return SvJump( std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                       std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ) );
    } // method

    SvJump getNextEnd( )
    {
        assert( hasNextEnd( ) );

        auto xTup = xTableIteratorEnd.get( );
        xTableIteratorEnd.next( );
        return SvJump( std::get<0>( xTup ), std::get<1>( xTup ), std::get<2>( xTup ), std::get<3>( xTup ),
                       std::get<4>( xTup ), std::get<5>( xTup ), std::get<6>( xTup ) );
    } // method

}; // class

class AllNucSeqFromSql : public Module<NucSeq, true>
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t> xQuery;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t>::Iterator xTableIterator;

  public:
    AllNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT read_table.sequence, read_table.id "
                  "FROM read_table " ),
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

class NucSeqFromSql : public Module<NucSeq, true>
{
    std::shared_ptr<SV_DB> pDb;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t> xQuery;
    CppSQLiteExtQueryStatement<NucSeqSql, uint32_t>::Iterator xTableIterator;

  public:
    NucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
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

  public:
    PairedNucSeqFromSql( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb )
        : pDb( pDb ),
          xQuery( *pDb->pDatabase,
                  "SELECT A.sequence, B.sequence, A.id, B.id "
                  "FROM read_table A, read_table B "
                  "INNER JOIN paired_read_table "
                  "ON paired_read_table.first_read == A.id "
                  "AND paired_read_table.second_read == B.id " ),
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
        pRet->back( )->iId = std::get<2>( xTup );
        pRet->push_back( std::get<1>( xTup ).pNucSeq );
        pRet->back( )->iId = std::get<3>( xTup );
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

class SvDbInserter : public Module<Container, false, ContainerVector<SvJump>, NucSeq>
{
    std::shared_ptr<SV_DB> pDb;
    SV_DB::SvJumpInserter xInserter;
    std::mutex xMutex;

  public:
    SvDbInserter( const ParameterSetManager& rParameters, std::shared_ptr<SV_DB> pDb, std::string sRunDesc )
        : pDb( pDb ), xInserter( pDb, "MA-SV", sRunDesc )
    {} // constructor

    std::shared_ptr<Container> execute( std::shared_ptr<ContainerVector<SvJump>> pJumps, std::shared_ptr<NucSeq> pRead )
    {
        std::lock_guard<std::mutex> xGuard( xMutex );

        SV_DB::SvJumpInserter::ReadContex xReadContext = xInserter.insertRead( pRead );
        for( SvJump& rJump : *pJumps )
            xReadContext.insertJump( rJump );
        return std::make_shared<Container>( );
        // end of score for xGuard
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
