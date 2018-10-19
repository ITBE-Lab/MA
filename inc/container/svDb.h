/**
 * @file svDb.h
 * @details
 * The database interface for the structural variant caller
 */

#include "util/sqlite3.h"

class SV_DB : public CppSQLite3DB
{
  private:
    CppSQLiteDBExtended xDatabase;


    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<std::string, // sequencer name
                                                     std::string // parameter setting (for the moment merely the mode)
                                                     >
        TP_SEQUENCER_TABLE;
    class SequencerTable : public TP_SEQUENCER_TABLE
    {
      public:
        SequencerTable( CppSQLiteDBExtended& rDb )
            : TP_SEQUENCER_TABLE( rDb, // the database where the table resides
                                  "sequencer_table", // name of the table in the database
                                  std::vector<std::string>{"name", "mode"} // column definitions of the table
              )
        {} // default constructor
    }; // class


    typedef CppSQLiteExtTableWithAutomaticPrimaryKey<int, // sequencer id
                                                     std::string, // read name
                                                     NucSeq // read sequence
                                                     >
        TP_READ_TABLE;
    class ReadTable : public TP_READ_TABLE
    {
      public:
        ReadTable( CppSQLiteDBExtended& rDb )
            : TP_READ_TABLE( rDb, // the database where the table resides
                             "sequencer_table", // name of the table in the database
                             std::vector<std::string>{"name", "mode"} // column definitions of the table
              )
        {} // default constructor
    }; // class
    CppSQLiteExtTable<std::string // read Name
                      >
        x;
    CppSQLiteExtTable<std::string // read Name
                      >
        xReadTable;

  public:
    SV_DB( std::string sName ) : xDatabase( "./", sName, eOPEN_DB )

    {

        // /* We create the insert statement for the table alignment results.
        //  */
        // CppSQLiteExtInsertStatement<int, // 1. analysis_id INTEGER
        //                             const char*, // 2. sequence_name TEXT
        //                             >
        //     xSQLInsertAlignmentResults( rxDatabase, "alignment_results" );
        //
        // xTestDB.vCreateTable( std::vector<std::string>{"number INTEGER", // Foreign Key The gene_id of the Gene
        //                                                "text TEXT"}, // vector initialization
        //                       "table1" // sTableName
        // );

        CppSQLiteExtTable<int, std::string> table1(
            xTestDB, // the database where the table resides
            "table1", // name of the table in the database
            std::vector<std::string>{"number", "text"}, // column definitions of the table
            false // true == we build automatically a column for the primary key
        );
        {
            CppSQLiteExtImmediateTransactionContext xTransactionContext( xTestDB );
            std::cout << "START INSERT" << std::endl;
            table1.xInsertRow( 12, "abc " );
            table1.xInsertRow( 14, "gfd " );
            std::cout << "END INSERT" << std::endl;
        }
        table1.vDump( std::cout );
        table1.vForAllTableRowsUnpackedDo( []( int i, std::string s ) { std::cout << s << std::endl; } );

        CppSQLiteExtQueryStatement<int> statement( xTestDB,
                                                   std::string( "SELECT number FROM " ).append( "table1" ).c_str( ) );
        statement.vExecuteAndForAllRowsUnpackedDo( []( int i ) { std::cout << i << std::endl; } );
    } // constructor

    SV_DB( const SV_DB& rOther ) = delete; // delete copy constructor

    class ReadInserter : private CppSQLite3Statement
    {}; // class

    void addSequencerType( )
    {} // method

}; // class