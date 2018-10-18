/**
 * @file svDb.h
 * @details
 * The database interface for the structural variant caller
 */

#include "CppSQLite3.h"

class SV_DB : public CppSQLite3DB
{
  private:
    void initTables( )
    {
        if( !tableExists( "Sequencers" ) )
            execDML( "CREATE TABLE Sequencers("
                     "ID INT PRIMARY    id              AUTOINCREMENT,"
                     "TEXT              name            NOT NULL"
                     ");" );
        if( !tableExists( "Reads" ) )
            execDML( "CREATE TABLE Reads("
                     "ID INT PRIMARY    id              AUTOINCREMENT,"
                     "INT               sequencerId     NOT NULL,"
                     "TEXT              sequence        NOT NULL"
                     ");" );
        if( !tableExists( "Socs" ) )
            execDML( "CREATE TABLE Socs("
                     "ID INT PRIMARY    id              AUTOINCREMENT,"
                     "INT               readId          NOT NULL,"
                     "BIGINT            start           NOT NULL,"
                     "BIGINT            end             NOT NULL,"
                     "INT               score           NOT NULL"
                     ");" );
    } // method
  public:
    SV_DB( std::string sName )
    {
        open( sName.c_str( ) );
        initTables( );
    } // constructor

    class ReadInserter: private CppSQLite3Statement
    {

    }; // class

    void addSequencerType()
    {

    } // method

}; // class