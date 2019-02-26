/**
 * @file fileWriter.h
 * @brief Writes alignments to a file.
 * @author Markus Schmidt
 */
#ifndef DB_WRITER_H
#define DB_WRITER_H
#ifdef WITH_POSTGRES

#include "container/alignment.h"
#include "libpq-fe.h"
#include "module/module.h"
#include "util/exception.h"

namespace libMA
{
class DbRunConnection
{
  private:
    PGconn* conn;

  public:
    class DbResult
    {
      public:
        PGresult* res;
        DbResult( )
        {} // constructor

        ~DbResult( )
        {
            PQclear( res );
        } // deconstructor

        const std::string get( int iRow, int iCol )
        {
            assert( iRow < PQntuples( res ) );
            assert( iCol < PQnfields( res ) );
            return std::string( PQgetvalue( res, iRow, iCol ) );
        } // method
    }; // class

    DbResult exec( std::string sTransaction )
    {
        DbResult xRet;
        xRet.res = PQexec( conn, sTransaction.c_str( ) );
        ExecStatusType xState = PQresultStatus( xRet.res );
        if( !( xState == PGRES_COMMAND_OK || xState == PGRES_TUPLES_OK ) )
        {
            std::cerr << PQerrorMessage( conn ) << std::endl;
            throw AnnotatedException( std::string( PQerrorMessage( conn ) ) );
        } // if
        return xRet;
    } // method

    DbRunConnection( std::string sConninfo )
    {
        conn = PQconnectdb( sConninfo.c_str( ) );

        /* Check to see that the backend connection was successfully made */
        if( PQstatus( conn ) != CONNECTION_OK )
            throw AnnotatedException( "Connection to database failed:" + std::string( PQerrorMessage( conn ) ) );

        exec( "BEGIN" );
    } // constructor

    ~DbRunConnection( )
    {
        exec( "COMMIT" );
        /* close the connection to the database and cleanup */
        PQfinish( conn );
    } // deconstructor
}; // class


/**
 * @brief Writes the alignment into a postgres database
 * @todo the query input here should actually be a n tuple input.
 */
class DbWriter : public Module<Container, false, NucSeq, ContainerVector<std::shared_ptr<Alignment>>, Pack>
{
  public:
    int32_t iRunId;
    DbRunConnection xConnection;
    /**
     * @brief creates a new FileWriter.
     * @details
     * if sFileName is "stdout" the writer will output to stdout instead of the file.
     * Otherwise sFileName is used as the filename to write to.
     * The file will be truncated is it already exists.
     */
    DbWriter( const ParameterSetManager& rParameters, std::string sConninfo, int32_t iRunId ) : iRunId( iRunId ), xConnection( sConninfo )
    {} // constructor

    virtual std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<NucSeq> pQuery,
                                                         std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
                                                             pAlignments,
                                                         std::shared_ptr<Pack>
                                                             pPack );

}; // class


/**
 * @brief @todo
 * @todo the query input here should actually be a n tuple input.
 */
class PairedDbWriter
    : public Module<Container, false, NucSeq, NucSeq, ContainerVector<std::shared_ptr<Alignment>>, Pack>
{
  public:
    int32_t iRunId;
    DbRunConnection xConnection;

    /**
     * @brief creates a new FileWriter.
     * @details
     * if sFileName is "stdout" the writer will output to stdout instead of the file.
     * Otherwise sFileName is used as the filename to write to.
     * The file will be truncated is it already exists.
     */
    PairedDbWriter( const ParameterSetManager& rParameters, std::string sConninfo, int32_t iRunId ) : iRunId( iRunId ), xConnection( sConninfo )
    {} // constructor

    virtual std::shared_ptr<Container>
        EXPORTED execute( std::shared_ptr<NucSeq> pQuery1, std::shared_ptr<NucSeq> pQuery2,
                          std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> pAlignments,
                          std::shared_ptr<Pack> pPack );

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
#ifdef WITH_BOOST
void exportDBWriter( );
#else
void exportDBWriter( py::module& rxPyModuleId );
#endif
#endif // WITH_PYTHON

#endif // WITH_POSTGRES


#endif // DB_WRITER_H