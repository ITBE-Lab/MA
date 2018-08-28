/** 
 * @file fileWriter.h
 * @brief Writes alignments to a file.
 * @author Markus Schmidt
 */
#ifndef DB_WRITER_H
#define DB_WRITER_H
#ifdef WITH_POSTGRES

#include "module/module.h"
#include "util/exception.h"
#include "container/alignment.h"
#include "libpq-fe.h"

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
                DbResult()
                {} // constructor

                ~DbResult()
                {
                    PQclear(res);
                }// deconstructor

                const std::string get(int iRow, int iCol)
                {
                    assert(iRow < PQntuples(res));
                    assert(iCol < PQnfields(res));
                    return std::string(PQgetvalue(res, iRow, iCol));
                }// method
            };// class

            DbResult exec(std::string sTransaction)
            {
                DbResult xRet;
                xRet.res = PQexec(conn, sTransaction.c_str());
                ExecStatusType xState = PQresultStatus(xRet.res);
                if ( !(xState == PGRES_COMMAND_OK || xState == PGRES_TUPLES_OK) )
                {
                    std::cerr << PQerrorMessage(conn) << std::endl;
                    throw AlignerException(std::string(PQerrorMessage(conn)));
                } // if
                return xRet;
            }// method

            DbRunConnection(std::string sConninfo)
            {
                conn = PQconnectdb(sConninfo.c_str());
                
                /* Check to see that the backend connection was successfully made */
                if (PQstatus(conn) != CONNECTION_OK)
                    throw AlignerException("Connection to database failed:" + std::string(PQerrorMessage(conn)));
            }// constructor

            ~DbRunConnection()
            {
                /* close the connection to the database and cleanup */
                PQfinish(conn);
            }// deconstructor
    };// class

    /**
     * @brief Writes SAM output.
     * @note flushing of the outstream; this must be done in the deconstructor of OutStream
     * 
     */
    class DbWriter: public Module
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
        DbWriter(std::string sConninfo, int32_t iRunId)
            :
            iRunId(iRunId),
            xConnection(sConninfo)
        {}//constructor

        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Nil
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(NucSeq)
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "DbWriter";
        }//function

        std::string getFullDesc() const
        {
            return "DbWriter";
        }//function

    };//class

}//namespace

#ifdef WITH_PYTHON
void exportDBWriter();
#endif // WITH_PYTHON

#endif // WITH_POSTGRES


#endif // DB_WRITER_H