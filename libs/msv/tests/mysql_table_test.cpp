#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "db_con_pool.h"
#include "msv/container/sv_db/tables/sequencer.h"
#include "mysql_con.h" // MySQL connector
#include "sql_api.h"
#include <cstdlib>
#include <iostream>
#include <type_traits>


int main( void )
{
#ifdef __GNUC__
    std::vector<std::future<void>> vFutures;
    {
        doNoExcept( [&] {
            SQLDBConPool<MySQLConDB> xDBPool( 32, json{{SCHEMA, {{NAME, "pooled_DB"}, {FLAGS, {DROP_ON_CLOSURE}}}}} );

            for( int i = 0; i < 32; i++ )
                // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
                vFutures.push_back( xDBPool.enqueue( []( auto pDBCon ) {
                    doNoExcept(
                        [&] {
                            // pDBCon is of type std::shared_ptr<ConnectionType>
                            using ConnectionType = typename std::remove_reference<decltype( *pDBCon )>::type;

                            libMSV::SequencerTable<ConnectionType> xSequencerTable( pDBCon );

                            std::cout << "Job executed in task: " << pDBCon->getTaskId( ) << std::endl;
                            // checkDB( pDBCon, pDBCon->getTaskId( ) );
                            pDBCon->doPoolSafe( [] { std::cout << "This print is pool safe ..." << std::endl; } );
                        },
                        "Problem during thread execution" );
                } ) );
        } );

    } // close the pool

    // Get all future exception safe
    for( auto& rFurture : vFutures )
        doNoExcept( [&] { rFurture.get( ); } );

    std::cout << "ALL WORK DONE ..." << std::endl;
#endif

    return EXIT_SUCCESS;
} /// main function
