#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "MySQL_con.h" // MySQL connector
#include "common.h"
#include "container/sv_db/svSchema.h"
#include "db_con_pool.h"
#include <cstdlib>
#include <iostream>


int main( void )
{
    std::vector<std::future<void>> vFutures;
    {
        doNoExcept( [&] {
            SQLDBConPool<MySQLConDB> xDBPool( 32, "Pooled_DB" );

            for( int i = 0; i < 32; i++ )
                // type behind auto: std::shared_ptr<SQLDBConPool<MySQLConDB>::PooledSQLDBCon>
                vFutures.push_back( xDBPool.enqueue( []( auto pDBCon ) {
                    doNoExcept(
                        [&] {
                            using ConnectionType = typename decltype( *pDBCon );
                            libMA::SV_Schema<ConnectionType> xSB_DB_SchemaView( pDBCon );
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


    return EXIT_SUCCESS;
} /// main function
