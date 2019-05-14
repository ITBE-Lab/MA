#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "container/container.h"
#include "module/module.h"
#include <cstdlib>
#include <iostream>

using namespace libMA;


int main( void )
{
    const ParameterSetManager xParams;
    auto pTest = std::make_shared<Test>( xParams );
    std::vector<std::shared_ptr<BasePledge>> vGraphSinks;
    auto pPromise = promiseMe( pTest );
    for( unsigned int i = 0; i < 32; i++ )
    {
        vGraphSinks.push_back( pPromise );
    } // for
    BasePledge::simultaneousGet( vGraphSinks );

    return EXIT_SUCCESS;
} /// main function