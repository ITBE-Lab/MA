#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#include "module/computeCoverage.h"
#include <cstdlib>
#include <iostream>

using namespace libMA;


int main( void )
{
#if 0
    ComputeCoverage::SvCallWrapper xCov( SvCall( 10, 100, 2, 2, false, 10 ) );

    auto pQuery = std::make_shared<NucSeq>( );

    ComputeCoverage::uiMaxDistance = 10;
    ComputeCoverage::uiAllowedOverlap = 1;

    std::vector<Seed> vSeeds = {Seed( 0, 4, 2, true ), Seed( 0, 6, 4, true ), Seed( 0, 6, 8, true ),
                                Seed( 0, 2, 12, true ), Seed( 0, 6, 16, true )};

    for( auto xSeed : vSeeds )
        xCov.addSeed( xSeed, pQuery, ComputeCoverage::uiMaxDistance, ComputeCoverage::uiAllowedOverlap );

    for( size_t uiX : {0, 1, 2, 3} )
        for( auto& rSeed : xCov.avCoverageAnalysis[ uiX ] )
            std::cout << uiX << " Seed(" << rSeed.start( ) << ", " << rSeed.size( ) << ", " << rSeed.start_ref( )
                      << ", " << ( rSeed.bOnForwStrand ? "True )" : " False )" ) << std::endl;


    std::vector<std::pair<size_t, nucSeqIndex>> vCoverageList;
    size_t uiCoverageSum = 0;

    xCov.fill_coverage_list( vCoverageList, uiCoverageSum, 0 );

    std::cout << "left-from " << std::endl;
    for( auto xP : vCoverageList )
        std::cout << xP.first << ", " << xP.second << std::endl;

    vCoverageList.clear();
    std::cout << "right-from " << std::endl;
    xCov.fill_coverage_list( vCoverageList, uiCoverageSum, 2 );
    for( auto xP : vCoverageList )
        std::cout << xP.first << ", " << xP.second << std::endl;

    xCov.compute_coverage();
    std::cout << "coverage: " << xCov.xCall.uiCoverage << std::endl;

    assert(xCov.xCall.uiCoverage == 1);
#endif
    return EXIT_SUCCESS;
} /// main function