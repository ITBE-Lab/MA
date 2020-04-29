

#include "ma/container/minimizer_index.h"
using namespace minimizer;

#ifdef WITH_PYTHON

void exportMinimizerIndex( libMS::SubmoduleOrganizer& xOrganizer )
{
#ifdef WITH_ZLIB
    // export the Index class
    py::class_<Index, libMS::Container, std::shared_ptr<Index>>( xOrganizer.container(), "MinimizerIndex" )
        .def( py::init<const ParameterSetManager&, std::string>( ) )
        .def( py::init<const ParameterSetManager&, std::vector<std::string>, std::vector<std::string>>( ) )
        .def( "dump", &Index::dump )
        .def( "seed", &Index::seed );
#endif
} // function
#endif