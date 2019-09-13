#include "container/minimizer_index.h"
using namespace minimizer;

#ifdef WITH_PYTHON

void exportMinimizerIndex( py::module& rxPyModuleId )
{
    // export the Index class
    py::class_<Index, std::shared_ptr<Index>>( rxPyModuleId, "MinimizerIndex" )
        .def( py::init<const ParameterSetManager&, std::string>( ) )
        .def( py::init<const ParameterSetManager&, std::vector<std::shared_ptr<libMA::NucSeq>>,
                       std::vector<std::string>>( ) )
        .def( "dump", &Index::dump )
        .def( "seed", &Index::seed );
} // function
#endif