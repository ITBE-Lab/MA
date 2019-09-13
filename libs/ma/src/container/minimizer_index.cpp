#include "container/minimizer_index.h"
using namespace minimizer;

#ifdef WITH_PYTHON
/**
 * @brief function used for time measuring in python code
 */
void create_minimizer_index(std::string sGenome, std::string sIndexPrefix)
{
    auto pParameter = std::make_shared<ParameterSetManager>();
    pParameter->pGlobalParameterSet->piNumberOfThreads->set(1);
    pParameter->pGlobalParameterSet->pbUseMaxHardareConcurrency->set(false);
    Index xIndex(pParameter, sGenome);
    xIndex.dump(sIndexPrefix);
} // method

void exportMinimizerIndex( py::module& rxPyModuleId )
{
    // export the Index class
    py::class_<Index, std::shared_ptr<Index>>( rxPyModuleId, "MinimizerIndex" )
        .def( py::init<std::shared_ptr<ParameterSetManager>, std::string>( ) )
        .def( "dump", &Index::dump );
    
    rxPyModuleId.def("create_minimizer_index", &create_minimizer_index);
} // function
#endif