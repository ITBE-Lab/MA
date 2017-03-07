#include "aligner.h"

//function called in order to export this module
void exportAligner()
{
       /* boost::python::class_<Aligner>("Aligner")
                .def("setData", &Aligner::setData)
                .def("addModule", &Aligner::addModule)
                .def("step", &Aligner::step)
                .def("steps", &Aligner::steps)
                .def("stepPossible", &Aligner::stepPossible)
                .def("test", &Aligner::test)
        ;*/
}

//this creates our main function
BOOST_PYTHON_MODULE(LAuS)
{
	exportAligner();
        exportModule();
        exportContainer();
        exportFM_index();
        exportSequence();
        exportSegmentation();
        exportPack();
        exportIntervalTree();
        exportGraphicalMethod();
}