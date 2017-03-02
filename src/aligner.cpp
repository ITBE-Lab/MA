#include "aligner.h"

//function called in order to export this module
void exportAligner()
{
        boost::python::class_<Aligner>("Aligner")
                .def("setData", &Aligner::setData)
                .def("addModule", &Aligner::addModule)
                .def("step", &Aligner::step)
                .def("steps", &Aligner::steps)
        ;
}

//this creates our main function
BOOST_PYTHON_MODULE(aligner)
{
	exportAligner();
        exportModule();
}