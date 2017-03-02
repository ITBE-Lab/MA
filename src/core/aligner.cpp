#include "core/aligner.h"

void exportAligner()
{
        boost::python::class_<Aligner>("Aligner")
                .def("setData", &Aligner::setData)
                .def("addModule", &Aligner::addModule)
                .def("step", &Aligner::step)
                .def("steps", &Aligner::steps)
        ;
}


BOOST_PYTHON_MODULE(aligner)
{
	exportAligner();
        exportModule();
}