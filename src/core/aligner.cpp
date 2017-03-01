#include "core/aligner.h"


BOOST_PYTHON_MODULE(aligner)
{
	boost::python::class_<Aligner>("Aligner")
        .def("setData", &Aligner::setData)
        .def("addModule", &Aligner::addModule)
        .def("step", &Aligner::step)
        .def("steps", &Aligner::steps)
        ;
}