#include <boost/python.hpp>
#include <Python.h>
#include <iostream>

char const* greetPython()
{
	return "Hi, Python!";

}

BOOST_PYTHON_MODULE(test)
{
	boost::python::def("greetPython", greetPython);
}