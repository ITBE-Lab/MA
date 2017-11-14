#include "nmwMultiple.h"

void exportNmwMultiple()
{
    //export the NmwMultiple class
	boost::python::class_<
            NmwMultiple, 
            boost::python::bases<CppModule>,
            std::shared_ptr<NmwMultiple>
        >(
        "NmwMultiple"
    )
        .def_readwrite("try_n_many", &NmwMultiple::tryNmany)
    ;
	boost::python::implicitly_convertible< 
		std::shared_ptr<NmwMultiple>,
		std::shared_ptr<CppModule> 
	>();
}//function