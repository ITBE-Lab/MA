#include "container/minimizersHash.h"
#include "module/minimizers.h"
using namespace libMA;


void exportMinimizersHash()
{
    
    boost::python::class_<
            MinimizersHash<Minimizers::w,Minimizers::k>, 
            boost::noncopyable,
            boost::python::bases<Container>, 
            std::shared_ptr<MinimizersHash<Minimizers::w,Minimizers::k>>
        >(
                "MinimizersHash"
            )
        .def("to_file", &MinimizersHash<Minimizers::w,Minimizers::k>::toFile)
        .def("from_file", &MinimizersHash<Minimizers::w,Minimizers::k>::fromFile)
        .staticmethod("from_file")
        .def("key_len", &MinimizersHash<Minimizers::w,Minimizers::k>::keyLen)
    ;

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible< 
        std::shared_ptr<MinimizersHash<Minimizers::w,Minimizers::k>>,
        std::shared_ptr<Container> 
    >();
    
    boost::python::class_<
            MinimizersVector<Minimizers::w,Minimizers::k>, 
            boost::python::bases<Container>, 
            std::shared_ptr<MinimizersVector<Minimizers::w,Minimizers::k>>
        >(
                "MinimizersVector"
            )
        .def("toHash", &MinimizersVector<Minimizers::w,Minimizers::k>::toHash)
    ;

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible< 
        std::shared_ptr<MinimizersVector<Minimizers::w,Minimizers::k>>,
        std::shared_ptr<Container> 
    >();

}