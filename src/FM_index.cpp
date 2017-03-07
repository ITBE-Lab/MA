#include "FM_index.h"

void exportFM_index()
{
    //export the FM_index class
	boost::python::class_<
        FM_IndexContainer, 
        boost::python::bases<Container>, 
        std::shared_ptr<FM_IndexContainer>
    >("FM_index")
        .def("load", &FM_IndexContainer::vLoadFM_Index)
        .def("exists", &FM_IndexContainer::packExistsOnFileSystem)
        .staticmethod("exists")
        .def("store", &FM_IndexContainer::vStoreFM_Index);
        //todo export the build function
    

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< std::shared_ptr<FM_IndexContainer>, std::shared_ptr<Container> >(); 
}//function