#include "FM_index.h"

void exportFM_index()
{
    //export the FM_index class
	boost::python::class_<FM_Index, boost::python::bases<Container>>("FM_index")
        .def("load", &FM_Index::vLoadFM_indexWrapper)
        .def("exists", &FM_Index::packExistsOnFileSystem)
        .def("store", &FM_Index::vStoreFM_IndexWrapper);
        //todo export the build function
    



}