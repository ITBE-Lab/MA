#include "getBestOnly.h"

void exportGetBestOnly()
{
    //export the GetBestOnly class
	boost::python::class_<
            GetBestOnly, 
            boost::python::bases<CppModule>,
            std::shared_ptr<GetBestOnly>
        >(
        "GetBestOnly"
    );
	boost::python::implicitly_convertible< 
		std::shared_ptr<GetBestOnly>,
		std::shared_ptr<CppModule> 
	>();
}//function