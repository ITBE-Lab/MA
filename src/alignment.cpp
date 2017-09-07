#include "alignment.h"

void exportAlignment()
{
    boost::python::class_<
            Alignment, 
            boost::python::bases<Container>, 
            std::shared_ptr<Alignment>
        >(
                "Alignment",
                "contains the final output of the aligner\n"
            );

        //tell boost python that pointers of these classes can be converted implicitly
        boost::python::implicitly_convertible< 
            std::shared_ptr<Alignment>,
            std::shared_ptr<Container> 
        >();

}