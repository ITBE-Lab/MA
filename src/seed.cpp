#include "seed.h"
using namespace libLAuS;

void exportSeed()
{

    exportInterval<nucSeqIndex>();

    //export the Seed class
    boost::python::class_<
            Seed,
            boost::python::bases<Interval<nucSeqIndex>>
        >("Seed")
        .def_readwrite("start_ref", &Seed::uiPosOnReference)
    ;

    //export the Seeds class
    boost::python::class_<
        Seeds, 
        boost::python::bases<Container>, 
        boost::python::bases<std::list<Seed>>, 
        std::shared_ptr<Seeds>
    >(
        "Seeds"
    )
    .def(boost::python::init<std::shared_ptr<Seeds>>())
    .def(boost::python::vector_indexing_suite<
            Seeds,
            /*
            *	true = noproxy this means that the content of the vector is already exposed by
            *	boost python. 
            *	if this is kept as false, Container would be exposed a second time.
            *	the two Containers would be different and not inter castable.
            */
            true
        >());

    //make vectors of container-pointers a thing
    iterable_converter()
        .from_python<Seeds>();
    
    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<
            std::shared_ptr<Seeds>,
            std::shared_ptr<Container>
        >(); 
}//function