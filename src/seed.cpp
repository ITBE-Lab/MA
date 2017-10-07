#include "seed.h"

void exportSeed()
{

    //export the Seed class
    boost::python::class_<
            Seed
        >(
            "Seed",
            "A single seed.\n",
            boost::python::no_init
        )
        .def(
                "start", 
                &Seed::start_boost1,
                "arg1: self\n"
                "returns: the start position on the query\n"
            )
        .def(
                "end", 
                &Seed::end_boost1,
                "arg1: self\n"
                "returns: the end position on the query\n"
            )
        .def(
                "start_ref", 
                &Seed::start_ref,
                "arg1: self\n"
                "returns: the start position on the reference\n"
            )
        .def(
                "end_ref", 
                &Seed::end_ref,
                "arg1: self\n"
                "returns: the end position on the reference\n"
            )
        .def(
                "size", 
                &Seed::size_boost1,
                "arg1: self\n"
                "returns: the size of the Seed\n"
            )
    ;

    //export the SeedIter class
    boost::python::class_<
            SeedIter
        >(
            "SeedIter",
            "Iterator for Seeds.\n",
            boost::python::no_init
        )
        .def(
                "__next__", 
                &SeedIter::next_boost,
                "arg1: self\n"
                "returns: the current element\n"
                "\n"
                "Return the current element.\n"
                "Move the iterator to the next element of the list.\n"
            )
    ;

    //export the Seeds class
    boost::python::class_<
        Seeds, 
        boost::python::bases<Container>, 
        boost::python::bases<std::list<Seed>>, 
        std::shared_ptr<Seeds>
    >(
        "Seeds",
        "Contains multiple strips of consideration.\n"
    )
    .def(
            "__len__",
            &Seeds::size,
            "arg1: self\n"
            "returns: list size\n"
        )
    .def(
            "__iter__",
            &Seeds::boost_python_iterator,
            "arg1: self\n"
            "returns: list size\n"
        )
    .def(
            "get_score",
            &Seeds::getScore,
            "arg1: self\n"
            "returns: sum off all scores within the list\n"
        )
    .def(
            "append",
            &Seeds::append,
            "arg1: self\n"
            "arg1: other Seed\n"
            "\n"
            "Appends a copy of the other list to this list.\n"
        );
    
    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<
            std::shared_ptr<Seeds>,
            std::shared_ptr<Container>
        >(); 

    //export the SeedsVector class
    boost::python::class_<
        SeedsVector,
        boost::python::bases<Container>,
        boost::python::bases<std::vector<std::shared_ptr<Seeds>>>,
        std::shared_ptr<SeedsVector>
    >(
        "SeedsVector",
        "Contains multiple seed lists.\n"
    )
    .def(boost::python::vector_indexing_suite<
            SeedsVector,
            /*
            *	true = noproxy this means that the content of the vector is already exposed by
            *	boost python. 
            *	if this is kept as false, SeedsVector would be exposed a second time.
            *	the two SeedsVector would be different and not inter castable.
            */
            true
        >());

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible< 
        std::shared_ptr<SeedsVector>, 
        std::shared_ptr<Container>
    >(); 
}//function