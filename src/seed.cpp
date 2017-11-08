#include "seed.h"


class SeedIter{
private:
    std::list<Seed>::iterator it;
    std::list<Seed>::iterator end;

public:
    SeedIter(std::list<Seed>* pList)
            :
        it(pList->begin()),
        end(pList->end())
    {}//constructor
    
    Seed next_boost()
    {
        if(it == end)
        {
            PyErr_SetNone(PyExc_StopIteration);
            boost::python::throw_error_already_set();
            return Seed(0,0,0);
        }//if
        Seed xRet = *it;
        ++it;
        return xRet;
    }//function
};//class


SeedIter boost_python_iterator(std::shared_ptr<Seeds> seeds)
{
    return SeedIter(seeds.get());
}//function

void exportSeed()
{

    //export the Seed class
    boost::python::class_<
            Seed
        >(
            "Seed",
            "A single seed.\n",
            boost::python::init<nucSeqIndex, nucSeqIndex, nucSeqIndex>()
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
                boost::python::with_custodian_and_ward_postcall<0,1>(),
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
    .def(boost::python::init<std::shared_ptr<Seeds>>())
    .def(
            "__len__",
            &Seeds::size,
            "arg1: self\n"
            "returns: list size\n"
        )
    .def(
            "__iter__",
            &boost_python_iterator,
            boost::python::with_custodian_and_ward_postcall<0,1>(),
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
            "arg1: other Seeds\n"
            "\n"
            "Appends a copy of the other list to this list.\n"
        )
    .def(
            "append",
            &Seeds::push_back,
            "arg1: self\n"
            "arg1: a Seed\n"
            "\n"
            "Appends a seed to this list.\n"
        );
    //boost::python::register_ptr_to_python< std::shared_ptr<Seeds> >();
    
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