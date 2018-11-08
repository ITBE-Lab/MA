/**
 * @file nucSeq.cpp
 * @author Arne Kutzner
 */
#include "container/nucSeq.h"
#include "util/pybind11.h"
using namespace libMA;

/* The translation table for columns.
 * Translates a single character into a 2-bit compressed code.
 */
const unsigned char NucSeq::xNucleotideTranslationTable[ 256 ] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // A == 0; C == 1;
                                                    // G == 2;
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // T == 3;
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // a == 0; c == 1; g == 2;
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // t == 3;
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4}; // predefined array

#ifdef WITH_PYTHON
#ifdef BOOST_PYTHON
void exportNucSeq( )
{
    // export the nucleotidesequence class
    boost::python::class_<NucSeq, boost::noncopyable, boost::python::bases<Container>,
                          std::shared_ptr<NucSeq>>(
        "NucSeq", "Holds a single nucleotide sequence.\n",
        boost::python::init<const char *>( "arg1: self\n"
                                           "arg2: string to initialize the sequence from\n" ) )
        .def( boost::python::init<const std::string>(
            "arg1: self\n"
            "arg2: string to initialize the sequence from\n" ) )
        .def( boost::python::init<>( "arg1: self\n" ) )
        .def( "at", &NucSeq::charAt )
        .def( "__getitem__", &NucSeq::charAt )
        .def( "append", &NucSeq::vAppend_boost )
        .def( "length", &NucSeq::length )
        .def( "__len__", &NucSeq::length )
        .def( "__str__", &NucSeq::toString )
    //.def(
    //        "reverse",
    //        &NucSeq::vReverse
    //    )
#if WITH_QUALITY
        .def( "quality", &NucSeq::getQuality )
#endif
        .def( "fastaq", &NucSeq::fastaq )
        .def_readwrite( "name", &NucSeq::sName );

    // tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<std::shared_ptr<NucSeq>, std::shared_ptr<Container>>( );

    // register return values of vectors of nucseqs
    boost::python::class_<std::vector<std::shared_ptr<NucSeq>>>( "VecRetNuc" )
        .def( boost::python::vector_indexing_suite<
              std::vector<std::shared_ptr<NucSeq>>,
              /*
               *    true = noproxy this means that the content of the vector is already exposed by
               *    boost python.
               *    if this is kept as false, StripOfConsideration would be exposed a second time.
               *    the two StripOfConsiderations would be different and not intercastable.
               *    => keep this as true
               */
              true>( ) );

} // function
#else
#endif
void exportNucSeq( py::module& rxPyModuleId )
{
    // export the nucleotidesequence class
    py::class_<NucSeq, Container, std::shared_ptr<NucSeq>>( rxPyModuleId, "NucSeq" )
        .def( py::init<>( ) ) // default constructor
        .def( py::init<const char*>( ) )
        .def( py::init<const std::string>( ) )
        .def( "at", &NucSeq::charAt )
        .def( "__getitem__", &NucSeq::charAt )
        .def( "append", &NucSeq::vAppend_boost )
        .def( "length", &NucSeq::length )
        .def( "__len__", &NucSeq::length )
        .def( "__str__", &NucSeq::toString )
#if WITH_QUALITY
        .def( "quality", &NucSeq::getQuality )
#endif
        .def( "fastaq", &NucSeq::fastaq )
        .def_readwrite( "name", &NucSeq::sName );

    // register return values of vectors of nucseqs
    py::bind_vector<std::vector<std::shared_ptr<NucSeq>>>( rxPyModuleId, "VecRetNuc" );
} // function
#endif