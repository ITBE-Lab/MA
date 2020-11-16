/**
 * @file nucSeq.cpp
 * @author Arne Kutzner
 */
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#endif

#include "ma/container/nucSeq.h"
#include "ms/util/pybind11.h"
using namespace libMA;

/* The translation table for columns.
 * Translates a single character into a 2-bit compressed code.
 */
const unsigned char NucSeq::xNucleotideTranslationTable[ 256 ] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4,
    4, 4, // A == 0; C == 1;
          // G == 2;
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // T == 3;
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // a == 0; c == 1; g == 2;
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // t == 3;
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}; // predefined array

#ifdef WITH_PYTHON
void exportNucSeq( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the nucleotidesequence class
    py::class_<NucSeq, libMS::Container, std::shared_ptr<NucSeq>>( xOrganizer.container( ), "NucSeq" )
        .def( py::init<>( ) ) // default constructor
        .def( py::init<const char*>( ) )
        .def( py::init<const std::string>( ) )
        .def( py::init<std::shared_ptr<NucSeq>, nucSeqIndex, nucSeqIndex>( ) )
        .def( "at", &NucSeq::charAt )
        .def( "__getitem__", &NucSeq::charAt )
        .def( "append", &NucSeq::vAppend_boost )
        .def( "length", &NucSeq::length )
        .def( "complement", &NucSeq::vSwitchAllBasePairsToComplement )
        .def( "reverse", &NucSeq::vReverseAll )
        .def( "__len__", &NucSeq::length )
        .def( "__str__", &NucSeq::toString )
#if WITH_QUALITY
        .def( "quality", &NucSeq::getQuality )
#endif
        .def( "fastaq", &NucSeq::fastaq )
        .def( "equals", &NucSeq::isEqual )
        .def( "check", &NucSeq::check )
        .def_readwrite( "name", &NucSeq::sName )
        .def_readwrite( "id", &NucSeq::iId );

    // register return values of vectors of nucseqs
    py::bind_vector<std::vector<std::shared_ptr<NucSeq>>>( xOrganizer.container( ), "VecRetNuc", "docstr" )
        .def( py::init<>( ) );

    // export the NucSeqSql class
    py::class_<NucSeqSql, std::shared_ptr<NucSeqSql>>( xOrganizer._container( ), "NucSeqSql" )
        .def( py::init<>( ) ) // default constructor
        .def( "fromBlob", &NucSeqSql::fromBlob )
        .def_readwrite( "seq", &NucSeqSql::pNucSeq );
} // function
#endif
