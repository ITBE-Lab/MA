#include "ma/module/sam_reader.h"
#include "ms/util/pybind11.h"

using namespace libMA;
using namespace libMS;

#ifdef WITH_PYTHON

void exportSamFileReader( libMS::SubmoduleOrganizer& xOrganizer )
{
    py::class_<ReadByName, libMS::Container, std::shared_ptr<ReadByName>>( xOrganizer.container( ), "ReadByName" )
        .def( py::init<>( ) ) // default constructor
        .def( "append", &ReadByName::append )
        .def( "__iter__", &ReadByName::iter );
    py::class_<SeedsByName, libMS::Container, std::shared_ptr<SeedsByName>>( xOrganizer.container( ), "SeedsByName" )
        .def( py::init<>( ) ) // default constructor
        .def( "append", &SeedsByName::append )
        .def( "mergeAll", &SeedsByName::mergeAll );
    exportModule<SamFileReader>( xOrganizer, "SamFileReader" );
    exportModule<GetSeedsByName>( xOrganizer, "GetSeedsByName" );
    exportModule<GetSeedSetCompByName>( xOrganizer, "GetSeedSetCompByName" );
    exportModule<GetSeedsByReadName>( xOrganizer, "GetSeedsByReadName" );
    exportModule<GetReadByName>( xOrganizer, "GetReadByName" );
} // function

#endif