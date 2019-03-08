#include "container/svDb.h"

using namespace libMA;

#ifdef WITH_PYTHON

#ifdef BOOST_PYTHON
void exportSoCDbWriter( )
{

    // export the SoCInserter class
    boost::python::class_<SV_DB, std::shared_ptr<SV_DB>, boost::noncopyable>(
        "SV_DB", boost::python::init<std::string, std::string>( ) );
    // export the SoCInserter class
    boost::python::class_<SV_DB::SoCInserter, std::shared_ptr<SV_DB::SoCInserter>, boost::noncopyable>(
        "SoCInserter", boost::python::init<std::shared_ptr<SV_DB>, std::string>( ) );

    // export the SoCDbWriter class
    exportModule<SoCDbWriter, std::shared_ptr<SV_DB::SoCInserter>>( "SoCDbWriter" );

    // export the NucSeqFromSql class
    exportModule<NucSeqFromSql, std::shared_ptr<SV_DB>, std::string>( "NucSeqFromSql" );
} // function
#else
void exportSoCDbWriter( py::module& rxPyModuleId )
{

    // export the SoCInserter class
    py::class_<SV_DB, std::shared_ptr<SV_DB>>( rxPyModuleId, "SV_DB" ) //
        .def( py::init<std::string, std::string>( ) );
    // export the SoCInserter class
    py::class_<SV_DB::SoCInserter, std::shared_ptr<SV_DB::SoCInserter>>( rxPyModuleId, "SoCInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, std::string>( ) )
        .def("insert_read", &SV_DB::SoCInserter::insertRead )
        .def("insert_paired_read", &SV_DB::SoCInserter::insertPairedRead );

    // export the SoCDbWriter class
    exportModule<SoCDbWriter, std::shared_ptr<SV_DB::SoCInserter>>( rxPyModuleId, "SoCDbWriter" );

    // export the NucSeqFromSql class
    exportModule<NucSeqFromSql, std::shared_ptr<SV_DB>, std::string>( rxPyModuleId, "NucSeqFromSql" );
    exportModule<PairedNucSeqFromSql, std::shared_ptr<SV_DB>, size_t, std::string>( rxPyModuleId,
                                                                                    "PairedNucSeqFromSql" );
} // function
#endif
#endif
