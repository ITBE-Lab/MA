#include "container/svDb.h"

using namespace libMA;

#ifdef WITH_PYTHON
void exportSoCDbWriter( )
{

    // export the SoCInserter class
    boost::python::class_<SV_DB::SoCInserter, std::shared_ptr<SV_DB::SoCInserter>, boost::noncopyable>(
        "SV_DB_SoCInserter", boost::python::no_init );
    // export the SoCInserter class
    boost::python::class_<SV_DB, std::shared_ptr<SV_DB>, boost::noncopyable>( "SV_DB",
                                                                              boost::python::init<std::string>( ) )
        .def( "add_sequencer", &SV_DB::addSequencerType );

    // export the SoCDbWriter class
    exportModule<SoCDbWriter, std::shared_ptr<SV_DB::SoCInserter>>( "SoCDbWriter" );
} // function
#endif