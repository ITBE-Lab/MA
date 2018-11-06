#include "container/svDb.h"

using namespace libMA;

#ifdef WITH_PYTHON

#ifdef BOOST_PYTHON
void exportSoCDbWriter( )
{

    // export the SoCInserter class
    boost::python::class_<SV_DB, std::shared_ptr<SV_DB>, boost::noncopyable>( "SV_DB",
                                                                              boost::python::init<std::string, std::string>( ) );
    // export the SoCInserter class
    boost::python::class_<SV_DB::SoCInserter, std::shared_ptr<SV_DB::SoCInserter>, boost::noncopyable>(
        "SoCInserter", boost::python::init<std::shared_ptr<SV_DB>, std::string>( ) );

    // export the SoCDbWriter class
    exportModule<SoCDbWriter, std::shared_ptr<SV_DB::SoCInserter>>( "SoCDbWriter" );

    // export the NucSeqFromSql class
    exportModule<NucSeqFromSql, std::shared_ptr<SV_DB>, std::string>( "NucSeqFromSql" );
} // function
#endif
#endif
