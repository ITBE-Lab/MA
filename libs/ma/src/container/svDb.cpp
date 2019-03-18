#include "container/svDb.h"

using namespace libMA;

SV_DB::PairedReadTable::PairedReadTable(
    std::shared_ptr<CppSQLiteDBExtended> pDatabase, std::shared_ptr<ReadTable> pReadTable )
    : TP_PAIRED_READ_TABLE( *pDatabase, // the database where the table resides
                            "paired_read_table", // name of the table in the database
                            // column definitions of the table
                            std::vector<std::string>{"first_read", "second_read"},
                            false, // do not generate automatic primary key...
                            // constraints for table
                            std::vector<std::string>{"FOREIGN KEY (first_read) REFERENCES read_table(id)", //
                                                     "FOREIGN KEY (second_read) REFERENCES read_table(id)", //
                                                     "PRIMARY KEY (first_read, second_read)"} ),
      pDatabase( pDatabase ),
      pReadTable( pReadTable )
{} // default constructor

#ifdef WITH_PYTHON

void exportSoCDbWriter( py::module& rxPyModuleId )
{

    // export the SoCInserter class
    py::class_<SV_DB, std::shared_ptr<SV_DB>>( rxPyModuleId, "SV_DB" ) //
        .def( py::init<std::string, std::string>( ) )
        .def("clear_soc_table" , &SV_DB::clearSocTable)
        .def("has_socs" , &SV_DB::hasSoCs);
    // export the SoCInserter class
    py::class_<SV_DB::ReadInserter, std::shared_ptr<SV_DB::ReadInserter>>( rxPyModuleId, "ReadInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, std::string>( ) )
        .def( "insert_read", &SV_DB::ReadInserter::insertRead )
        .def( "insert_paired_read", &SV_DB::ReadInserter::insertPairedRead );

    py::class_<SV_DB::SoCInserter, std::shared_ptr<SV_DB::SoCInserter>>( rxPyModuleId, "SoCInserter" )
        .def( py::init<std::shared_ptr<SV_DB>>( ) );


    py::class_<SV_DB::SvInserter::LineContex>( rxPyModuleId, "SvInserterLineContext" )
        .def("insert_support", &SV_DB::SvInserter::LineContex::insertSupport);
    // export the SvInserter class
    py::class_<SV_DB::SvInserter, std::shared_ptr<SV_DB::SvInserter>>( rxPyModuleId, "SvInserter" )
        .def( py::init<std::shared_ptr<SV_DB>, std::string, std::string>( ) )
        .def( "insert_sv_line", &SV_DB::SvInserter::insertSvLine );
        //.def( "connect_sv_lines", &SV_DB::SvInserter::connectSvLines );

    // export the SoCDbWriter class
    exportModule<SoCDbWriter, std::shared_ptr<SV_DB::SoCInserter>>( rxPyModuleId, "SoCDbWriter" );

    // export the NucSeqFromSql classes
    exportModule<AllNucSeqFromSql, std::shared_ptr<SV_DB>>( rxPyModuleId, "AllNucSeqFromSql" );
    exportModule<NucSeqFromSql, std::shared_ptr<SV_DB>>( rxPyModuleId, "NucSeqFromSql" );
    exportModule<PairedNucSeqFromSql, std::shared_ptr<SV_DB>>( rxPyModuleId, "PairedNucSeqFromSql" );
} // function
#endif
