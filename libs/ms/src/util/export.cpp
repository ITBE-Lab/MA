/**
 * @file export.cpp
 * @author Markus Schmidt
 */
#include "util/export.h"
#include "util/execution-context.h"
#include "util/parameter.h"

using namespace libMS;


#ifdef WITH_PYTHON

/*
 * The exposing of the ParameterSetManager and it's components has been reworked in SV branch
 * -> here only some quickfixes, so that python does not immedeately segfault if someone tries to use this...
 */

template <typename TP_VALUE> void exportAlignerParameter( py::module& rxPyModuleId, std::string sName )
{
    py::class_< //
        AlignerParameter<TP_VALUE>, //
        AlignerParameterBase, //
        std::shared_ptr<AlignerParameter<TP_VALUE>> //
        >( rxPyModuleId, sName.c_str( ) ) //
        .def( py::init<const std::string,
                       const std::string,
                       const std::string,
                       const size_t,
                       const std::string,
                       const TP_VALUE>( ) )
        .def( "set", &AlignerParameter<TP_VALUE>::set ) //
        .def( "get", &AlignerParameter<TP_VALUE>::get_py );
} // function

/**
 * Parameter set manager export is here, so that we do not need to include pybind11 in parameter.h
 */
void exportParameter( py::module& rxPyModuleId )
{
    py::class_<AlignerParameterBase, std::shared_ptr<AlignerParameterBase>>( rxPyModuleId, "AlignerParameterBase" ) //
        .def( "mirror", &AlignerParameterBase::mirror ) //
        .def_readonly( "name", &AlignerParameterBase::sName ) //
        .def_readonly( "description", &AlignerParameterBase::sDescription );

    exportAlignerParameter<int>( rxPyModuleId, "AlignerParameterInt" );
    exportAlignerParameter<short>( rxPyModuleId, "AlignerParameterShort" );
    exportAlignerParameter<bool>( rxPyModuleId, "AlignerParameterBool" );
    exportAlignerParameter<double>( rxPyModuleId, "AlignerParameterDouble" );
    exportAlignerParameter<float>( rxPyModuleId, "AlignerParameterFloat" );
    py::class_<AlignerParameter<AlignerParameterBase::ChoicesType>,
               AlignerParameterBase,
               std::shared_ptr<AlignerParameter<AlignerParameterBase::ChoicesType>>>( rxPyModuleId,
                                                                                      "AlignerParameterChoice" ) //
        .def( "set", &AlignerParameter<AlignerParameterBase::ChoicesType>::set ) //
        .def( "get", &AlignerParameter<AlignerParameterBase::ChoicesType>::get_py );
    // exportAlignerParameter<fs::path>( rxPyModuleId, "AlignerParameterFilePath" ); @todo

    // Export ParameterSetBase Class
    py::class_<ParameterSetBase, std::shared_ptr<ParameterSetBase>>( rxPyModuleId, "ParameterSetBase" ) //
        .def( "mirror", &ParameterSetBase::mirror ) //
        .def( "register_parameter", &ParameterSetBase::registerParameter ) //
        .def( "unregister_parameter", &ParameterSetBase::unregisterParameter ) //
        .def( "by_name", &ParameterSetBase::byName ) //
        .def( "by_short", &ParameterSetBase::byShort );


    // Export Presetting Class
    py::class_<Presetting, ParameterSetBase, std::shared_ptr<Presetting>>( rxPyModuleId, "Presetting" ) //
        .def( py::init<std::string>( ) ) //
        .def_readonly( "name", &Presetting::sName );

    // Export Presetting Class
    py::class_<GeneralParameter, ParameterSetBase, std::shared_ptr<GeneralParameter>>( rxPyModuleId,
                                                                                       "GeneralSettings" ) //
        .def( py::init<>( ) );

    // Export ParameterSetManager Class
    py::class_<ParameterSetManager>( rxPyModuleId, "ParameterSetManager" ) //
        .def( py::init<>( ) ) //
        .def_readwrite( "global_settings", &ParameterSetManager::pGeneralParameterSet ) //
        .def( "get", &ParameterSetManager::get )
        .def( "add_setting", &ParameterSetManager::addSetting )
        .def( "get_num_threads", &ParameterSetManager::getNumThreads )
        .def( "set_selected", &ParameterSetManager::setSelected )
        .def( "get_selected", &ParameterSetManager::getSelected_py )
        .def( "by_name", &ParameterSetManager::byName )
        .def( "by_short", &ParameterSetManager::byShort );
} // function

void exportExecutionContext( py::module& rxPyModuleId )
{
    py::class_<GenomeManager>( rxPyModuleId, "GenomeManager" ) //
        .def( "load_genome", &GenomeManager::loadGenome );
    py::class_<ReadsManager>( rxPyModuleId, "ReadsManager" ) //
        .def_readwrite( "primary_queries", &ReadsManager::vsPrimaryQueryFullFileName ) //
        .def_readwrite( "mate_queries", &ReadsManager::vsMateQueryFullFileName );
    py::class_<OutputManager>( rxPyModuleId, "OutputManager" );

    py::class_<ExecutionContext>( rxPyModuleId, "ExecutionContext" ) //
        .def( py::init<>( ) ) //
        .def( "do_align", &ExecutionContext::doAlignCallbackLess ) //
        // def_readonly required since xParameterSetManager is non-copyable
        .def_readonly( "parameter_set_manager", &ExecutionContext::xParameterSetManager ) //
        .def_readwrite( "genome_manager", &ExecutionContext::xGenomeManager ) //
        .def_readwrite( "reads_manager", &ExecutionContext::xReadsManager ) //
        .def_readwrite( "output_manager", &ExecutionContext::xOutputManager );
} // function

PYBIND11_MODULE( libMS, libMsModule )
{
    DEBUG_3( std::cout.setf( std::ios::unitbuf ); )
    exportParameter( libMsModule );
    exportContainer( libMsModule );
    exportModuleClass( libMsModule );
    exportPoolContainer( libMsModule );
} // function

#endif
