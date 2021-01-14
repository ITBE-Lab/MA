/**
 * @file export.cpp
 * @author Markus Schmidt
 */
#include "ms/container/container.h"
#include "ms/container/sv_db/pool_container.h"
#include "ms/module/module.h"
#include "ms/util/parameter.h"

using namespace libMS;


#ifdef WITH_PYTHON

/*
 * The exposing of the ParameterSetManager and it's components has been reworked in SV branch
 * -> here only some quickfixes, so that python does not immedeately segfault if someone tries to use this...
 */

template <typename TP_VALUE> void exportAlignerParameter( SubmoduleOrganizer& xOrganizer, std::string sName )
{
    py::class_< //
        AlignerParameter<TP_VALUE>, //
        AlignerParameterBase, //
        std::shared_ptr<AlignerParameter<TP_VALUE>> //
        >( xOrganizer.util( ), sName.c_str( ) ) //
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
void exportParameter( SubmoduleOrganizer& xOrganizer )
{
    py::class_<AlignerParameterBase, std::shared_ptr<AlignerParameterBase>>( xOrganizer._util( ),
                                                                             "AlignerParameterBase" ) //
        .def( "mirror", &AlignerParameterBase::mirror ) //
        .def_readonly( "name", &AlignerParameterBase::sName ) //
        .def_readonly( "description", &AlignerParameterBase::sDescription );

    exportAlignerParameter<int>( xOrganizer, "AlignerParameterInt" );
    exportAlignerParameter<short>( xOrganizer, "AlignerParameterShort" );
    exportAlignerParameter<bool>( xOrganizer, "AlignerParameterBool" );
    exportAlignerParameter<double>( xOrganizer, "AlignerParameterDouble" );
    exportAlignerParameter<float>( xOrganizer, "AlignerParameterFloat" );
    py::class_<AlignerParameter<AlignerParameterBase::ChoicesType>,
               AlignerParameterBase,
               std::shared_ptr<AlignerParameter<AlignerParameterBase::ChoicesType>>>( xOrganizer.util( ),
                                                                                      "AlignerParameterChoice" ) //
        .def( "set", &AlignerParameter<AlignerParameterBase::ChoicesType>::set ) //
        .def( "get", &AlignerParameter<AlignerParameterBase::ChoicesType>::get_py );
    // exportAlignerParameter<fs::path>( rxPyModuleId, "AlignerParameterFilePath" ); @todo

    // Export ParameterSetBase Class
    py::class_<ParameterSetBase, std::shared_ptr<ParameterSetBase>>( xOrganizer.util( ), "ParameterSetBase" ) //
        .def( "mirror", &ParameterSetBase::mirror ) //
        .def( "register_parameter", &ParameterSetBase::registerParameter ) //
        .def( "unregister_parameter", &ParameterSetBase::unregisterParameter ) //
        .def( "by_name", &ParameterSetBase::byName ) //
        .def( "by_short", &ParameterSetBase::byShort );


    // Export Presetting Class
    py::class_<Presetting, ParameterSetBase, std::shared_ptr<Presetting>>( xOrganizer.util( ), "Presetting" ) //
        .def( py::init<std::string>( ) ) //
        .def_readonly( "name", &Presetting::sName );

    // Export Presetting Class
    py::class_<GeneralParameter, ParameterSetBase, std::shared_ptr<GeneralParameter>>( xOrganizer.util( ),
                                                                                       "GeneralSettings" ) //
        .def( py::init<>( ) );

    // Export ParameterSetManager Class
    py::class_<ParameterSetManager>( xOrganizer.util( ), "ParameterSetManager" ) //
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


PYBIND11_MODULE( libMS, libMsModule )
{
    SubmoduleOrganizer xOrganizer( libMsModule );
    DEBUG_3( std::cout.setf( std::ios::unitbuf ); )
    exportParameter( xOrganizer );
    exportContainer( xOrganizer );
    exportModuleClass( xOrganizer );
#ifdef WITH_DB
    exportPoolContainer( xOrganizer );
#endif
} // function

#endif
