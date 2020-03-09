cmake_minimum_required (VERSION 3.8)

project("MSV Lib" VERSION 1.1.1 DESCRIPTION "The Modular Structural Variant Calling Library")

ms_add_module( MSV )
target_link_libraries ( libMSV LINK_PUBLIC libMA )
target_compile_definitions( libMSV PUBLIC DLL_PORT_MA=dllimport )