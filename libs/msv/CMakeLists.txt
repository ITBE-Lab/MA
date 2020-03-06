cmake_minimum_required (VERSION 3.8)

project("MSV Lib" VERSION 1.1.1 DESCRIPTION "The Modular Structural Variant Calling Library")

ms_add_module( MSV )
target_link_libraries ( MSV LINK_PUBLIC MA )
target_compile_definitions( MSV PUBLIC DLL_PORT_MA=dllimport )