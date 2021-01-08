// InstructionSet.cpp
// Compile by using: cl /EHsc /W4 InstructionSet.cpp
// processor: x86, x64
// Uses the __cpuid intrinsic to get information about
// CPU extended instruction set support.

#include "cpu_info.h"

/* Initialize static member data */
const CPU_Info::CPU_Info_Internal CPU_Info::xCPU_Rep;
std::string sFilePrefix = "/MAdata/tmp/.CIGARMemoryManager";
size_t uiKswHashTableGbMinSize = 8;

/* Export vendor and brand as const strings */
const std::string& CPU_Info::sVendor = xCPU_Rep.sVendor;
const std::string& CPU_Info::sBrand = xCPU_Rep.sBrand;

// Print out supported instruction set extensions
void show_cpu_info( )
{
    CPU_Info::showCpuInfo( );
} // function
