/* Uses the __cpuid intrinsic to get information about current platform.
 * CPU extended instruction set support.
 * For processors of types: x86, x64
 */

#pragma once

#include <array>
#include <bitset>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#if defined( __GNUG__ ) || defined( __clang__ )
#include "cpuid.h"
#include <string.h>
#elif defined( _MSC_VER )
#include <intrin.h>
#endif

class CPU_Info
{
  private:
    class CPU_Info_Internal
    {
      public:
        /* attributes */
        std::string sVendor; // processor vendor
        std::string sBrand; // processor brand
        bool isIntel_;
        bool isAMD_;
        std::bitset<32> f_1_ECX_; // processor register
        std::bitset<32> f_1_EDX_; // processor register
        std::bitset<32> f_7_EBX_; // processor register
        std::bitset<32> f_7_ECX_; // processor register
        std::bitset<32> f_81_ECX_; // processor register
        std::bitset<32> f_81_EDX_; // processor register

        CPU_Info_Internal( )
            : isIntel_{false},
              isAMD_{false},
              f_1_ECX_{0},
              f_1_EDX_{0},
              f_7_EBX_{0},
              f_7_ECX_{0},
              f_81_ECX_{0},
              f_81_EDX_{0}
        {
            std::vector<std::array<int, 4>> vData = {};
            std::vector<std::array<int, 4>> vExtData = {};
            std::array<int, 4> aCpuInfo;

            /* Calling get_cpuid with 0x0 as the function_id argument
             * gets the number of the highest valid function ID.
             */
            get_cpuid( aCpuInfo.data( ), 0 );
            auto nIds_ = aCpuInfo[ 0 ];

            for( int i = 0; i <= nIds_; ++i )
            {
                get_cpuidex( aCpuInfo.data( ), i, 0 );
                vData.push_back( aCpuInfo );
            } // for

            /* Capture vendor string */
            static_assert( sizeof( int ) == 4 ); // implies 0x20 % sizeof( int ) == 0
            int vendor[ 0x20 / sizeof( int ) ];
            memset( vendor, 0, sizeof( vendor ) );
            *reinterpret_cast<int*>( vendor ) = vData[ 0 ][ 1 ];
            *reinterpret_cast<int*>( vendor + 1 ) = vData[ 0 ][ 3 ];
            *reinterpret_cast<int*>( vendor + 2 ) = vData[ 0 ][ 2 ];
            sVendor = reinterpret_cast<char*>( vendor );
            if( sVendor == "GenuineIntel" )
            {
                isIntel_ = true;
            } // if
            else if( sVendor == "AuthenticAMD" )
            {
                isAMD_ = true;
            } // else if

            /* load bitset with flags for function 0x00000001 */
            if( nIds_ >= 1 )
            {
                f_1_ECX_ = vData[ 1 ][ 2 ];
                f_1_EDX_ = vData[ 1 ][ 3 ];
            }

            /* load bitset with flags for function 0x00000007 */
            if( nIds_ >= 7 )
            {
                f_7_EBX_ = vData[ 7 ][ 1 ];
                f_7_ECX_ = vData[ 7 ][ 2 ];
            }

            /* Calling get_cpuid with 0x80000000 as the function_id argument
             * gets the number of the highest valid extended ID.
             */
            get_cpuid( aCpuInfo.data( ), 0x80000000 );
            uint32_t nExIds_ = aCpuInfo[ 0 ];

            char aBrand[ 0x40 ];
            memset( aBrand, 0, sizeof( aBrand ) );

            for( uint32_t i = 0x80000000; i <= nExIds_; ++i )
            {
                get_cpuidex( aCpuInfo.data( ), i, 0 ); // read 16 byte form CPU reg
                vExtData.push_back( aCpuInfo );
            } // for

            /* load bitset with flags for function 0x80000001   */
            if( nExIds_ >= 0x80000001 )
            {
                f_81_ECX_ = vExtData[ 1 ][ 2 ];
                f_81_EDX_ = vExtData[ 1 ][ 3 ];
            } // if

            /* Interpret CPU brand string if reported */
            if( nExIds_ >= 0x80000004 )
            {
                memcpy( aBrand, vExtData[ 2 ].data( ), sizeof( aCpuInfo ) );
                memcpy( aBrand + 16, vExtData[ 3 ].data( ), sizeof( aCpuInfo ) );
                memcpy( aBrand + 32, vExtData[ 4 ].data( ), sizeof( aCpuInfo ) );
                sBrand = aBrand;
            } // if
        }; // constructor

        inline void get_cpuid( int cpuInfo[ 4 ], int function_id )
        {
#if defined( __GNUG__ ) || defined( __clang__ )
            __cpuid( function_id, cpuInfo[ 0 ], cpuInfo[ 1 ], cpuInfo[ 2 ], cpuInfo[ 3 ] );
#elif defined( _MSC_VER )
            __cpuid( cpuInfo, function_id );
#endif
        } // method

        inline void get_cpuidex( int cpuInfo[ 4 ], int function_id, int subfunction_id )
        {
#if defined( __GNUG__ ) || defined( __clang__ )
            __cpuid_count( function_id, subfunction_id, cpuInfo[ 0 ], cpuInfo[ 1 ], cpuInfo[ 2 ], cpuInfo[ 3 ] );
#elif defined( _MSC_VER )
            __cpuidex( cpuInfo, function_id, subfunction_id );
#endif
        } // method
    }; // class InstructionSet_Internal (inner class)

    static const CPU_Info_Internal xCPU_Rep;

  public:
    /* Getters */
    static bool SSE3( void )
    {
        return xCPU_Rep.f_1_ECX_[ 0 ];
    }
    static bool PCLMULQDQ( void )
    {
        return xCPU_Rep.f_1_ECX_[ 1 ];
    }
    static bool MONITOR( void )
    {
        return xCPU_Rep.f_1_ECX_[ 3 ];
    }
    static bool SSSE3( void )
    {
        return xCPU_Rep.f_1_ECX_[ 9 ];
    }
    static bool FMA( void )
    {
        return xCPU_Rep.f_1_ECX_[ 12 ];
    }
    static bool CMPXCHG16B( void )
    {
        return xCPU_Rep.f_1_ECX_[ 13 ];
    }
    static bool SSE41( void )
    {
        return xCPU_Rep.f_1_ECX_[ 19 ];
    }
    static bool SSE42( void )
    {
        return xCPU_Rep.f_1_ECX_[ 20 ];
    }
    static bool MOVBE( void )
    {
        return xCPU_Rep.f_1_ECX_[ 22 ];
    }
    static bool POPCNT( void )
    {
        return xCPU_Rep.f_1_ECX_[ 23 ];
    }
    static bool AES( void )
    {
        return xCPU_Rep.f_1_ECX_[ 25 ];
    }
    static bool XSAVE( void )
    {
        return xCPU_Rep.f_1_ECX_[ 26 ];
    }
    static bool OSXSAVE( void )
    {
        return xCPU_Rep.f_1_ECX_[ 27 ];
    }
    static bool AVX( void )
    {
        return xCPU_Rep.f_1_ECX_[ 28 ];
    }
    static bool F16C( void )
    {
        return xCPU_Rep.f_1_ECX_[ 29 ];
    }
    static bool RDRAND( void )
    {
        return xCPU_Rep.f_1_ECX_[ 30 ];
    }

    static bool MSR( void )
    {
        return xCPU_Rep.f_1_EDX_[ 5 ];
    }
    static bool CX8( void )
    {
        return xCPU_Rep.f_1_EDX_[ 8 ];
    }
    static bool SEP( void )
    {
        return xCPU_Rep.f_1_EDX_[ 11 ];
    }
    static bool CMOV( void )
    {
        return xCPU_Rep.f_1_EDX_[ 15 ];
    }
    static bool CLFSH( void )
    {
        return xCPU_Rep.f_1_EDX_[ 19 ];
    }
    static bool MMX( void )
    {
        return xCPU_Rep.f_1_EDX_[ 23 ];
    }
    static bool FXSR( void )
    {
        return xCPU_Rep.f_1_EDX_[ 24 ];
    }
    static bool SSE( void )
    {
        return xCPU_Rep.f_1_EDX_[ 25 ];
    }
    static bool SSE2( void )
    {
        return xCPU_Rep.f_1_EDX_[ 26 ];
    }

    static bool FSGSBASE( void )
    {
        return xCPU_Rep.f_7_EBX_[ 0 ];
    }
    static bool BMI1( void )
    {
        return xCPU_Rep.f_7_EBX_[ 3 ];
    }
    static bool HLE( void )
    {
        return xCPU_Rep.isIntel_ && xCPU_Rep.f_7_EBX_[ 4 ];
    }
    static bool AVX2( void )
    {
        return xCPU_Rep.f_7_EBX_[ 5 ];
    }
    static bool BMI2( void )
    {
        return xCPU_Rep.f_7_EBX_[ 8 ];
    }
    static bool ERMS( void )
    {
        return xCPU_Rep.f_7_EBX_[ 9 ];
    }
    static bool INVPCID( void )
    {
        return xCPU_Rep.f_7_EBX_[ 10 ];
    }
    static bool RTM( void )
    {
        return xCPU_Rep.isIntel_ && xCPU_Rep.f_7_EBX_[ 11 ];
    }
    static bool AVX512F( void )
    {
        return xCPU_Rep.f_7_EBX_[ 16 ];
    }
    static bool RDSEED( void )
    {
        return xCPU_Rep.f_7_EBX_[ 18 ];
    }
    static bool ADX( void )
    {
        return xCPU_Rep.f_7_EBX_[ 19 ];
    }
    static bool AVX512PF( void )
    {
        return xCPU_Rep.f_7_EBX_[ 26 ];
    }
    static bool AVX512ER( void )
    {
        return xCPU_Rep.f_7_EBX_[ 27 ];
    }
    static bool AVX512CD( void )
    {
        return xCPU_Rep.f_7_EBX_[ 28 ];
    }
    static bool SHA( void )
    {
        return xCPU_Rep.f_7_EBX_[ 29 ];
    }

    static bool PREFETCHWT1( void )
    {
        return xCPU_Rep.f_7_ECX_[ 0 ];
    }

    static bool LAHF( void )
    {
        return xCPU_Rep.f_81_ECX_[ 0 ];
    }
    static bool LZCNT( void )
    {
        return xCPU_Rep.isIntel_ && xCPU_Rep.f_81_ECX_[ 5 ];
    }
    static bool ABM( void )
    {
        return xCPU_Rep.isAMD_ && xCPU_Rep.f_81_ECX_[ 5 ];
    }
    static bool SSE4a( void )
    {
        return xCPU_Rep.isAMD_ && xCPU_Rep.f_81_ECX_[ 6 ];
    }
    static bool XOP( void )
    {
        return xCPU_Rep.isAMD_ && xCPU_Rep.f_81_ECX_[ 11 ];
    }
    static bool TBM( void )
    {
        return xCPU_Rep.isAMD_ && xCPU_Rep.f_81_ECX_[ 21 ];
    }

    static bool SYSCALL( void )
    {
        return xCPU_Rep.isIntel_ && xCPU_Rep.f_81_EDX_[ 11 ];
    }
    static bool MMXEXT( void )
    {
        return xCPU_Rep.isAMD_ && xCPU_Rep.f_81_EDX_[ 22 ];
    }
    static bool RDTSCP( void )
    {
        return xCPU_Rep.isIntel_ && xCPU_Rep.f_81_EDX_[ 27 ];
    }
    static bool _3DNOWEXT( void )
    {
        return xCPU_Rep.isAMD_ && xCPU_Rep.f_81_EDX_[ 30 ];
    }
    static bool _3DNOW( void )
    {
        return xCPU_Rep.isAMD_ && xCPU_Rep.f_81_EDX_[ 31 ];
    }

    static const std::string& sVendor; // export vendor
    static const std::string& sBrand; // export brand

    /* Dump a summary about the characteristics of the current CPU */
    static void showCpuInfo( std::ostream& outstream = std::cout )
    {
        auto support_message = [&outstream]( std::string isa_feature, bool is_supported ) {
            outstream << isa_feature << ( is_supported ? " supported" : " not supported" ) << std::endl;
        };

        std::cout << CPU_Info::sVendor << std::endl;
        std::cout << CPU_Info::sBrand << std::endl;

        support_message( "3DNOW", CPU_Info::_3DNOW( ) );
        support_message( "3DNOWEXT", CPU_Info::_3DNOWEXT( ) );
        support_message( "ABM", CPU_Info::ABM( ) );
        support_message( "ADX", CPU_Info::ADX( ) );
        support_message( "AES", CPU_Info::AES( ) );
        support_message( "AVX", CPU_Info::AVX( ) );
        support_message( "AVX2", CPU_Info::AVX2( ) );
        support_message( "AVX512CD", CPU_Info::AVX512CD( ) );
        support_message( "AVX512ER", CPU_Info::AVX512ER( ) );
        support_message( "AVX512F", CPU_Info::AVX512F( ) );
        support_message( "AVX512PF", CPU_Info::AVX512PF( ) );
        support_message( "BMI1", CPU_Info::BMI1( ) );
        support_message( "BMI2", CPU_Info::BMI2( ) );
        support_message( "CLFSH", CPU_Info::CLFSH( ) );
        support_message( "CMPXCHG16B", CPU_Info::CMPXCHG16B( ) );
        support_message( "CX8", CPU_Info::CX8( ) );
        support_message( "ERMS", CPU_Info::ERMS( ) );
        support_message( "F16C", CPU_Info::F16C( ) );
        support_message( "FMA", CPU_Info::FMA( ) );
        support_message( "FSGSBASE", CPU_Info::FSGSBASE( ) );
        support_message( "FXSR", CPU_Info::FXSR( ) );
        support_message( "HLE", CPU_Info::HLE( ) );
        support_message( "INVPCID", CPU_Info::INVPCID( ) );
        support_message( "LAHF", CPU_Info::LAHF( ) );
        support_message( "LZCNT", CPU_Info::LZCNT( ) );
        support_message( "MMX", CPU_Info::MMX( ) );
        support_message( "MMXEXT", CPU_Info::MMXEXT( ) );
        support_message( "MONITOR", CPU_Info::MONITOR( ) );
        support_message( "MOVBE", CPU_Info::MOVBE( ) );
        support_message( "MSR", CPU_Info::MSR( ) );
        support_message( "OSXSAVE", CPU_Info::OSXSAVE( ) );
        support_message( "PCLMULQDQ", CPU_Info::PCLMULQDQ( ) );
        support_message( "POPCNT", CPU_Info::POPCNT( ) );
        support_message( "PREFETCHWT1", CPU_Info::PREFETCHWT1( ) );
        support_message( "RDRAND", CPU_Info::RDRAND( ) );
        support_message( "RDSEED", CPU_Info::RDSEED( ) );
        support_message( "RDTSCP", CPU_Info::RDTSCP( ) );
        support_message( "RTM", CPU_Info::RTM( ) );
        support_message( "SEP", CPU_Info::SEP( ) );
        support_message( "SHA", CPU_Info::SHA( ) );
        support_message( "SSE", CPU_Info::SSE( ) );
        support_message( "SSE2", CPU_Info::SSE2( ) );
        support_message( "SSE3", CPU_Info::SSE3( ) );
        support_message( "SSE4.1", CPU_Info::SSE41( ) );
        support_message( "SSE4.2", CPU_Info::SSE42( ) );
        support_message( "SSE4a", CPU_Info::SSE4a( ) );
        support_message( "SSSE3", CPU_Info::SSSE3( ) );
        support_message( "SYSCALL", CPU_Info::SYSCALL( ) );
        support_message( "TBM", CPU_Info::TBM( ) );
        support_message( "XOP", CPU_Info::XOP( ) );
        support_message( "XSAVE", CPU_Info::XSAVE( ) );
    } // method
}; // class InstructionSet (outer class)

#define DP_CPU_ENFORCE_SSE 0

/* CPU depended dispatcher
 */
template <typename TP_RETURN, typename TP_FUN, typename... TP_ARGS>
TP_RETURN dispatchbyCPU( std::bitset<1> xFlags, TP_FUN&& fun_AVX2, TP_FUN&& fun_SSE41, TP_ARGS&&... args )
{
    if( CPU_Info::AVX2( ) && !xFlags.test( DP_CPU_ENFORCE_SSE ) )
    {
        //// std::cout << "Dispatch to AVX2 code " << std::endl;
        return std::bind( std::forward<TP_FUN>( fun_AVX2 ), std::forward<TP_ARGS>( args )... )( );
    } // if
    else if( CPU_Info::SSE41( ) )
    {
        //// std::cout << "Dispatch to SSE41 code " << std::endl;
        return std::bind( std::forward<TP_FUN>( fun_SSE41 ), std::forward<TP_ARGS>( args )... )( );
    } // else if
    else
    {
        std::cout << "CPU:" << CPU_Info::sVendor << " " << CPU_Info::sBrand << " is not supported by this build."
                  << std::endl;
        exit( 0 );
    } // else
} // function