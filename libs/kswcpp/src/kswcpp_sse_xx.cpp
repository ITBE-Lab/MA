/* kswcpp for Intel and AMD CPUs that support SSE_XX-family instrcutions.
 * File must be compiled with appropriate switches:
 *	MSVC: For 64 bit, no flags required.
 *  GCC: Requires flag: -msse4.1
 */
#if defined( __GNUC__ )
#ifndef __SSE4_1__
#error "Wrong architecture flag."
#endif
#endif

#include "cpu_info.h" // information about current platform
#include "sse_xx.h" // include must be before kswcpp_core.h include
// DEBUG #include "scalar.h"
#include "kswcpp.h" // general kswcpp header
#include "kswcpp_core.h"

/* The kswcpp dispatcher for the SSE-family of instruction sets.
 * INSTR_SET: Family of SSE-instruction supported by the current processor.
 */
template <unsigned int INSTR_SET>
void kswcpp_sse_xx_instr_set( int qlen, // query length -- TO DO: size_t
                              const uint8_t* query, // query sequence
                              int tlen, // reference length -- TO DO: size_t
                              const uint8_t* target, // reference sequence
#ifdef USE_CPP_PARAM
                              const KswCppParam<5>& xParam,
#else
                              int8_t m, const int8_t* mat, int8_t q, int8_t e, int8_t q2, int8_t e2, // scoring
#endif
                              int w, // bandwidth
                              int zdrop, // zdrop value
                              int flag, // flags for DP control
                              kswcpp_extz_t* ez, // outcome of alignment
                              AlignedMemoryManager& rxMemManager )
{
    /* Dispatch according to possible max and min scores */
    if( !xParam.riskOfOverflow<int16_t>( std::max( qlen, tlen ) ) )
    { /* int16_t scoring */
        kswcpp_core<int16_t, // type used for scoring
                    Vec_m128i_8_int16_t<INSTR_SET>, // SIMD vector-type used for scoring (128 bit)
                    int8_t, // Signed type for differences
                    uint8_t, // Unsigned type for differences
                    // DEBUG Vec_scalar<16, int8_t>
                    Vec_m128i_16_int8_t<INSTR_SET> // SIMD-vector type for differences (128 bit)
                    >( qlen, query, tlen, target,
#ifdef USE_CPP_PARAM
                       xParam,
#else
                       m, mat, q, e, q2, e2,
#endif
                       w, zdrop, flag, ez, rxMemManager );
    } // if
    else if( !xParam.riskOfOverflow<int32_t>( std::max( qlen, tlen ) ) )
    { /* int32_t scoring */
        kswcpp_core<int32_t, // type used for scoring
                    Vec_m128i_4_int32_t<INSTR_SET>, // SIMD vector-type used for scoring (128 bit)
                    int8_t, // Signed type for differences
                    uint8_t, // Unsigned type for differences
                    Vec_m128i_16_int8_t<INSTR_SET> // SIMD-vector type for differences (128 bit)
                    >( qlen, query, tlen, target,
#ifdef USE_CPP_PARAM
                       xParam,
#else
                       m, mat, q, e, q2, e2,
#endif
                       w, zdrop, flag, ez, rxMemManager );
    } // if
    else
    {} // TO DO: Throw exception that indicates an overflow.
} // function

/* The kswcpp dispatcher for the SSE-family of instruction sets.
 * Dispatches according to the SSE-instruction set supported by the current processor.
 * (This code cannot be inlined, because it has to be compiled with specific flags.)
 */
void kswcpp_sse_xx( int qlen, // query length -- TO DO: size_t
                    const uint8_t* query, // query sequence
                    int tlen, // reference length -- TO DO: size_t
                    const uint8_t* target, // reference sequence
#ifdef USE_CPP_PARAM
                    const KswCppParam<5>& xParam,
#else
                    int8_t m, const int8_t* mat, int8_t q, int8_t e, int8_t q2, int8_t e2, // scoring
#endif
                    int w, // bandwidth
                    int zdrop, // zdrop value
                    int flag, // flags for DP control
                    kswcpp_extz_t* ez, // outcome of alignment
                    AlignedMemoryManager& rxMemManager )
{
    if( CPU_Info::SSE41( ) )
    { /* Platform supports SSE 4.1 (includes SSSE3, SSE2) */
        // std::cout << "CPU_Info::SSE41" << std::endl;
        kswcpp_sse_xx_instr_set<41>( qlen, query, tlen, target,
#ifdef USE_CPP_PARAM
                                     xParam,
#else
                                     m, mat, q, e, q2, e2,
#endif
                                     w, zdrop, flag, ez, rxMemManager );
    } // if
    else if( CPU_Info::SSSE3( ) )
    { /* Platform supports SSSE3 (includes SSE2) */
        kswcpp_sse_xx_instr_set<31>( qlen, query, tlen, target,
#ifdef USE_CPP_PARAM
                                     xParam,
#else
                                     m, mat, q, e, q2, e2,
#endif
                                     w, zdrop, flag, ez, rxMemManager );
    } // else if
    else
    { /* SSE2 only ... */
        kswcpp_sse_xx_instr_set<20>( qlen, query, tlen, target,
#ifdef USE_CPP_PARAM
                                     xParam,
#else
                                     m, mat, q, e, q2, e2,
#endif
                                     w, zdrop, flag, ez, rxMemManager );
    } // else
} // function