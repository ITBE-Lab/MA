/* kswcpp for environments that support AVX2. 
 * File must be compiled with appropriate flags:
 *	MSVC: /arch:AVX2
 *  GCC: -mavx2
 */
#if defined( __GNUC__ ) || defined( _MSC_VER )
	#ifndef  __AVX2__
		// #error "Wrong architecture flag during compilation."
	#endif
#else
	#error "Compiler unsupported."
#endif

#include "kswcpp.h" // general kswcpp header
#include "avx2.h" // include must be before kswcpp_core.h include
//#include "scalar.h" // include must be before kswcpp_core.h include // DEBUG
#include "kswcpp_core.h" // uses AVX2 instructions now

/* AVX2 variant of kswcpp.
 * Ensure that your CPU supports the appropriate instruction set!
 */
void kswcpp_avx2(
	int qlen,				// query length -- TO DO: size_t
	const uint8_t *query,	// query sequence
	int tlen,				// reference length -- TO DO: size_t
	const uint8_t *target,	// reference sequence
#ifdef USE_CPP_PARAM
	const KswCppParam<5> &xParam,
#else
	int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2, // scoring
#endif
	int w,					// bandwidth
	int zdrop,				// zdrop value
	int flag,				// flags for DP control
	kswcpp_extz_t *ez,		// outcome of alignment
	AlignedMemoryManager &rxMemManager )
{
	/* Dispatch according to possible max and min scores */
	if( !xParam.riskOfOverflow<int16_t>( std::max( qlen, tlen ) ) )
	{	/* int16_t scoring */
		kswcpp_core<
			int16_t, // type used for scoring        
			Vec_m256i_16_int16_t, // SIMD vector-type used for scoring (128 bit)
			int8_t, // Signed type for differences
			uint8_t, // Unsigned type for differences
            // Vec_scalar<64, int8_t> >
			Vec_m256i_32_int8_t> // SIMD-vector type for differences (128 bit); AVX 2: Vec_m256i_32_int8_t
			( qlen, query, tlen, target,
#ifdef USE_CPP_PARAM
			  xParam,
#else
			  m, mat, q, e, q2, e2,
#endif  
			  w, zdrop, flag,
			  ez,
			  rxMemManager );
	} // if
	else if( !xParam.riskOfOverflow<int32_t>( std::max( qlen, tlen ) ) )
	{	/* int32_t scoring */
		kswcpp_core<
			int32_t, // type used for scoring
			Vec_m256i_8_int32_t, // SIMD vector-type used for scoring (128 bit)
			int8_t, // Signed type for differences
			uint8_t, // Unsigned type for differences
			Vec_m256i_32_int8_t> // SIMD-vector type for differences (128 bit); AVX 2: Vec_m256i_32_int8_t
			( qlen, query, tlen, target,
#ifdef USE_CPP_PARAM
			  xParam,
#else
			  m, mat, q, e, q2, e2,
#endif  
			  w, zdrop, flag,
			  ez,
			  rxMemManager );
	} // else if
	else
	{}// TO DO: Throw exception that indicates an overflow.
} // function