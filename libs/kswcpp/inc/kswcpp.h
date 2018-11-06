/* TO DO:
 * CHECK THE HORIZONTAL MAX WITH AVX2! (It might be broken)
 */
#pragma once
#if defined(_MSC_VER)
	#define NOMINMAX
#endif

#include <stdint.h>
#include <vector>
#include <array>
#include <iostream>
#include <limits>
#include <algorithm>
#include "kswcpp_mem.h"
#include "cpu_info.h"

#define USE_CPP_PARAM

/* Flags for DP control */
#define KSW_EZ_SCORE_ONLY  0x01 // don't record alignment path/cigar
#define KSW_EZ_RIGHT       0x02 // right-align gaps
#define KSW_EZ_GENERIC_SC  0x04 // switch to matrix based scoring
#define KSW_EZ_APPROX_MAX  0x08 // approximate max; this is faster with sse
#define KSW_EZ_APPROX_DROP 0x10 // approximate Z-drop; faster with sse
#define KSW_EZ_EXTZ_ONLY   0x40 // only perform extension
#define KSW_EZ_REV_CIGAR   0x80 // reverse CIGAR in the output
#define KSW_EZ_FORCE_SSE   0x100 // for SSEX.X instruction use, even if AVX2 is possible


struct kswcpp_extz_t 
{
	uint32_t max:31, zdropped:1;
	int max_q, max_t;      // max extension coordinate
	int mqe, mqe_t;        // max score when reaching the end of query
	int mte, mte_q;        // max score when reaching the end of target
	int score;             // max score reaching both ends; may be KSW_NEG_INF
	int m_cigar, n_cigar;
	int reach_end;
	uint32_t *cigar;		// cigar is represented sequence of uint32_t values
}; // struct

/* Parameter for controlling kswcpp */
template<unsigned int ALPHABET_SIZE = 5>
class KswCppParam
{
public:
	/* Attributes */
	std::array<int8_t, ALPHABET_SIZE * ALPHABET_SIZE> mat; 
	const int8_t m = ALPHABET_SIZE; // alphabet size; last character has to be wildcard! (N)
	int8_t q = 4; int8_t e = 2; int8_t q2 = 24; int8_t e2 = 1; // gap scoring

	int8_t max_sc; // (ksw2) maximum match score 
	int8_t min_sc; // (ksw2) minimum mismatch score 
	int8_t iOverallMinScr; // overall minimum score
	
	/* Default constructor */
	KswCppParam() :
		max_sc( 2 ), min_sc( -4 )
	{
		fillScoringMatrix();
	} // default constructor

	/* Constructor */
	KswCppParam(
		int8_t iMatchSc, // match score
		int8_t iMisMatchSc, // mismatch score
		int8_t iGap, // affine gap, gap open costs costs
		int8_t iExtend, // affine gap, gap extension costs
		int8_t iGap2, // affine gap, gap open costs (2-piece)
		int8_t iExtend2 // affine gap, gap extension costs (2-piece) 
	) :
		q( iGap ), e( iExtend ), q2( iGap2 ), e2( iExtend2 ),
		max_sc( iMatchSc < 0 ? -iMatchSc : iMatchSc ),
		min_sc( iMisMatchSc > 0 ? -iMisMatchSc : iMisMatchSc )
	{
		fillScoringMatrix();
		this->iOverallMinScr = std::min<int>({ min_sc, -q, -e, -q2, -e2 });
	} // constructor


	void fillScoringMatrix( void )
	{
		for( auto i = 0U; i < ALPHABET_SIZE - 1; ++i ) {
			for( auto j = 0U; j < ALPHABET_SIZE - 1; ++j )
				mat[i * ALPHABET_SIZE + j] = (i == j) ? this->max_sc : this->min_sc;
			mat[i * ALPHABET_SIZE + ALPHABET_SIZE - 1] = 0;
		} // for
		for( auto j = 0U; j < ALPHABET_SIZE; ++j )
			mat[(ALPHABET_SIZE - 1) * ALPHABET_SIZE + j] = 0;
	} // method


	/* Check for the given datatype and input size for a possible overflow with respect
	 * to the scoring values.
	 */
	template<typename TP_SCR>
	bool riskOfOverflow( int64_t iSize ) const
	{
		//DBG std::cout << "overall min: " << (int)this->iOverallMinScr
		//DBG 	<< "\noverall max: " << (int)this->max_sc
		//DBG 	<< "\n type min: " << std::numeric_limits<TP_SCR>::min()
		//DBG 	<< "\n type max: " << std::numeric_limits<TP_SCR>::max()
		//DBG 	<< "\n type size: " << sizeof( TP_SCR )
		//DBG 	<< "\npossible max score: " << iSize * this->max_sc
		//DBG 	<< "\npossible min score: " << iSize * this->iOverallMinScr
		//DBG 	<< std::endl;

		int64_t i64MaxTp = static_cast<int64_t>(std::numeric_limits<TP_SCR>::max());
		int64_t i64MinTp = static_cast<int64_t>(std::numeric_limits<TP_SCR>::min());
		return (iSize * this->iOverallMinScr < i64MinTp) || (iSize * this->max_sc > i64MaxTp);
	} // method


	/* Sets the attributes max_sc and min_sc according to the given matrix*/
	void setMaxMinFromMatrix( void )
	{
		this->max_sc = mat[0]; this->min_sc = mat[1];
		for( auto iItr = 1; iItr < ALPHABET_SIZE * ALPHABET_SIZE; ++iItr )
		{
			this->max_sc = this->max_sc > mat[iItr] ? this->max_sc : mat[iItr];
			this->min_sc = this->min_sc < mat[iItr] ? this->min_sc : mat[iItr];
		} // for
	} // method
}; // class 

/* kswcpp using SSE_XX instructions */
void kswcpp_sse_xx(
	int qlen,				// query length -- TO DO: size_t
	const uint8_t *query,	// query sequence
	int tlen,				// reference length -- TO DO: size_t
	const uint8_t *target,	// reference sequence
#ifdef USE_CPP_PARAM
	const KswCppParam<5> &xPara,
#else
	int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2, // scoring
#endif
	int w,					// bandwidth
	int zdrop,				// zdrop value
	int flag,				// flags for DP control
	kswcpp_extz_t *ez,		// outcome of alignment
	AlignedMemoryManager &rxMemManager );


/* kswcpp using AVX2 instructions*/
void kswcpp_avx2(
	int qlen,				// query length -- TO DO: size_t
	const uint8_t *query,	// query sequence
	int tlen,				// reference length -- TO DO: size_t
	const uint8_t *target,	// reference sequence
#ifdef USE_CPP_PARAM
	const KswCppParam<5> &xPara,
#else
	int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2, // scoring
#endif
	int w,					// bandwidth
	int zdrop,				// zdrop value
	int flag,				// flags for DP control
	kswcpp_extz_t *ez,		// outcome of alignment
	AlignedMemoryManager &rxMemManager );

/* kswcpp main dispatcher */
static inline void kswcpp_dispatch(
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
    auto vFlags = std::bitset<1>( 1 << DP_CPU_ENFORCE_SSE );
     //<1>( 1 << DP_CPU_ENFORCE_SSE );
	dispatchbyCPU<void>
	(	vFlags,
		kswcpp_avx2,
		kswcpp_sse_xx,
		qlen, query, tlen, target,
#ifdef USE_CPP_PARAM
		xParam,
#else
		m, mat, q, e, q2, e2,
#endif  
		w, zdrop, flag,
		ez,
		rxMemManager
	); // dispatch
} // inlined static function