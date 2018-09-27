#pragma once
#ifdef __INTELLISENSE__
#include "cuda_runtime_api.h"
#endif
/* The below configuration parameter effect the kernel merely */

/* Use the optimized indexing scheme */
#define OPTIMIZED_INDEXING ( 1 )

/* Use registers instead of shared memory as cache.
 * Boost runtime tremendously. (However, maximum STRIPE_WIDTH is smaller than with shared memory)
 */
#define USE_REGISTER_CACHE ( 1 )

/* The kernel copies the scoring profile first to the shared memory. */
#define SCORING_PROFILE_IN_SHARED_MEMORY ( 1 )

/* Scores in constant memory is rather disadvantageous */
#define SCORES_IN_CONSTANT_MEMORY ( 0 )

/* Use the checksum optimization */
#define USE_CHECKSUM ( 1 )

/* Use the checksum optimization */
#define USE_MAXZERO ( 0 )

template <typename TP> __device__ TP gpuMax( TP first, TP second )
{
    return max( first, second );
} // device function

template <> __device__ float gpuMax( float first, float second )
{
    /* __device__ ​ float fmaxf ( float  x, float  y )
     *  Determine the maximum numeric value of the arguments.
     */
    return fmaxf( first, second );
} // device function


/* Non-branching max(x, 0)
 * Instead of branching we have a comparison and multiplication.
 */
template <typename TP> __device__ TP maxZero( TP iValue )
{
    return ( iValue > static_cast<TP>( 0 ) ) * iValue;
} // device function

#if 1
/* Stripe Kernel That works for packed shorts. (two short values in one int value)
 * SEGMENT size should be a power of 2.
 * WARNING: int is a 32-bit value on GPU side, so max is 2^32 NT references.
 * IDEA: Reference into the texture memory for fast reading?
 * SIMD on CUDA:
 * http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__INTRINSIC__SIMD.html#group__CUDA__MATH__INTRINSIC__SIMD
 * Optimization: Keep iH_ValueLeft - (GAP_INIT + GAP_EXT) in the pHFM_Vector, this avoids on
 * arithmetic operation. TO DO:
 */
template <unsigned int TILE_HEIGHT, // only necessary in the case USE_REGISTER_CACHE == 0
          unsigned int STRIPE_WIDTH, // width of the stripe
          bool LAZY_MODE> // In lazy mode the kernel stops at equal checksums
__global__ void
tileKernel16(
    uint2 *pvHF_VectorIn, // Vector has size of reference sequence H1, H2, F1, F2
    uint2 *pvHF_VectorOut, // Vector has size of reference sequence
    unsigned int *pvM_VectorIn, // Vector has size of reference sequence M1(short), M2(short),
    int *pvC_Vector, // Checksum considers left side and right side together
    uint2 *pvHE_CacheIn, // H1, H2, E1, E2 size: (gridDim.x * blockDim.x) * STRIPE_WIDTH equal to
                         // uiNumberOfSegments * STRIPE_WIDTH
    uint2 *pvHE_CacheOut, // H1, H2, E1, E2 size: like pvHE_CacheInt
    char *pRefAnchor, // Address of reference on device. Organized as pair of 4 bit.
#if( SCORES_IN_CONSTANT_MEMORY != 1 )
    const int4 *i4_Scores, // Scoring Profile (should be moved to constant memory) Optimization: 16
                           // bit scores should be enough!
#endif
    char *pContinuationFlagVector, // Only used in the lazy-mode to communicate the need of a
                                   // further kernel call.
    const unsigned int uiSegmentSize // equal to TILE_HEIGHT
)
{ /* !! NUM_OF_SEGMENTS must be equal to gridDim.x * blockDim.x !! */
    const unsigned int &tid = threadIdx.x;
    const unsigned int &bid = blockIdx.x;

    unsigned int uiSegmentId = ( bid * blockDim.x + tid );
    unsigned int uiNumberOfSegments = gridDim.x * blockDim.x;
    unsigned int uiTotalSize = gridDim.x * blockDim.x * uiSegmentSize;

    /* Gap opening and gap extension as two short values packed into one unsigned int */
    const unsigned int iv2GapOpenScore = ( ( GAP_INIT + GAP_EXT ) << 16 ) |
                                         ( GAP_INIT + GAP_EXT ); // init-gap 16-bit, init-gap 16 bit
    const unsigned int iv2GapExtScore =
        ( GAP_EXT << 16 ) | GAP_EXT; // gap-ext 16-bit, gap-ext 16 bit

#if( USE_CHECKSUM == 1 )
    /* The Flush-Mode can be only active in LAZY-MODE.
     * Possible optimization: IN LAZY_MODE first check the re-computation flag of the block before
     * and depending on this flag set the Flush-Mode over here.
     */
    bool bFlushMode = false;
#endif

#if( SCORING_PROFILE_IN_SHARED_MEMORY == 1 )
    /* WARNING: uiSegmentSize MUST be >= STRIPE_WIDTH or you are in big trouble! */
    __shared__ int4 aProfileCache[ STRIPE_WIDTH ];

    if( tid < STRIPE_WIDTH )
    {
        aProfileCache[ tid ] = i4_Scores[ tid ];
    } // if
#endif
    /* IMPORTANT: All threads have to wait until all score data are fully loaded */
    __syncthreads( );

    /* Each thread of the block owns one row in all matrices.
     * blockDim.x is equal to TILE_HEIGHT.
     * OBSERVATION: Long vectors can result in bad runtimes.
     */
#if( USE_REGISTER_CACHE == 1 )
    uint2 HE_Cache[ STRIPE_WIDTH ];
#else
    __shared__ uint2 HE_Cache[ TILE_HEIGHT * STRIPE_WIDTH ];
#endif
    for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
    { /* The threads write parallel to consecutive addresses.
       * This avoids banking conflicts.
       */
        if( !LAZY_MODE )
        {
#if( USE_REGISTER_CACHE == 1 )
            HE_Cache[ uiColIndex ].x = 0; // H1 = 0 and H2 = 0
            HE_Cache[ uiColIndex ].y = 0; // E1 = 0 and E2 = 0
#else
            int uiIndex = blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_HEIGHT
            HE_Cache[ uiIndex ].x = 0;
            HE_Cache[ uiIndex ].y = 0;
#endif
        } // if
        else
        { /* In the lazy-mode we have to reload the HE-cache*/
            unsigned int uiBasicIndex = ( uiNumberOfSegments * uiColIndex );
            if( uiSegmentId == 0 )
            { /* Segment 0 special, because there is no previous segment */
                unsigned int uiGlobalMemIndex = uiBasicIndex + ( uiNumberOfSegments - 1 );
                HE_Cache[ uiColIndex ].x = pvHE_CacheIn[ uiGlobalMemIndex ].x >> 16;
                HE_Cache[ uiColIndex ].y = pvHE_CacheIn[ uiGlobalMemIndex ].y >> 16;
            } // if
            else
            { /* Get the cached H-values and E-values of the previous block(segment).
               * These values were computed as outcome of the previous call of the kernel.
               */
                unsigned int uiGlobalMemIndex = uiBasicIndex + ( uiSegmentId - 1 );
#if( USE_REGISTER_CACHE == 1 )
                HE_Cache[ uiColIndex ] =
                    pvHE_CacheIn[ uiGlobalMemIndex ]; // get H1, H2, E1, E2 from previous block
#else
                int uiCacheIndex =
                    blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_HEIGHT
                HE_Cache[ uiCacheIndex ] = pvHE_CacheIn[ uiGlobalMemIndex ];
#endif
            } // (uiSegmentId != 0)
        } // else (in LAZY_MODE)
    } // for (column iteration for cache initialization)

    /* In the beginning of each iteration it keeps up-left H-value of the first column.
     * Threads read consecutive addresses in parallel (but misaligned)
     */
    unsigned int iH_ValueLeftUp; // H1-LeftUp(Upper 16 bit), H2-LeftUp(Lower 16 bit)
    if( uiSegmentId > 0 )
    {
        iH_ValueLeftUp = pvHF_VectorIn[ uiTotalSize - uiNumberOfSegments + ( uiSegmentId - 1 ) ].x;
    } // if
    else
    { /* The very first segment is special here and has to be initialized always by 0 and H from the
         final row*/
        iH_ValueLeftUp = pvHF_VectorIn[ uiTotalSize - 1 ].x >> 16;
    } // id

    /* The kernel now iterates over the rows of the segment.*/
#if( OPTIMIZED_INDEXING == 1 )
    unsigned int uiIndex = uiSegmentId;
    unsigned int uiRowCounter;
    for( uiRowCounter = 0; uiRowCounter < uiSegmentSize; ++uiRowCounter )
#else
    unsigned int uiSegmentStart = uiSegmentSize * uiSegmentId;
    for( unsigned int uiIndex = uiSegmentStart; uiIndex < uiSegmentStart + uiSegmentSize;
         ++uiIndex )
    //// for (unsigned int uiRow = uiSegmentStart; uiRow < uiSegmentStart + uiSegmentSize; ++uiRow)
#endif
    {
#if( USE_CHECKSUM == 1 )
        if( LAZY_MODE )
        { /* In flush-mode it is not necessary to compute rows anymore.
           * We know that all remaining rows of the current segment are already correctly computed.
           */
            if( bFlushMode )
            {
                continue;
            } // if
        } // if
#endif
        /* The 'left-values' for the current row are in the ingoing HF-Vector */
        const uint2 i2_HF_Cell =
            pvHF_VectorIn[ uiIndex ]; // H1(short), H2(short), E1(short), E2(short)
        unsigned int iH_ValueLeft = i2_HF_Cell.x; // H1-left(short), H2-left(short)
        unsigned int iF_ValueLeft = i2_HF_Cell.y; // F1-left(short), F2-left(short)

        /* Read the overall maximum score for the current row so far */
        unsigned int iH_ValueMax =
            pvM_VectorIn[ uiIndex ]; // M1(short), M2(short) previously i2_MC_CellIn.x;

        /* Get the checksum(s) for the current row*/
        //// const int2 i2_C_CellIn = pvC_Vector[uiIndex]; // C1 (.x), C2 (.y)

        /* Coalesced Read. So very efficient.
         * The symbol comprises two reference positions simultaneously.
         */
        char iSymbol = pRefAnchor[ uiIndex ]; // (c1 (higher nibble), c2 (lower nibble))
        char iRightSymbol = iSymbol & 0x0F; // extract the lower nibble of iSymbol
        char iLeftSymbol = iSymbol >> 4; // extract the higher nibble of iSymbol

        unsigned int iH_Value;
        unsigned int iE_Value;
        unsigned int iF_Value;

#if( USE_CHECKSUM == 1 )
        /* The iHE_Checksum keeps the sum of all H-values and E-values of the current row.
         * HE-checksum can be negative. So, we need a signed int type.
         * OPTIMIZATION: Use a modulo checksum on vector side.
         */
        int iHE_Checksum;
        iHE_Checksum = 0;
#endif

        /* Inner loop that is iterating over the columns.
         * Maximum size of a stripe is 2^16.
         * This loop can use CUDA-SIMD instructions for optimization.
         */
        for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
        {
#if( USE_REGISTER_CACHE == 1 )
            unsigned int iH_ValueUp = HE_Cache[ uiColIndex ].x; // H1-Up, H2-Up
            unsigned int iE_ValueUp = HE_Cache[ uiColIndex ].y; // E1-Up, E2-Up
#else
            /* We organize the matrix in column layout for avoiding banking conflicts */
            int uiCellIndex = ( blockDim.x * uiColIndex ) + tid;

            /* Take the existing values in the cache as H-value and E-value of the previous row */
            int &iH_ValueUp = HE_Cache[ uiCellIndex ].x;
            int &iE_ValueUp = HE_Cache[ uiCellIndex ].y;
#endif
            /* Compute provisional E-value and H-value of the current cell.
             * This values can be later improved by the lazy-e loop (kernel).
             * H-value keeps the value of the previous iteration.
             * Implements: iE_Value = gpuMax( iH_ValueUp - (GAP_INIT + GAP_EXT), iE_ValueUp -
             * GAP_EXT );
             * __vsub2 performs per-halfword (un)signed substraction, with wrap-around.
             * __vmaxs2 performs per-halfword signed maximum computation.
             */
            iE_Value = __vmaxs2( __vsub2( iH_ValueUp, iv2GapOpenScore ),
                                 __vsub2( iE_ValueUp, iv2GapExtScore ) );

            /* The below statement results in the catch of a single int from global memory.
             * iSymbol comprises two reference symbols simultaneously (see PPT for more info).
             * We get S1, S2.
             */
#if( SCORING_PROFILE_IN_SHARED_MEMORY == 1 )
            /* The scores are still stored as 32-bit values.
             * We do two lookups in the scoring profile (one for left symbol, one for right symbol)
             * and combine the two scores in one 32 bit integer.
             */
            const unsigned int iv2_Score = ( ( reinterpret_cast<const unsigned int *>(
                                               &aProfileCache[ uiColIndex ] ) )[ iRightSymbol ] ) &
                                               0x0000FFFF |
                                           ( ( reinterpret_cast<const unsigned int *>(
                                                 &aProfileCache[ uiColIndex ] ) )[ iLeftSymbol ]
                                             << 16 );
#else
            int iScore = ( reinterpret_cast<const int *>( &i4_Scores[ uiColIndex ] ) )[ iSymbol ];
#endif
            /* iH_ValueLeft keeps currently up-left H-value.
             * Implements: iH_Value = gpuMax( iH_ValueLeftUp + iScore, iE_Value );
             * __vadd2 performs per-halfword (un)signed addition, with wrap-around: a + b.
             */
            iH_Value = __vmaxs2( __vadd2( iH_ValueLeftUp, iv2_Score ), iE_Value );

            /* Update the F-value for the current cell.
             * Implements: iF_Value = gpuMax( iH_ValueLeft - (GAP_INIT + GAP_EXT), iF_ValueLeft -
             * GAP_EXT );
             */
            iF_Value = __vmaxs2( __vsub2( iH_ValueLeft, iv2GapOpenScore ),
                                 __vsub2( iF_ValueLeft, iv2GapExtScore ) );

            /* 1. Check whether the F-value delivers a better H-value.
             * 2. H-value must be greater than 0
             * Implementation 1: iH_Value = gpuMax( iH_Value, iF_Value );
             * Implementation 2: iH_Value = gpuMax( iH_Value, 0 );
             */
            iH_Value = __vmaxs2( iH_Value, iF_Value );
            iH_Value = __vmaxs2( iH_Value, 0 );

            /* Update the maximum for the current row.
             * Implementation: iH_ValueMax = gpuMax( iH_Value, iH_ValueMax );
             */
            iH_ValueMax = __vmaxs2( iH_Value, iH_ValueMax );

            /* Prepare next iteration */
            iH_ValueLeftUp =
                iH_ValueUp; // the H-value-up will be the H-value-left-up in the next iteration
            iH_ValueLeft =
                iH_Value; // the current H-value will be the H-value left in the next iteration
            iF_ValueLeft =
                iF_Value; // the current E-value will be the E-value left in the next iteration

#if( USE_CHECKSUM == 1 )
            /* Checksum Update (Remember Checksums must stay 32 bit values, in order to avoid
             * overflows) Implementation: iHE_Checksum += iH_Value + iE_Value; Idea: One 32-bit
             * checksum for left side and right side simultaneously. Idea: Implement the checksum by
             * using one large 64-bit long value
             */
            iHE_Checksum += iH_Value >> 16; // left side H-value
            iHE_Checksum += iE_Value >> 16; // left side E-value

            iHE_Checksum += iH_Value & 0x0000FFFF; // right side H-value
            iHE_Checksum += iE_Value & 0x0000FFFF; // right side E-value
#endif

#if( USE_REGISTER_CACHE == 1 )
            HE_Cache[ uiColIndex ].x = iH_Value; // Keep the H-Value in the cache
            HE_Cache[ uiColIndex ].y = iE_Value; // Keep the E-Value in the cache
#else
            HE_Cache[ uiCellIndex ].x = iH_Value; // log the provisional H-Value in the cache
            HE_Cache[ uiCellIndex ].y = iE_Value; // log the provisional E-Value in the cache
#endif
        } // for column

        /* In the next row iteration the H-value-left-up is the current H-left-value (of the first
         * column). Setting up-left via the left H-value safes one access to global memory.
         */
        iH_ValueLeftUp = i2_HF_Cell.x;

        uint2 i2_HF_CellOut;
        /* Prepare the int 4 value for logging the values of the last column into the pHFM_Vector.
         */
        i2_HF_CellOut.x = iH_Value; // H-value of last column (written to HF-output vector)
        i2_HF_CellOut.y = iF_Value; // F-value of last column (written to HF-output vector)

        /* Update the cell in the vector with the values of the current column */
        pvHF_VectorOut[ uiIndex ] = i2_HF_CellOut;

#if( USE_CHECKSUM == 1 )
        /* Possible small optimization: Write one int2 over here */
        if( LAZY_MODE )
        { /* Flag value 0 : Stop, Flag value 1: continue.
           * If the checksum is equal for this row, than all H-values and E-values stayed unaltered
           * and we can simply flush the remaining H-values and F-values.
           * 32-Bit code: bFlushMode = (iHE_Checksum == i2_MC_CellIn.y);
           */
            int &i2_CachedChecksum = pvC_Vector[ uiIndex ];
            bFlushMode = iHE_Checksum == i2_CachedChecksum;
        } // if

        /* Update the checksum */
        pvC_Vector[ uiIndex ] = iHE_Checksum;
#endif
        /* Update the row maximum*/
        pvM_VectorIn[ uiIndex ] = iH_ValueMax; // Write Max 1, Max 2

#if( OPTIMIZED_INDEXING == 1 )
        uiIndex += gridDim.x * blockDim.x; // increment the index by segment-size
#endif // OPTIMIZED_INDEXING
    } // for (row counting)

    /* Lead out ... */
#if( USE_CHECKSUM == 1 )
    /* Evaluate out current working mode */
    if( LAZY_MODE && bFlushMode )
    { /* Copy the existing in-cache in global memory to the out-cache in global memory */
        for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
        { /* The threads write parallel to consecutive addresses.
           * This avoids banking conflicts.
           */
            unsigned int uiGlobalMemIndex = ( uiNumberOfSegments * uiColIndex ) + uiSegmentId;
            pvHE_CacheOut[ uiGlobalMemIndex ] = pvHE_CacheIn[ uiGlobalMemIndex ];
        } // for
    } // if
    else
    { /* If we did not set the flush-mode, the computation recognized differences continuing to the
       * last row. In this case we have to write the thread-cache to the global out-cache.
       */
        for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
        { /* The threads write parallel to consecutive addresses.
           * This avoids banking conflicts.
           */
            unsigned int uiGlobalMemIndex = ( uiNumberOfSegments * uiColIndex ) + uiSegmentId;
#if( USE_REGISTER_CACHE == 1 )
            pvHE_CacheOut[ uiGlobalMemIndex ] = HE_Cache[ uiColIndex ];
#else
            int uiCacheIndex =
                blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_PER_BLOCK
            pvHE_CacheOut[ uiGlobalMemIndex ] = HE_Cache[ uiCacheIndex ];
#endif
        } // for
    } // else (bFlushMode)
    if( LAZY_MODE ) // static compile switch
    { /* Set the re-computation flag in the continuation vector */
        pContinuationFlagVector[ uiSegmentId ] = !bFlushMode;
    } // if (LAZY_MODE)
#else // (NOT) USE_CHECKSUM
    for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
    { /* The threads write parallel to consecutive addresses.
       * This avoids banking conflicts.
       */
        unsigned int uiGlobalMemIndex = ( uiNumberOfSegments * uiColIndex ) + uiSegmentId;
        pvHE_CacheOut[ uiGlobalMemIndex ].x = HE_Cache[ uiColIndex ].x;
        pvHE_CacheOut[ uiGlobalMemIndex ].y = HE_Cache[ uiColIndex ].y;
    } // for

    // Check for termination in lazy mode
    if( LAZY_MODE )
    {
        pContinuationFlagVector[ uiSegmentId ] = 0;

        for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
        { /* The threads write parallel to consecutive addresses.
           * This avoids banking conflicts.
           */
            int uiGlobalMemIndex = ( uiNumberOfSegments * uiColIndex ) + uiSegmentId;
#if( USE_REGISTER_CACHE == 1 )
            /* Important: We must check-value and E-value and H-value */
            char bE_differend =
                ( HE_Cache[ uiColIndex ].x != pvHE_CacheIn[ uiGlobalMemIndex ].x ) ||
                ( HE_Cache[ uiColIndex ].y != pvHE_CacheIn[ uiGlobalMemIndex ].y );
#else
            unsigned int uiCacheIndex =
                blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_HEIGHT
            char bE_differend =
                ( HE_Cache[ uiCacheIndex ].x != pvHE_CacheIn[ uiGlobalMemIndex ].x ) ||
                ( HE_Cache[ uiCacheIndex ].y != pvHE_CacheIn[ uiGlobalMemIndex ].y );
#endif
            if( bE_differend )
            { /* We found at least one different E. So, we have to compute once again. */
                pContinuationFlagVector[ uiSegmentId ] = 1;
                break;
            } // if
        } // for
    } // if LAZY_MODE
#endif // USE_CHECKSUM
} // CUDA kernel
#endif // PROPOSED 16-bit SCORING KERNEL
#if 0
/* TILE Kernel with flexible datatypes.
 * SEGMENT size should be a power of 2.
 * WARNING: int is a 32-bit value on GPU side, so max is 2^32 NT references.
 * IDEA: Reference into the texture memory for fast reading?
 * SIMD on CUDA: http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__INTRINSIC__SIMD.html#group__CUDA__MATH__INTRINSIC__SIMD
 * Optimization: Keep iH_ValueLeft - (GAP_INIT + GAP_EXT) in the pHFM_Vector, this avoids on arithmetic operation.
 * TO DO: 
 */
template<
	typename T_PAIR, // type of a pair of scores (int2 for int scores, short2 for short scores, char4 for char scores)
	typename T_SCORE, // type of a single score value (must be chosen in accordance with T_PAIR)
	unsigned int TILE_HEIGHT, // only necessary in the case USE_REGISTER_CACHE == 0
	unsigned int STRIPE_WIDTH,  // width of the stripe
	bool LAZY_MODE>
__global__
void tileKernel( T_PAIR* pvHF_VectorIn, // Vector has size of reference sequence
				 T_PAIR* pvHF_VectorOut, // Vector has size of reference sequence
				 int2* pvM_VectorIn, // Vector has size of reference sequence (must first stay int2 because of the checksum)
				 T_PAIR* pvHE_CacheIn, // size: (gridDim.x * blockDim.x) * STRIPE_WIDTH equal to uiNumberOfSegments * STRIPE_WIDTH
				 T_PAIR* pvHE_CacheOut, // size: like pvHE_CacheInt
				 char* pRefAnchor, // Address of reference on device (GPU side)
#if( SCORES_IN_CONSTANT_MEMORY != 1 )
				 const int4* i4_Scores, // Scoring Profile (should be moved to constant memory) Optimization: 16 bit scores should be enough!
#endif
				 char* pContinuationFlagVector, // Only used in the lazy-mode to communicate the need of a further kernel call.
				 const unsigned int uiSegmentSize // equal to TILE_HEIGHT
)
{	/* !! NUM_OF_SEGMENTS must be equal to gridDim.x * blockDim.x !! */
	const unsigned int &tid = threadIdx.x;
	const unsigned int &bid = blockIdx.x;

	unsigned int uiSegmentId = (bid * blockDim.x + tid);
	unsigned int uiNumberOfSegments = gridDim.x * blockDim.x;
	unsigned int uiTotalSize = gridDim.x * blockDim.x * uiSegmentSize;

#if( USE_CHECKSUM == 1 )
	/* The Flush-Mode can be only active in LAZY-MODE.
	* Possible optimization: IN LAZY_MODE first check the re-computation flag of the block before and 
	* depending on this flag set the Flush-Mode over here.
	*/
	bool bFlushMode = false;
#endif

#if( SCORING_PROFILE_IN_SHARED_MEMORY == 1 )
	/* WARNING: uiSegmentSize MUST be >= STRIPE_WIDTH or you are in big trouble! */
	__shared__ int4 aProfileCache[STRIPE_WIDTH];

	if( tid < STRIPE_WIDTH )
	{
		aProfileCache[tid] = i4_Scores[tid];
	} // if
#endif
	  /* IMPORTANT: All threads have to wait until all score data are fully loaded */
	__syncthreads();

	/* Each thread of the block owns one row in all matrices.
	* blockDim.x is equal to TILE_HEIGHT.
	* OBSERVATION: Long vectors can result in bad runtimes.
	*/
#if( USE_REGISTER_CACHE == 1 )
	T_PAIR HE_Cache[STRIPE_WIDTH];
#else
	__shared__ T_PAIR HE_Cache[TILE_HEIGHT * STRIPE_WIDTH];
#endif
	for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
	{	/* The threads write parallel to consecutive addresses.
		 * This avoids banking conflicts.
		 */
		if( !LAZY_MODE || ( uiSegmentId == 0 ) )
		{	/* Clear the cache in the first call (non-lazy).
			 * The first segment is special and has to be always initialized with 0.
			 */
#if( USE_REGISTER_CACHE == 1 )
			HE_Cache[uiColIndex].x = 0;
			HE_Cache[uiColIndex].y = 0; // -(GAP_INIT + GAP_EXT); //  0;
#else
			int uiIndex = blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_HEIGHT
			HE_Cache[uiIndex].x = 0;
			HE_Cache[uiIndex].y = -(GAP_INIT + GAP_EXT); //  0;
#endif
		} // if
		else
		{	/* Get the cached H-values and E-values of the previous block(segment).
			 * These values were computed as outcome of the previous call of the kernel.
			 */
			unsigned int uiGlobalMemIndex = (uiNumberOfSegments * uiColIndex) + (uiSegmentId - 1);
#if( USE_REGISTER_CACHE == 1 )
			HE_Cache[uiColIndex] = pvHE_CacheIn[uiGlobalMemIndex];
#else
			int uiCacheIndex = blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_HEIGHT
			HE_Cache[uiCacheIndex] = pvHE_CacheIn[uiGlobalMemIndex];
#endif
		} // else 
	} // for (column iteration for cache initialization)

	  /* In the beginning of each iteration it keeps up-left H-value of the first column.
	  * Threads read consecutive addresses in parallel (but misaligned)
	  */
	T_SCORE iH_ValueLeftUp = uiSegmentId > 0 ? pvHF_VectorIn[uiTotalSize - uiNumberOfSegments + (uiSegmentId - 1)].x
											 : 0; /* The very first segment is special here and has to be initialized always by 0 */

			 /* The kernel now iterates over the rows of the segment.*/
#if( OPTIMIZED_INDEXING == 1 )
	unsigned int uiIndex = uiSegmentId;
	unsigned int uiRowCounter;
	for( uiRowCounter = 0; uiRowCounter < uiSegmentSize; ++uiRowCounter )
#else
	unsigned int uiSegmentStart = uiSegmentSize * uiSegmentId;
	for( unsigned int uiIndex = uiSegmentStart; uiIndex < uiSegmentStart + uiSegmentSize; ++uiIndex )
		//// for (unsigned int uiRow = uiSegmentStart; uiRow < uiSegmentStart + uiSegmentSize; ++uiRow)
#endif
	{
#if( USE_CHECKSUM == 1 )
		if( LAZY_MODE )
		{	/* In flush-mode it is not necessary to compute rows anymore.
			* We know that all remaining rows of the current segment are already correctly computed.
			*/
			if( bFlushMode )
			{ 
				continue;
			} // if
		} // if
#endif
		  /* The 'left-values' for the current row are in the ingoing HF-Vector */
		const T_PAIR i2_HF_Cell = pvHF_VectorIn[uiIndex];
		const int2 i2_MC_CellIn = pvM_VectorIn[uiIndex];

		T_SCORE iH_ValueLeft = i2_HF_Cell.x; // First component of 4-tuple is H-value
		T_SCORE iF_ValueLeft = i2_HF_Cell.y; // Second component of 4-tuple is F-value

										 /* If we write i2_ME_Cell.x, we are in trouble */
		/* Read the overall maximum score for the current row so far */
		T_SCORE iH_ValueMax = i2_MC_CellIn.x;

		/*  Coalesced Read. So very efficient */
		char iSymbol = pRefAnchor[uiIndex];

		T_SCORE iH_Value;
		T_SCORE iE_Value;
		T_SCORE iF_Value;

#if( USE_CHECKSUM == 1 )
		/* The iHE_Checksum keeps the sum of all H-values and E-values of the current row.
		* HE-checksum can be negative. So, we need a signed int type.
		*/
		int iHE_Checksum = 0;
#endif

		/* Inner loop that is iterating over the columns.
		 * Maximum size of a stripe is 2^16.
		 * This loop can use CUDA-SIMD instructions for optimization. 
		 */
		for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
		{
#if( USE_REGISTER_CACHE == 1 )
			T_SCORE iH_ValueUp = HE_Cache[uiColIndex].x;
			T_SCORE iE_ValueUp = HE_Cache[uiColIndex].y;
#else
			/* We organize the matrix in column layout for avoiding banking conflicts */
			int uiCellIndex = (blockDim.x * uiColIndex) + tid;

			/* Take the existing values in the cache as H-value and E-value of the previous row */
			int &iH_ValueUp = HE_Cache[uiCellIndex].x;
			int &iE_ValueUp = HE_Cache[uiCellIndex].y;
#endif
			/* Compute provisional E-value and H-value of the current cell.
			* This values can be later improved by the lazy-e loop (kernel).
			* H-value keeps the value of the previous iteration.
			*/
			iE_Value = gpuMax( iH_ValueUp - (GAP_INIT + GAP_EXT), iE_ValueUp - GAP_EXT );

			/* The below statement results in the catch of a single int from global memory. */
#if( SCORING_PROFILE_IN_SHARED_MEMORY == 1 )
			T_SCORE iScore = static_cast<T_SCORE>(((reinterpret_cast<const int *>(&aProfileCache[uiColIndex]))[iSymbol]));
#else
			int iScore = (reinterpret_cast<const int *>(&i4_Scores[uiColIndex]))[iSymbol];
#endif
			/* iH_ValueLeft keeps currently up-left H-value */
			iH_Value = gpuMax( iH_ValueLeftUp + iScore, iE_Value );

			/* Update the F-value for the current cell */
			iF_Value = gpuMax( iF_ValueLeft - GAP_EXT_HORIZONAL, iH_ValueLeft - (GAP_INIT + GAP_EXT_HORIZONAL) );

			/* 1. Check whether the F-value delivers a better H-value.
			* 2. H-value must be greater than 0
			*/
			iH_Value = gpuMax( iH_Value, iF_Value );
			iH_Value = gpuMax( iH_Value, 0 ); // Use arithmetic over here

											  /* Update the maximum for the current row */
			iH_ValueMax = gpuMax( iH_Value, iH_ValueMax );

			/* Prepare next iteration */
			iH_ValueLeftUp = iH_ValueUp; // the H-value-up will be the H-value-left-up in the next iteration
			iH_ValueLeft = iH_Value; // the current H-value will be the H-value left in the next iteration
			iF_ValueLeft = iF_Value; // the current E-value will be the E-value left in the next iteration

#if( USE_CHECKSUM == 1 )
									 /* Checksum Update */
			iHE_Checksum += iH_Value + iE_Value;
#endif

#if( USE_REGISTER_CACHE == 1 )
			HE_Cache[uiColIndex].x = iH_Value; // Keep the H-Value in the cache 
			HE_Cache[uiColIndex].y = iE_Value; // Keep the E-Value in the cache
#else
			HE_Cache[uiCellIndex].x = iH_Value; // log the provisional H-Value in the cache 
			HE_Cache[uiCellIndex].y = iE_Value; // log the provisional E-Value in the cache
#endif
		} // for column

		  /* In the next row iteration the H-value-left-up is the current H-left-value (of the first column).
		  * Setting up-left via the left H-value safes one access to global memory.
		  */
		iH_ValueLeftUp = i2_HF_Cell.x;

		T_PAIR i2_HF_CellOut;
		/* Prepare the int 4 value for logging the values of the last column into the pHFM_Vector.
		*/
		i2_HF_CellOut.x = iH_Value; // i4_Cell.x = iH_Value; // H-value of last column
		i2_HF_CellOut.y = iF_Value; // i4_Cell.y = iF_Value; // F-value of last column

									/* Update the cell in the vector with the values of the current column */
		pvHF_VectorOut[uiIndex] = i2_HF_CellOut;

#if( USE_CHECKSUM == 1 )
		/* Possible small optimization: Write one int2 over here */
		if( LAZY_MODE )
		{	/* Flag value 0 : Stop, Flag value 1: continue.
			* If the checksum is equal for this row, than all H-values and E-values stayed unaltered
			* and we can simply flush the remaining H-values and F-values.
			*/
			bFlushMode = (iHE_Checksum == i2_MC_CellIn.y);
		} // if

		  /* Update the checksum */
		pvM_VectorIn[uiIndex].y = iHE_Checksum;
#endif
		/* Update the row maximum*/
		pvM_VectorIn[uiIndex].x = iH_ValueMax;

#if( OPTIMIZED_INDEXING == 1 )
		uiIndex += gridDim.x * blockDim.x; // increment the index by segment-size
#endif // OPTIMIZED_INDEXING
	} // for (row counting)

	  /* Lead out ... */
#if( USE_CHECKSUM == 1 )
	  /* Evaluate out current working mode */
	if( LAZY_MODE && bFlushMode )
	{	/* Copy the existing in-cache in global memory to the out-cache in global memory */
		for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
		{  /* The threads write parallel to consecutive addresses.
		   * This avoids banking conflicts.
		   */
			unsigned int uiGlobalMemIndex = (uiNumberOfSegments * uiColIndex) + uiSegmentId;
			pvHE_CacheOut[uiGlobalMemIndex] = pvHE_CacheIn[uiGlobalMemIndex];
		} // for
	} // if
	else
	{	/* If we did not set the flush-mode, the computation recognized differences continuing to the last row.
		* In this case we have to write the thread-cache to the global out-cache.
		*/
		for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
		{  /* The threads write parallel to consecutive addresses.
		   * This avoids banking conflicts.
		   */
			unsigned int uiGlobalMemIndex = (uiNumberOfSegments * uiColIndex) + uiSegmentId;
#if( USE_REGISTER_CACHE == 1 )
			pvHE_CacheOut[uiGlobalMemIndex] = HE_Cache[uiColIndex];
#else
			int uiCacheIndex = blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_PER_BLOCK
			pvHE_CacheOut[uiGlobalMemIndex] = HE_Cache[uiCacheIndex];
#endif
		} // for
	} // else (bFlushMode)
	if( LAZY_MODE ) // static compile switch
	{	/* Set the re-computation flag in the continuation vector */
		pContinuationFlagVector[uiSegmentId] = !bFlushMode;
	} // if (LAZY_MODE)
#else // (NOT) USE_CHECKSUM 
	for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
	{  /* The threads write parallel to consecutive addresses.
	   * This avoids banking conflicts.
	   */
		unsigned int uiGlobalMemIndex = (uiNumberOfSegments * uiColIndex) + uiSegmentId;
		pvHE_CacheOut[uiGlobalMemIndex].x = HE_Cache[uiColIndex].x;
		pvHE_CacheOut[uiGlobalMemIndex].y = HE_Cache[uiColIndex].y;
	} // for

	  // Check for termination in lazy mode
	if( LAZY_MODE )
	{	
		pContinuationFlagVector[uiSegmentId] = 0;

		for( unsigned short uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
		{  /* The threads write parallel to consecutive addresses.
		   * This avoids banking conflicts.
		   */
			int uiGlobalMemIndex = (uiNumberOfSegments * uiColIndex) + uiSegmentId;
#if( USE_REGISTER_CACHE == 1 )
			/* Important: We must check-value and E-value and H-value */
			char bE_differend =	   (HE_Cache[uiColIndex].x != pvHE_CacheIn[uiGlobalMemIndex].x)
				|| (HE_Cache[uiColIndex].y != pvHE_CacheIn[uiGlobalMemIndex].y);
#else
			unsigned int uiCacheIndex = blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_HEIGHT
			char bE_differend =	   (HE_Cache[uiCacheIndex].x != pvHE_CacheIn[uiGlobalMemIndex].x)
				|| (HE_Cache[uiCacheIndex].y != pvHE_CacheIn[uiGlobalMemIndex].y);
#endif		
			if( bE_differend )
			{	/* We found at least one different E. So, we have to compute once again. */
				pContinuationFlagVector[uiSegmentId] = 1;
				break;
			} // if
		} // for
	} // if LAZY_MODE
#endif // USE_CHECKSUM	
} // CUDA kernel
#endif // WORGING KERNEL WITH FLEXIBLE TYPES
#if 1
/* The HE-cache has to be cleared before the first kernel-call.
 * SEGMENT size should be a power of 2.
 * WARNING: int is a 32-bit value on GPU side, so max is 2^32 NT references.
 * IDEA: Reference into the texture memory for fast reading?
 * SIMD on CUDA:
 * http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__INTRINSIC__SIMD.html#group__CUDA__MATH__INTRINSIC__SIMD
 * Optimization: Keep iH_ValueLeft - (GAP_INIT + GAP_EXT) in the pHFM_Vector, this avoids on
 * arithmetic operation.
 */
template <typename SCORE_TP4, typename SCORE_TP2, typename SCORE_TP, typename CHECKSUM_TP,
          unsigned int TILE_PER_BLOCK, unsigned int STRIPE_WIDTH, bool LAZY_MODE,
          bool FULL_STRIPE_WIDTH> // ignore parameter uiWidth, use full stripe width instead.
__global__ void
tileKernel( SCORE_TP2 *pvHF_VectorIn, // Vector has size of reference sequence
            SCORE_TP2 *pvHF_VectorOut, // Vector has size of reference sequence
            SCORE_TP *pvM_VectorIn, // Vector has size of reference sequence
            CHECKSUM_TP *pvC_Vector, // Checksum vector
            const SCORE_TP2 *pvHE_CacheInSegment0, // size: STRIPE_WIDTH
            const SCORE_TP xLeftUpH_Value, // left up H-value for reconstruction
            SCORE_TP2 *pvHE_CacheIn, // size: (gridDim.x * blockDim.x) * STRIPE_WIDTH equal to
                                     // uiNumberOfSegments * STRIPE_WIDTH
            SCORE_TP2 *pvHE_CacheOut, // size: like pvHE_CacheInt
            char *pRefAnchor, // Address of reference on device (GPU side)
#if( SCORES_IN_CONSTANT_MEMORY != 1 )
            const SCORE_TP4 *i4_Scores, // Scoring Profile (should be moved to constant memory)
#endif
            char *pContinuationFlagVector, // Only used in the lazy-mode to communicate the need of
                                           // a further kernel call.
            const unsigned int
                uiSegmentSize, // uiNumberOfSegments * uiSegmentSize is equal to the total size
            const unsigned int uiWidth // number of columns computed with current stripe (
                                       // assert(uiWidth <= STRIPE_WIDTH) )
)
{ /* !! NUM_OF_SEGMENTS must be equal to gridDim.x * blockDim.x !! */
    const unsigned int &tid = threadIdx.x;
    const unsigned int &bid = blockIdx.x;

    unsigned int uiSegmentId = ( bid * blockDim.x + tid );
    unsigned int uiNumberOfSegments = gridDim.x * blockDim.x;
    unsigned int uiTotalSize =
        gridDim.x * blockDim.x * uiSegmentSize; // uiNumberOfSegments * uiSegmentSize

#if( USE_CHECKSUM == 1 )
    /* The Flush-Mode can be only active in LAZY-MODE.
     * Possible optimization: IN LAZY_MODE first check the re-computation flag of the block before
     * and depending on this flag set the Flush-Mode over here.
     */
    bool bFlushMode = false;
#endif

#if( SCORING_PROFILE_IN_SHARED_MEMORY == 1 )
    /* WARNING: uiSegmentSize MUST be >= STRIPE_WIDTH or you are in big trouble! */
    __shared__ SCORE_TP4 aProfileCache[ STRIPE_WIDTH ];

    if( tid < STRIPE_WIDTH )
    {
        aProfileCache[ tid ] = i4_Scores[ tid ];
    } // if
#endif
    /* IMPORTANT: All threads have to wait until all score data are fully loaded */
    __syncthreads( );

    /* Each thread of the block owns one row in all matrices.
     * blockDim.x is equal to TILE_PER_BLOCK.
     * OBSERVATION: Long vectors can result in bad runtimes.
     */
#if( USE_REGISTER_CACHE == 1 )
    SCORE_TP2 HE_Cache[ STRIPE_WIDTH ];
#else
    __shared__ SCORE_TP2 HE_Cache[ TILE_PER_BLOCK * STRIPE_WIDTH ];
#endif
    /* Initialize the HE_Cache appropriately */
    for( unsigned int uiColIndex = 0; uiColIndex < ( FULL_STRIPE_WIDTH ? STRIPE_WIDTH : uiWidth );
         ++uiColIndex )
    { /* The threads write parallel to consecutive addresses.
       * This avoids banking conflicts.
       */
        //// if( !LAZY_MODE || ( uiSegmentId == 0 ) )
        if( uiSegmentId == 0 )
        { /* The first segment is currently special and gets it values form a separated vector.
           * (This solution is not optimal but does it for the moment.)
           */
#if( USE_REGISTER_CACHE == 1 )
            //// HE_Cache[uiColIndex].x = 0;
            //// HE_Cache[uiColIndex].y = 0;
            HE_Cache[ uiColIndex ] = pvHE_CacheInSegment0[ uiColIndex ];
#else
            int uiIndex = blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_PER_BLOCK
            HE_Cache[ uiIndex ].x = 0;
            HE_Cache[ uiIndex ].y = -( GAP_INIT + GAP_EXT ); //  0;
#endif
        } // if
        else
        { /* Get the cached H-values and E-values of the previous block(segment).
           * These values were computed as outcome of the previous call of the kernel.
           */
            unsigned int uiGlobalMemIndex =
                ( uiNumberOfSegments * uiColIndex ) + ( uiSegmentId - 1 );
#if( USE_REGISTER_CACHE == 1 )
            HE_Cache[ uiColIndex ] = pvHE_CacheIn[ uiGlobalMemIndex ];
#else
            int uiCacheIndex =
                blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_PER_BLOCK
            HE_Cache[ uiCacheIndex ] = pvHE_CacheIn[ uiGlobalMemIndex ];
#endif
        } // else
    } // for (column iteration for cache initialization)

    /* In the beginning of each iteration it keeps up-left H-value of the first column.
     * Threads read consecutive addresses in parallel (but misaligned)
     */
    SCORE_TP iH_ValueLeftUp =
        uiSegmentId > 0 ? pvHF_VectorIn[ uiTotalSize - uiNumberOfSegments + ( uiSegmentId - 1 ) ].x
                        : xLeftUpH_Value; //  0; /* The very first segment is special here and has
                                          //  to be initialized always by 0 */

    /* The kernel now iterates over the rows of the segment.*/
#if( OPTIMIZED_INDEXING == 1 )
    unsigned int uiIndex = uiSegmentId;
    unsigned int uiRowCounter;
    for( uiRowCounter = 0; uiRowCounter < uiSegmentSize; ++uiRowCounter )
#else
    unsigned int uiSegmentStart = uiSegmentSize * uiSegmentId;
    for( unsigned int uiIndex = uiSegmentStart; uiIndex < uiSegmentStart + uiSegmentSize;
         ++uiIndex )
    //// for (unsigned int uiRow = uiSegmentStart; uiRow < uiSegmentStart + uiSegmentSize; ++uiRow)
#endif
    {
#if( USE_CHECKSUM == 1 )
        if( LAZY_MODE )
        { /* In flush-mode it is not necessary to compute rows anymore.
           * We know that all remaining rows of the current segment are already correctly computed.
           */
            if( bFlushMode )
            {
                continue;
            } // if
        } // if
#endif
        /* The 'left-values' for the current row are in the ingoing HF-Vector */
        SCORE_TP2 i2_HF_Cell = pvHF_VectorIn[ uiIndex ];

        ////? SCORE_TP iH_ValueLeft = i2_HF_Cell.x; // First component of 4-tuple is H-value
        ////? SCORE_TP iF_ValueLeft = i2_HF_Cell.y; // Second component of 4-tuple is F-value

        const SCORE_TP iH_ValueLeftBackup = i2_HF_Cell.x; // added Jan. 27

        /* If we write i2_ME_Cell.x, we are in trouble */
        SCORE_TP &iH_ValueMax = pvM_VectorIn[ uiIndex ]; // added reference Jan. 27

        /* Coalesced Read. So very efficient */
        char iSymbol = pRefAnchor[ uiIndex ];

        SCORE_TP &iH_Value = i2_HF_Cell.x; // added reference Jan. 27
        //// SCORE_TP iE_Value;
        SCORE_TP &iF_Value = i2_HF_Cell.y; // added reference Jan. 27

#if( USE_CHECKSUM == 1 )
        /* The iHE_Checksum keeps the sum of all H-values and E-values of the current row.
         * HE-checksum can be negative. So, we need a signed int type.
         * The checksum should stay 32 bit.
         */
        CHECKSUM_TP iHE_Checksum = 0;
#endif

        /* Inner loop that is iterating over the columns.
         * A tremendous performance boost comes from unrolling this loop by the compiler.
         */
        for( unsigned int uiColIndex = 0;
             uiColIndex < ( FULL_STRIPE_WIDTH ? STRIPE_WIDTH : uiWidth ); ++uiColIndex )
        {
#if( USE_REGISTER_CACHE == 1 )
            SCORE_TP &iH_ValueUp = HE_Cache[ uiColIndex ].x; // added reference Jan. 27
            SCORE_TP &iE_ValueUp = HE_Cache[ uiColIndex ].y; // added reference Jan. 27
#else
            /* We organize the matrix in column layout for avoiding banking conflicts */
            int uiCellIndex = ( blockDim.x * uiColIndex ) + tid;

            /* Take the existing values in the cache as H-value and E-value of the previous row */
            int &iH_ValueUp = HE_Cache[ uiCellIndex ].x;
            int &iE_ValueUp = HE_Cache[ uiCellIndex ].y;
#endif
            /* Compute provisional E-value and H-value of the current cell.
             * This values can be later improved by the lazy-e loop (kernel).
             * H-value keeps the value of the previous iteration.
             */
            SCORE_TP iE_Value = gpuMax( iH_ValueUp - ( GAP_INIT + GAP_EXT ), iE_ValueUp - GAP_EXT );

            /* The below statement results in the catch of a single int from global memory. */
#if( SCORING_PROFILE_IN_SHARED_MEMORY == 1 )
            SCORE_TP iScore =
                ( reinterpret_cast<const SCORE_TP *>( &aProfileCache[ uiColIndex ] ) )[ iSymbol ];
#else
            SCORE_TP iScore =
                ( reinterpret_cast<const SCORE_TP *>( &i4_Scores[ uiColIndex ] ) )[ iSymbol ];
#endif
            /* Update the F-value for the current cell.
             * At this moment iH_Value still contains the H-value of the left cell.
             */
            iF_Value =
                gpuMax( iF_Value - GAP_EXT_HORIZONAL, iH_Value - ( GAP_INIT + GAP_EXT_HORIZONAL ) );
            //// before iF_Value = gpuMax( iF_ValueLeft - GAP_EXT_HORIZONAL, iH_ValueLeft -
            ///(GAP_INIT + GAP_EXT_HORIZONAL) );

            /* iH_ValueLeft keeps currently up-left H-value.
             * By this assignment we destroy the knowledge about the left H-value. (This is not
             * needed any more)
             */
            iH_Value = gpuMax( iH_ValueLeftUp + iScore, iE_Value );

            /* 1. Check whether the F-value delivers a better H-value.
             * 2. H-value must be greater than 0
             */
            iH_Value = gpuMax( iH_Value, iF_Value );

#if( USE_MAXZERO == 1 )
            iH_Value = maxZero( iH_Value ); // Use arithmetic over here
#else
            iH_Value = gpuMax( iH_Value, static_cast<SCORE_TP>( 0 ) ); // Use arithmetic over here
#endif // USE_MAXZERO

            /* Update the maximum for the current row */
            iH_ValueMax = gpuMax( iH_Value, iH_ValueMax );

            /* Prepare next iteration */
            iH_ValueLeftUp =
                iH_ValueUp; // the H-value-up will be the H-value-left-up in the next iteration
            //// ?? iH_ValueLeft = iH_Value; // the current H-value will be the H-value left in the
            /// next iteration / ?? iF_ValueLeft = iF_Value; // the current E-value will be the
            /// E-value left in the next iteration

#if( USE_CHECKSUM == 1 )
            /* Checksum Update.
             * The checksum data-type must have been chosen large enough for avoiding overflows.
             */
            iHE_Checksum += static_cast<CHECKSUM_TP>( iH_Value + iE_Value );
#endif

#if( USE_REGISTER_CACHE == 1 )
            HE_Cache[ uiColIndex ].x = iH_Value; // Keep the H-Value in the cache
            HE_Cache[ uiColIndex ].y = iE_Value; // Keep the E-Value in the cache
#else
            HE_Cache[ uiCellIndex ].x = iH_Value; // log the provisional H-Value in the cache
            HE_Cache[ uiCellIndex ].y = iE_Value; // log the provisional E-Value in the cache
#endif
        } // for column

        /* In the next row iteration the H-value-left-up is the current H-left-value (of the first
         * column). Setting up-left via the left H-value safes one access to global memory.
         */
        iH_ValueLeftUp = iH_ValueLeftBackup; // Jan. 27
        //// Jan. 27 iH_ValueLeftUp = i2_HF_Cell.x;

        //// Jan. 27 SCORE_TP2 i2_HF_CellOut;
        //// Jan. 27 /* Prepare the int 4 value for logging the values of the last column into the
        /// pHFM_Vector. / Jan. 27  */ / Jan. 27 i2_HF_CellOut.x = iH_Value; // i4_Cell.x =
        /// iH_Value;
        ///// H-value of last column / Jan. 27 i2_HF_CellOut.y = iF_Value; // i4_Cell.y = iF_Value;
        ///// F-value of last column / Jan. 27 / Jan. 27 /* Update the cell in the vector with the
        /// values of the current column */ / Jan. 27 pvHF_VectorOut[uiIndex] = i2_HF_CellOut;

        pvHF_VectorOut[ uiIndex ] = i2_HF_Cell;

#if( USE_CHECKSUM == 1 )
        /* Possible small optimization: Write one int2 over here */
        if( LAZY_MODE )
        { /* Flag value 0 : Stop, Flag value 1: continue.
           * If the checksum is equal for this row, than all H-values and E-values stayed unaltered
           * and we can simply flush the remaining H-values and F-values.
           */
            //// **** bFlushMode = (iHE_Checksum == i2_MC_CellIn.y);
            bFlushMode = ( iHE_Checksum == pvC_Vector[ uiIndex ] );
        } // if

        /* Update the checksum */
        //// **** pvM_VectorIn[uiIndex].y = iHE_Checksum;
        pvC_Vector[ uiIndex ] = iHE_Checksum;
#endif
        /* Update the row maximum*/
        //// **** pvM_VectorIn[uiIndex].x = iH_ValueMax;
        pvM_VectorIn[ uiIndex ] = iH_ValueMax;

#if( OPTIMIZED_INDEXING == 1 )
        uiIndex += gridDim.x * blockDim.x; // increment the index by segment-size
#endif // OPTIMIZED_INDEXING
    } // for (row counting)

    /* Lead out ... */
#if( USE_CHECKSUM == 1 )
    /* Evaluate out current working mode */
    if( LAZY_MODE && bFlushMode )
    { /* Copy the existing in-cache in global memory to the out-cache in global memory */
        for( int uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
        { /* The threads write parallel to consecutive addresses.
           * This avoids banking conflicts.
           */
            unsigned int uiGlobalMemIndex = ( uiNumberOfSegments * uiColIndex ) + uiSegmentId;
            pvHE_CacheOut[ uiGlobalMemIndex ] = pvHE_CacheIn[ uiGlobalMemIndex ];
        } // for
    } // if
    else
    { /* If we did not set the flush-mode, the computation recognized differences continuing to the
       * last row. In this case we have to write the thread-cache to the global out-cache.
       */
        for( int uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
        { /* The threads write parallel to consecutive addresses.
           * This avoids banking conflicts.
           */
            unsigned int uiGlobalMemIndex = ( uiNumberOfSegments * uiColIndex ) + uiSegmentId;
#if( USE_REGISTER_CACHE == 1 )
            pvHE_CacheOut[ uiGlobalMemIndex ] = HE_Cache[ uiColIndex ];
#else
            int uiCacheIndex =
                blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_PER_BLOCK
            pvHE_CacheOut[ uiGlobalMemIndex ] = HE_Cache[ uiCacheIndex ];
#endif
        } // for
    } // else (bFlushMode)
    if( LAZY_MODE ) // static compile switch
    { /* Set the re-computation flag in the continuation vector */
        pContinuationFlagVector[ uiSegmentId ] = !bFlushMode;
    } // if (LAZY_MODE)
#else // (NOT) USE_CHECKSUM
    for( unsigned int uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
    { /* The threads write parallel to consecutive addresses.
       * This avoids banking conflicts.
       */
        unsigned int uiGlobalMemIndex = ( uiNumberOfSegments * uiColIndex ) + uiSegmentId;
        pvHE_CacheOut[ uiGlobalMemIndex ].x = HE_Cache[ uiColIndex ].x;
        pvHE_CacheOut[ uiGlobalMemIndex ].y = HE_Cache[ uiColIndex ].y;
    } // for

    // Check for termination in lazy mode
    if( LAZY_MODE )
    {
        pContinuationFlagVector[ uiSegmentId ] = 0;

        for( int uiColIndex = 0; uiColIndex < STRIPE_WIDTH; ++uiColIndex )
        { /* The threads write parallel to consecutive addresses.
           * This avoids banking conflicts.
           */
            unsigned int uiGlobalMemIndex = ( uiNumberOfSegments * uiColIndex ) + uiSegmentId;
#if( USE_REGISTER_CACHE == 1 )
            /* Important: We must check-value and E-value and H-value */
            char bE_differend =
                ( HE_Cache[ uiColIndex ].x != pvHE_CacheIn[ uiGlobalMemIndex ].x ) ||
                ( HE_Cache[ uiColIndex ].y != pvHE_CacheIn[ uiGlobalMemIndex ].y );
#else
            unsigned int uiCacheIndex =
                blockDim.x * uiColIndex + tid; // blockDim.x is equal to TILE_PER_BLOCK
            char bE_differend =
                ( HE_Cache[ uiCacheIndex ].x != pvHE_CacheIn[ uiGlobalMemIndex ].x ) ||
                ( HE_Cache[ uiCacheIndex ].y != pvHE_CacheIn[ uiGlobalMemIndex ].y );
#endif
            if( bE_differend )
            { /* We found at least one different E. So, we have to compute once again. */
                pContinuationFlagVector[ uiSegmentId ] = 1;
                break;
            } // if
        } // for
    } // if LAZY_MODE
#endif // USE_CHECKSUM
} // CUDA kernel
#endif // WORKING KERNEL FOR 32 bit