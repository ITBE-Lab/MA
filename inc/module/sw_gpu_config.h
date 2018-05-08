#pragma once

#define DO_TESTS ( 0 )

#define COMPARE_SSE2 ( 0 )

#define USE_TILES ( 1 )

#define CARD ( 1080 )

#define USE_UNIFIED_MEMORY ( 0 )

#define ALLOW_ABITRARY_REF_SIZE ( 1 )

#define USE_THREADPOOL ( 1 )

#if ( DO_TESTS == 1 )
	#define BLOCK_SIZE ( 4 ) // 32 Number of threads in a block. should be multiple of wrap-size (16)
	#define GRID_SIZE ( 4 ) // 64 should be similar to number of cores
	#define TILE_WIDTH ( 4 )
	#define CHUNK_SIZE ( BLOCK_SIZE * GRID_SIZE * 2 ) // must be a multiple of BLOCK_SIZE * GRID_SIZE
	
	//// #define QUERY_SIZE ( 61 ) // ( 1024 + 512 ) // size query
	//// #define REFERENCE_SIZE ( 15 ) // (BLOCK_SIZE * GRID_SIZE * 2 * 32) - 11 ) // size reference
#else
	#if ( CARD == 1080)
		/* 64, 128 are good values for GTX 1080i (TILE_WIDTH 16) */
		#define BLOCK_SIZE ( 128 ) // 32 Number of threads in a block. should be multiple of wrap-size (16)
		#define GRID_SIZE ( 64 ) // 64 should be similar to number of cores
		#define TILE_WIDTH ( 48 ) // 16 or 32 (depends on the architecture)
        #define CHUNK_SIZE ( ((size_t)BLOCK_SIZE) * GRID_SIZE * 2048 * 16 )
		// #define CHUNK_SIZE ( ((size_t)BLOCK_SIZE) * GRID_SIZE * 2048 * 16 ) // must be a multiple of BLOCK_SIZE * GRID_SIZE
		
		//// #define QUERY_SIZE ( 1008 ) // ( 1024 + 512 ) // size query 21 * 48 = 1008 / 768
		//// #define REFERENCE_SIZE ( ((size_t)BLOCK_SIZE) * GRID_SIZE * 2048 * 16 * 15 )   // 2048 * 8 ) // 2048 * 5 ) // maximum 4096 * 6
	#else
		/* 64, 128 are good values for GTX 1080i (TILE_WIDTH 16) */
		#define BLOCK_SIZE ( 32 ) // 32 Number of threads in a block. should be multiple of wrap-size (16)
		#define GRID_SIZE ( 64 ) // 64 should be similar to number of cores
		#define TILE_WIDTH ( 16 ) // 16 or 32 (depends on the architecture)
		#define CHUNK_SIZE ( BLOCK_SIZE * GRID_SIZE * 1024 ) // must be a multiple of BLOCK_SIZE * GRID_SIZE
		
		//// #define QUERY_SIZE ( 2044 ) // ( 1024 + 512 ) // size query 21 * 48 = 1008 / 768
		//// #define REFERENCE_SIZE ( BLOCK_SIZE * GRID_SIZE * 2048 * 2 * 1 )   // 2048 * 8 ) // 2048 * 5 ) // maximum 4096 * 6
	#endif
#endif



