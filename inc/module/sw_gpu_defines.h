#pragma once

#define GAP_INIT ( 4 )
#define GAP_EXT ( 1 )
#define GAP_EXT_HORIZONAL ( 1 )

#define NUM_OF_SYMBOLS ( 4 )

#define CUERR { cudaError_t err;			\
				if ((err = cudaGetLastError()) != cudaSuccess) {		\
  					printf("CUDA error: %s : %s, line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); exit(0); } \
			  } 
