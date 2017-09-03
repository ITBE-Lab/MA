#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <errno.h>

#include <boost/log/trivial.hpp>
#include "QSufSort.h"

/* C++ includes
 */
#include <algorithm>
#include <vector>
#include <memory> 
#include <array>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

#ifdef _MSC_VER
	#define inline __inline
#endif

#ifdef _MSC_VER
	#define __func__ __FUNCTION__
#endif

typedef uint64_t bgint_t;
typedef int64_t sbgint_t;

#define ALPHABET_SIZE				4
#define BIT_PER_CHAR				2
#define CHAR_PER_WORD				16
#define CHAR_PER_BYTE				4

#define BITS_IN_WORD 32
#define BITS_IN_BYTE 8
#define BYTES_IN_WORD 4

#define ALL_ONE_MASK 0xFFFFFFFF
#define DNA_OCC_CNT_TABLE_SIZE_IN_WORD	65536

#define BITS_PER_OCC_VALUE			16
#define OCC_VALUE_PER_WORD			2
#define OCC_INTERVAL				256
#define OCC_INTERVAL_MAJOR			65536

#define TRUE    1
#define FALSE   0

#define BWTINC_INSERT_SORT_NUM_ITEM 7

#define MIN_AVAILABLE_WORD 0x10000

#define average(value1, value2)					( ((value1) & (value2)) + ((value1) ^ (value2)) / 2 )

#define med3(a, b, c)							( a<b ? (b<c ? b : a<c ? c : a) : (b>c ? b : a>c ? c : a))
#define swap(a, b, t);							t = a; a = b; b = t;
#define truncateLeft(value, offset)				( (value) << (offset) >> (offset) )
#define truncateRight(value, offset)			( (value) >> (offset) << (offset) )
#define DNA_OCC_SUM_EXCEPTION(sum)			((sum & 0xfefefeff) == 0)

typedef struct BWT {
	bgint_t textLength;					// length of the text
	bgint_t inverseSa0;					// SA-1[0] representing the substring of all preceding characters
	bgint_t *cumulativeFreq;			// cumulative frequency
	unsigned int *bwtCode;				// BWT code
	unsigned int *occValue;				// Occurrence values stored explicitly
	bgint_t *occValueMajor;				// Occurrence values stored explicitly
	unsigned int *decodeTable;			// For decoding BWT by table lookup
	bgint_t bwtSizeInWord;				// Temporary variable to hold the memory allocated
	bgint_t occSizeInWord;				// Temporary variable to hold the memory allocated
	bgint_t occMajorSizeInWord;			// Temporary variable to hold the memory allocated
} BWT;

typedef struct BWTInc {
	BWT *bwt;
	unsigned int numberOfIterationDone;
	bgint_t *cumulativeCountInCurrentBuild;
	bgint_t availableWord;
	bgint_t buildSize;
	bgint_t initialMaxBuildSize;
	bgint_t incMaxBuildSize;
	unsigned int firstCharInLastIteration;
	unsigned int *workingMemory;
	unsigned int *packedText;
	unsigned char *textBuffer;
	unsigned int *packedShift;
} BWTInc;

/* New function for interfacing with the FM_Index class
 */
 std::tuple<
            bgint_t,
            bgint_t,
            std::vector<bgint_t>,
            std::vector<unsigned int >
        > bwtLarge( const char *pcFileNamePack );