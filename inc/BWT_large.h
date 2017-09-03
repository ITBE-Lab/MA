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


/* New function for interfacing with the FM_Index class
 */
 std::tuple<
            bgint_t,
            bgint_t,
            std::vector<bgint_t>,
            std::vector<unsigned int >
        > bwtLarge( const char *pcFileNamePack );