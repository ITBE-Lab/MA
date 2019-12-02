/**
 * @file bwt_large.h
 * @brief BWT-Index Construction for large inputs
 * @author Arne Kutzner
 */
#pragma once

#include "container/qSufSort.h"
#include "util/support.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <errno.h>
#include <stdio.h>
#include <string.h>


/* C++ includes
 */
#include <algorithm>
#include <array>
#include <memory>
/// @endcond

#ifdef USE_MALLOC_WRAPPERS
#include "malloc_wrap.h"
#endif

#ifdef _MSC_VER
#define inline __inline
#endif

#ifdef _MSC_VER
#define __func__ __FUNCTION__
#endif

typedef uint64_t bgint_t;
typedef int64_t sbgint_t;


/**
 * @brief BWT-Index Construction for large inputs
 */
std::tuple<bgint_t, bgint_t, std::vector<bgint_t>, std::vector<unsigned int>>
    EXPORTED bwtLarge( const char *pcFileNamePack );