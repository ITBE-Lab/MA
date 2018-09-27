/**
 * @file debug.h
 * @brief provides some handy defines for debug purposes
 */

#ifndef DEBUG_H
#define DEBUG_H

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif // DEBUG_LEVEL

#if DEBUG_LEVEL == 0
// disable instertions
#define NDEBUG
#endif

#if DEBUG_LEVEL >= 1
#define DEBUG( x ) x
#define DEBUG_PARAM( x ) , x
#else // DEBUG_LEVEL
#define DEBUG_PARAM( x )
#define DEBUG( x )
#endif // DEBUG_LEVEL

#if DEBUG_LEVEL >= 2
#define DEBUG_2( x ) x
#else // DEBUG_LEVEL
#define DEBUG_2( x )
#endif // DEBUG_LEVEL

#if DEBUG_LEVEL >= 3
#define DEBUG_3( x ) x
#else // DEBUG_LEVEL
#define DEBUG_3( x )
#endif // DEBUG_LEVEL

#endif

#include <assert.h>