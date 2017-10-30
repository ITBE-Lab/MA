/** 
 * @file debug.h
 * @brief provides some handy defines for debug purposes
 */

#ifndef DEBUG_H
#define DEBUG_H

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 3
#endif //DEBUG_LEVEL

#if DEBUG_LEVEL >= 1
#define DEBUG(x) x
#else //DEBUG_LEVEL
#define DEBUG(x)
#endif //DEBUG_LEVEL

#if DEBUG_LEVEL >= 2
#define DEBUG_2(x) x
#else //DEBUG_LEVEL
#define DEBUG_2(x)
#endif //DEBUG_LEVEL

#if DEBUG_LEVEL >= 3
#define DEBUG_3(x) x
#else //DEBUG_LEVEL
#define DEBUG_3(x)
#endif //DEBUG_LEVEL

#endif