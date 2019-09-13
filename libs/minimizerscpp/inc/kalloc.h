#ifndef _KALLOC_H_
#define _KALLOC_H_

#include <stddef.h> /* for size_t */

#ifdef __cplusplus
extern "C" {
#endif

#define kmalloc(km, size) malloc(size)
#define krealloc(km, ptr, size) realloc(ptr, size)
#define kcalloc(km, count, size) calloc(count, size)
#define kfree(km, ptr) free(ptr)

//void *kmalloc(void *km, size_t size);
//void *krealloc(void *km, void *ptr, size_t size);
//void *kcalloc(void *km, size_t count, size_t size);
//void kfree(void *km, void *ptr);


#ifdef __cplusplus
}
#endif

#endif
