#pragma once

#include <stdlib.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

typedef unsigned char ubyte_t;

int is_sa(const ubyte_t *T, int *SA, int n);

int is_bwt(ubyte_t *T, int n);