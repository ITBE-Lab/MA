#pragma once

#include <stdlib.h>

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

typedef unsigned char ubyte_t;


static void getCounts(const unsigned char *T, int *C, int n, int k, int cs);

static void getBuckets(const int *C, int *B, int k, int end);

static void induceSA(const unsigned char *T, int *SA, int *C, int *B, int n, int k, int cs);

static int sais_main(const unsigned char *T, int *SA, int fs, int n, int k, int cs);

int is_sa(const ubyte_t *T, int *SA, int n);

int is_bwt(ubyte_t *T, int n);