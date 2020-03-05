#ifndef MMPRIV2_H
#define MMPRIV2_H

#include "bseq.h"
#include "minimap.h"
#include <assert.h>

#define MM_PARENT_UNSET ( -1 )
#define MM_PARENT_TMP_PRI ( -2 )

#define MM_DBG_NO_KALLOC 0x1
#define MM_DBG_PRINT_QNAME 0x2
#define MM_DBG_PRINT_SEED 0x4
#define MM_DBG_PRINT_ALN_SEQ 0x8

#define MM_SEED_LONG_JOIN ( 1ULL << 40 )
#define MM_SEED_IGNORE ( 1ULL << 41 )
#define MM_SEED_TANDEM ( 1ULL << 42 )
#define MM_SEED_SELF ( 1ULL << 43 )

#define MM_SEED_SEG_SHIFT 48
#define MM_SEED_SEG_MASK ( 0xffULL << ( MM_SEED_SEG_SHIFT ) )

#ifndef kroundup32
#define kroundup32( x )                                                                                                \
    ( --( x ), ( x ) |= ( x ) >> 1, ( x ) |= ( x ) >> 2, ( x ) |= ( x ) >> 4, ( x ) |= ( x ) >> 8,                     \
      ( x ) |= ( x ) >> 16, ++( x ) )
#endif

#define mm_seq4_set( s, i, c ) ( ( s )[ ( i ) >> 3 ] |= ( uint32_t )( c ) << ( ( (i)&7 ) << 2 ) )
#define mm_seq4_get( s, i ) ( ( s )[ ( i ) >> 3 ] >> ( ( (i)&7 ) << 2 ) & 0xf )

#define MALLOC( type, len ) ( (type*)malloc( ( len ) * sizeof( type ) ) )
#define CALLOC( type, len ) ( (type*)calloc( ( len ), sizeof( type ) ) )

#ifdef __cplusplus
extern "C" {
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t
{
    unsigned l, m;
    char* s;
} kstring_t;
#endif

typedef struct
{
    int n_u, n_a;
    uint64_t* u;
    mm128_t* a;
} mm_seg_t;

double cputime( void );
double realtime( void );
long peakrss( void );

void radix_sort_128x( mm128_t* beg, mm128_t* end );
void radix_sort_64( uint64_t* beg, uint64_t* end );
uint32_t ks_ksmall_uint32_t( size_t n, uint32_t arr[], size_t kk );

void mm_sketch( void* km, const char* str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v* p );


const uint64_t* mm_idx_get( const mm_idx_t* mi, uint64_t minier, int* n );
int32_t mm_idx_cal_max_occ( const mm_idx_t* mi, float f );
mm128_t* mm_chain_dp( int max_dist_x, int max_dist_y, int bw, int max_skip, int min_cnt, int min_sc, int is_cdna,
                      int n_segs, int64_t n, mm128_t* a, int* n_u_, uint64_t** _u, void* km );


FILE* mm_split_init( const char* prefix, const mm_idx_t* mi );
mm_idx_t* mm_split_merge_prep( const char* prefix, int n_splits, FILE** fp, uint32_t* n_seq_part );
void mm_split_rm_tmp( const char* prefix, int n_splits );

void mm_err_puts( const char* str );
void mm_err_fwrite( const void* p, size_t size, size_t nitems, FILE* fp );
void mm_err_fread( void* p, size_t size, size_t nitems, FILE* fp );

#ifdef __cplusplus
}
#endif

#endif
