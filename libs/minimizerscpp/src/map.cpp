#include "bseq.h"
#include "kalloc.h"
#include "khash.h"
#include "kthread.h"
#include "kvec.h"
#include "mmpriv.h"
#include "sdust.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

static int mm_dust_minier( void* km, int n, mm128_t* a, int l_seq, const char* seq, int sdust_thres )
{
    int n_dreg, j, k, u = 0;
    const uint64_t* dreg;
    sdust_buf_t* sdb;
    if( sdust_thres <= 0 )
        return n;
    sdb = sdust_buf_init( km );
    dreg = sdust_core( (const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb );
    for( j = k = 0; j < n; ++j )
    { // squeeze out minimizers that significantly overlap with LCRs
        int32_t qpos = (uint32_t)a[ j ].y >> 1, span = a[ j ].x & 0xff;
        int32_t s = qpos - ( span - 1 ), e = s + span;
        while( u < n_dreg && (int32_t)dreg[ u ] <= s )
            ++u;
        if( u < n_dreg && ( int32_t )( dreg[ u ] >> 32 ) < e )
        {
            int v, l = 0;
            for( v = u; v < n_dreg && ( int32_t )( dreg[ v ] >> 32 ) < e; ++v )
            { // iterate over LCRs overlapping this minimizer
                int ss = s > ( int32_t )( dreg[ v ] >> 32 ) ? s : dreg[ v ] >> 32;
                int ee = e < (int32_t)dreg[ v ] ? e : (uint32_t)dreg[ v ];
                l += ee - ss;
            }
            if( l <= span >> 1 )
                a[ k++ ] = a[ j ]; // keep the minimizer if less than half of it falls in masked region
        }
        else
            a[ k++ ] = a[ j ];
    }
    sdust_buf_destroy( sdb );
    return k; // the new size
}

/*
 * @param mi         minimap2 index
 * @param _n_m       number of returned mathces
 */
static void collect_minimizers( void* km, const mm_mapopt_t* opt, const mm_idx_t* mi, int n_segs, const int* qlens,
                                const char** seqs, mm128_v* mv, void ( *mm_filter )( mm128_v*, void* ),
                                void* pFilterArg )
{
    int i, n, sum = 0;
    mv->n = 0;
    for( i = n = 0; i < n_segs; ++i )
    {
        size_t j;
        mm_sketch( km, seqs[ i ], qlens[ i ], mi->w, mi->k, i, mi->flag & MM_I_HPC, mv );
        ( *mm_filter )( mv, pFilterArg );
        for( j = n; j < mv->n; ++j )
            mv->a[ j ].y += sum << 1;
        if( opt->sdust_thres > 0 ) // mask low-complexity minimizers
            mv->n = n + mm_dust_minier( km, mv->n - n, mv->a + n, qlens[ i ], seqs[ i ], opt->sdust_thres );
        sum += qlens[ i ], n = mv->n;
    }
}

#include "ksort.h"
#define heap_lt( a, b ) ( ( a ).x > ( b ).x )
KSORT_INIT( heap, mm128_t, heap_lt )

typedef struct
{
    uint32_t n;
    uint32_t q_pos, q_span;
    uint32_t seg_id : 31, is_tandem : 1;
    const uint64_t* cr;
} mm_match_t;

size_t uiXSkipped = 0;
size_t uiXRetrived = 0;
/**
 *
 * @param rep_len returns the number of nucleotides that are not covered by minimizers
 */
static mm_match_t* collect_matches( void* km, int* _n_m, int max_occ, const mm_idx_t* mi, const mm128_v* mv,
                                    int64_t* n_a, int* rep_len, int* n_mini_pos, uint64_t** mini_pos )
{
    int rep_st = 0;
    int rep_en = 0;
    mm_match_t* m;
    *n_mini_pos = 0;
    *mini_pos = (uint64_t*)kmalloc( km, mv->n * sizeof( uint64_t ) );
    m = (mm_match_t*)kmalloc( km, mv->n * sizeof( mm_match_t ) );

    *rep_len = 0;
    int n_m = 0;
    *n_a = 0;
    for( size_t i = 0; i < mv->n; ++i )
    {
        const uint64_t* cr;
        mm128_t* p = &mv->a[ i ];
        uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
        int t;
        cr = mm_idx_get( mi, p->x >> 8, &t );
        // fprintf(stdout, "t: %d max_occ: %d\n", t, max_occ);
        if( t >= max_occ )
        {
            uiXSkipped += t;
            // fprintf(stdout, "XXXX\n");
            int end = ( q_pos >> 1 ) + 1;
            int start = end - q_span;
            if( start > rep_en )
            {
                *rep_len += rep_en - rep_st;
                rep_st = start;
                rep_en = end;
            }
            else
                rep_en = end;
        }
        else
        {
            uiXRetrived += t;
            // fprintf(stdout, "found entry\n");
            mm_match_t* q = &m[ n_m++ ];
            q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
            q->is_tandem = 0;
            if( i > 0 && p->x >> 8 == mv->a[ i - 1 ].x >> 8 )
                q->is_tandem = 1;
            if( i < mv->n - 1 && p->x >> 8 == mv->a[ i + 1 ].x >> 8 )
                q->is_tandem = 1;
            *n_a += q->n;
            ( *mini_pos )[ ( *n_mini_pos )++ ] = (uint64_t)q_span << 32 | q_pos >> 1;
        }
    }
    *rep_len += rep_en - rep_st;
    *_n_m = n_m;
    return m;
}
#if 0
/**
 *
 * @param rep_len returns the number of nucleotides that are not covered by minimizers
 */
static mm_match_t* collect_matches_adaptive_filter( void* km, int* _n_m, int max_occ, const mm_idx_t* mi,
                                                    const mm128_v* mv, int64_t* n_a, int* rep_len, int* n_mini_pos,
                                                    uint64_t** mini_pos )
{
    int max_rep_len = mi->k + mi->w;
    int rep_st = -1;
    int rep_en = -1;
    mm_match_t* m;
    *n_mini_pos = 0;
    *mini_pos = (uint64_t*)kmalloc( km, mv->n * sizeof( uint64_t ) );
    m = (mm_match_t*)kmalloc( km, mv->n * sizeof( mm_match_t ) );

    *rep_len = 0;
    int n_m = 0;
    *n_a = 0;

    int best_t = 0;
    const uint64_t* best_cr = NULL;
    size_t best_i = 0;

    for( size_t j = 0; j < mv->n; ++j )
    {
        // so that we can overwrite i without affecting the current loop position
        size_t i = j;
        const uint64_t* cr;
        mm128_t* p = &mv->a[ i ];
        uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
        int t;
        cr = mm_idx_get( mi, p->x >> 8, &t );
        if(t == 0)
            continue;
        assert( cr != NULL );
        // if the minimizer is too ambiguous
        if( t >= max_occ )
        {
            // get the end of the current minimizer
            int end = ( q_pos >> 1 ) + 1;
            // get the start of the current minimizer
            int start = end - q_span;
            // if we are starting a new hole on the query
            // (hole := region without seeds due to ambiguity)
            if( start > rep_en )
            {
                if( rep_en > rep_st )
                    *rep_len += rep_en - rep_st;
                // set the start and end of the hole
                rep_st = start;
                rep_en = end;

                // since this is the first seed of the hole set this as the best skipped seed
                best_i = i;
                best_cr = cr;
                best_t = t;

                // do not extract seeds for this minimizer
                // fprintf( stdout, "a\n" );
                continue;
            }
            // if we are extending a hole
            else
            {
                // adjust the end of the hole
                rep_en = end;
                //fprintf( stdout, "e\n" );

                // if the hole has become longer than allowed
                if( start - rep_st > max_rep_len )
                {
                    // extract the best minimizer in the hole
                    assert( best_i < j );
                    i = best_i;
                    cr = best_cr;
                    t = best_t;
                    if(t == 0)
                        continue;

                    // find the best minimizer after the extracted one:
                    best_t = 0;
                    for( size_t k = i + 1; k <= j; ++k )
                    {
                        const uint64_t* cr;
                        mm128_t* p = &mv->a[ k ];
                        int t;
                        cr = mm_idx_get( mi, p->x >> 8, &t );
                        if(t == 0)
                            continue;
                        assert( cr != NULL );
                        if( best_t == 0 || t < best_t )
                        {
                            // remember this minimizer
                            best_i = k;
                            best_cr = cr;
                            best_t = t;
                        } // if
                    } // for

                    assert( cr != NULL );
                    p = &mv->a[ i ];
                    q_pos = (uint32_t)p->y, q_span = p->x & 0xff;

                    // increment rep_len by the length of the gap before the best minimizer in the hole
                    int rep_en_curr = ( q_pos >> 1 ) + 1 - q_span;
                    if( rep_en_curr > rep_st )
                        *rep_len += rep_en_curr - rep_st;
                    // set the start of the hole to the end of currently extracted minimizer
                    rep_st = ( q_pos >> 1 ) + 1;
                    // NO continue: do extract the current minimizer (overwritten with the best minimizer in the hole)
                    //fprintf( stdout, "d\n" );
                } // if
                // the hole is not too long yet...
                else
                {
                    // if the current minimizer is the best one we have seen within this hole
                    if( t < best_t )
                    {
                        // remember this minimizer
                        best_i = i;
                        best_cr = cr;
                        best_t = t;
                        //fprintf( stdout, "c\n" );
                    } // if
                    // don't extract the current minimizer
                    //fprintf( stdout, "b\n" );
                    continue;
                } // if
            } // else
        } // if
        // else
        //   fprintf( stdout, "f\n" );

        if( t > 0 )
        {
            assert( cr != NULL );

            mm_match_t* q = &m[ n_m++ ];
            q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
            q->is_tandem = 0;
            if( i > 0 && p->x >> 8 == mv->a[ i - 1 ].x >> 8 )
                q->is_tandem = 1;
            if( i < mv->n - 1 && p->x >> 8 == mv->a[ i + 1 ].x >> 8 )
                q->is_tandem = 1;
            *n_a += q->n;
            ( *mini_pos )[ ( *n_mini_pos )++ ] = (uint64_t)q_span << 32 | q_pos >> 1;
        } // if
    } // for
    if( rep_en > rep_st )
        *rep_len += rep_en - rep_st;
    *_n_m = n_m;
    return m;
}
#endif

static inline int skip_seed( int flag, uint64_t r, const mm_match_t* q, const char* qname, int qlen, const mm_idx_t* mi,
                             int* is_self )
{
    *is_self = 0;
    if( qname && ( flag & ( MM_F_NO_DIAG | MM_F_NO_DUAL ) ) )
    {
        // fprintf(stdout, "markus liegt falsch\n");
        const mm_idx_seq_t* s = &mi->seq[ r >> 32 ];
        int cmp;
        cmp = strcmp( qname, s->name );
        if( ( flag & MM_F_NO_DIAG ) && cmp == 0 && (int)s->len == qlen )
        {
            if( (uint32_t)r >> 1 == ( q->q_pos >> 1 ) )
                return 1; // avoid the diagnonal anchors
            if( ( r & 1 ) == ( q->q_pos & 1 ) )
                *is_self = 1; // this flag is used to avoid spurious extension on self chain
        }
        if( ( flag & MM_F_NO_DUAL ) && cmp > 0 ) // all-vs-all mode: map once
            return 1;
    }
    if( flag & ( MM_F_FOR_ONLY | MM_F_REV_ONLY ) )
    {
        // fprintf(stdout, "markus liegt falsch\n");
        if( ( r & 1 ) == ( q->q_pos & 1 ) )
        { // forward strand
            if( flag & MM_F_REV_ONLY )
                return 1;
        }
        else
        {
            if( flag & MM_F_FOR_ONLY )
                return 1;
        }
    }
    return 0;
}

static mm128_t* collect_seed_hits_heap( void* km, const mm_mapopt_t* opt, int max_occ, const mm_idx_t* mi,
                                        const char* qname, const mm128_v* mv, int qlen, int64_t* n_a, int* rep_len,
                                        int* n_mini_pos, uint64_t** mini_pos )
{
    int i, n_m, heap_size = 0;
    int64_t j, n_for = 0, n_rev = 0;
    mm_match_t* m;
    mm128_t *a, *heap;

    m = collect_matches( km, &n_m, max_occ, mi, mv, n_a, rep_len, n_mini_pos, mini_pos );

    heap = (mm128_t*)kmalloc( km, n_m * sizeof( mm128_t ) );
    a = (mm128_t*)kmalloc( km, *n_a * sizeof( mm128_t ) );

    for( i = 0, heap_size = 0; i < n_m; ++i )
    {
        if( m[ i ].n > 0 )
        {
            heap[ heap_size ].x = m[ i ].cr[ 0 ];
            heap[ heap_size ].y = (uint64_t)i << 32;
            ++heap_size;
        }
    }
    ks_heapmake_heap( heap_size, heap );
    while( heap_size > 0 )
    {
        mm_match_t* q = &m[ heap->y >> 32 ];
        mm128_t* p;
        uint64_t r = heap->x;
        int32_t is_self, rpos = (uint32_t)r >> 1;
        if( !skip_seed( opt->flag, r, q, qname, qlen, mi, &is_self ) )
        {
            if( ( r & 1 ) == ( q->q_pos & 1 ) )
            { // forward strand
                p = &a[ n_for++ ];
                p->x = ( r & 0xffffffff00000000ULL ) | rpos;
                p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
            }
            else
            { // reverse strand
                p = &a[ ( *n_a ) - ( ++n_rev ) ];
                p->x = 1ULL << 63 | ( r & 0xffffffff00000000ULL ) | rpos;
                p->y = (uint64_t)q->q_span << 32 | ( qlen - ( ( q->q_pos >> 1 ) + 1 - q->q_span ) - 1 );
            }
            p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
            if( q->is_tandem )
                p->y |= MM_SEED_TANDEM;
            if( is_self )
                p->y |= MM_SEED_SELF;
        }
        // update the heap
        if( (uint32_t)heap->y < q->n - 1 )
        {
            ++heap[ 0 ].y;
            heap[ 0 ].x = m[ heap[ 0 ].y >> 32 ].cr[ (uint32_t)heap[ 0 ].y ];
        }
        else
        {
            heap[ 0 ] = heap[ heap_size - 1 ];
            --heap_size;
        }
        ks_heapdown_heap( 0, heap_size, heap );
    }
    kfree( km, m );
    kfree( km, heap );

    // reverse anchors on the reverse strand, as they are in the descending order
    for( j = 0; j<n_rev>> 1; ++j )
    {
        mm128_t t = a[ ( *n_a ) - 1 - j ];
        a[ ( *n_a ) - 1 - j ] = a[ ( *n_a ) - ( n_rev - j ) ];
        a[ ( *n_a ) - ( n_rev - j ) ] = t;
    }
    if( *n_a > n_for + n_rev )
    {
        memmove( a + n_for, a + ( *n_a ) - n_rev, n_rev * sizeof( mm128_t ) );
        *n_a = n_for + n_rev;
    }
    return a;
}

static mm128_t* collect_seed_hits( void* km, const mm_mapopt_t* opt, int max_occ, const mm_idx_t* mi, const char* qname,
                                   const mm128_v* mv, int qlen, int64_t* n_a, int* rep_len, int* n_mini_pos,
                                   uint64_t** mini_pos )
{
    int i, n_m;
    mm_match_t* m;
    mm128_t* a;
    m = collect_matches( km, &n_m, max_occ, mi, mv, n_a, rep_len, n_mini_pos, mini_pos );
    a = (mm128_t*)kmalloc( km, *n_a * sizeof( mm128_t ) ); // result vector holding seeds
    for( i = 0, *n_a = 0; i < n_m; ++i )
    {
        mm_match_t* q = &m[ i ];
        const uint64_t* r = q->cr;
        uint32_t k;
        for( k = 0; k < q->n; ++k )
        {
            int32_t is_self, rpos = (uint32_t)r[ k ] >> 1;
            mm128_t* p;
            if( skip_seed( opt->flag, r[ k ], q, qname, qlen, mi, &is_self ) )
                continue;
            p = &a[ ( *n_a )++ ];
            if( ( r[ k ] & 1 ) == ( q->q_pos & 1 ) )
            { // forward strand
                p->x = ( r[ k ] & 0xffffffff00000000ULL ) | rpos;
                p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
            }
            else
            { // reverse strand
                p->x = 1ULL << 63 | ( r[ k ] & 0xffffffff00000000ULL ) | rpos;
                p->y = (uint64_t)q->q_span << 32 | ( qlen - ( ( q->q_pos >> 1 ) + 1 - q->q_span ) - 1 );
            }
            p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
            if( q->is_tandem )
                p->y |= MM_SEED_TANDEM;
            if( is_self )
                p->y |= MM_SEED_SELF;
        }
    }
    kfree( km, m );
    // radix_sort_128x( a, a + ( *n_a ) ); // sorting the result vector
    return a;
}

/**
 * by markus
 */
mm128_t* collect_seeds( const mm_idx_t* mi, int n_segs, const int* qlens, const char** seqs, mm_tbuf_t* b,
                        const mm_mapopt_t* opt, const char* qname, int64_t* n_a, void ( *mm_filter )( mm128_v*, void* ),
                        void* pFilterArg )
{
    int i, rep_len, qlen_sum, n_mini_pos;
    uint64_t* mini_pos;
    mm128_v mv = {0, 0, 0};

    for( i = 0, qlen_sum = 0; i < n_segs; ++i )
        qlen_sum += qlens[ i ];

    if( qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG )
    {
        *n_a = -1;
        fprintf( stderr, "%d %d\n", qlen_sum, n_segs );
        return NULL;
    } // if
    // fprintf(stderr, "qlen_sum = %d\n", qlen_sum);

    mm128_t* a;

    collect_minimizers( b->km, opt, mi, n_segs, qlens, seqs, &mv, mm_filter, pFilterArg );
    // fprintf(stderr, "n=%lu m=%lu\n", mv.n, mv.m);
    if( opt->flag & MM_F_HEAP_SORT )
        a = collect_seed_hits_heap( b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, n_a, &rep_len, &n_mini_pos,
                                    &mini_pos );
    else
        a = collect_seed_hits( b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, n_a, &rep_len, &n_mini_pos,
                               &mini_pos );
#if 0
    if(rep_len != 0)
        fprintf( stderr, "RS\t%d\n", rep_len );
    for( i = 0; i < *n_a; ++i )
        fprintf( stderr, "SD\t%s\t%d\t%c\t%d\t%d\t%d\n", mi->seq[ a[ i ].x << 1 >> 33 ].name, (int32_t)a[ i ].x,
                 "+-"[ a[ i ].x >> 63 ], (int32_t)a[ i ].y, ( int32_t )( a[ i ].y >> 32 & 0xff ),
                 i == 0
                     ? 0
                     : ( (int32_t)a[ i ].y - (int32_t)a[ i - 1 ].y ) - ( (int32_t)a[ i ].x - (int32_t)a[ i - 1 ].x ) );
#endif
    kfree( b->km, mv.a );
    kfree( b->km, mini_pos );
    return a;
}
