#include "kthread.h"
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

#define USE_STL
#ifdef USE_STL
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>
#else
#include <pthread.h>
#endif

#if( defined( WIN32 ) || defined( _WIN32 ) ) && defined( _MSC_VER )
#define __sync_fetch_and_add( ptr, addend ) _InterlockedExchangeAdd( (long*)ptr, addend ) // AK: void -> long
#endif

/************
 * kt_for() *
 ************/

struct kt_for_t;

typedef struct
{
    struct kt_for_t* t;
    long i;
} ktf_worker_t;

typedef struct kt_for_t
{
    int n_threads;
    long n;

#ifdef USE_STL
    std::vector<ktf_worker_t> w;
#else
    ktf_worker_t* w;
#endif

    void ( *func )( void*, long, int );
    void* data;
} kt_for_t;

static inline long steal_work( kt_for_t* t )
{
    int i, min_i = -1;
    long k, min = LONG_MAX;
    for( i = 0; i < t->n_threads; ++i )
        if( min > t->w[ i ].i )
            min = t->w[ i ].i, min_i = i;
    k = __sync_fetch_and_add( &t->w[ min_i ].i, t->n_threads ); // do atomic addition
    return k >= t->n ? -1 : k;
} // function

static void* ktf_worker( void* data )
{
    ktf_worker_t* w = (ktf_worker_t*)data; // cast the argument to the appropriate type
    long i;
    for( ;; )
    {
        i = __sync_fetch_and_add( &w->i, w->t->n_threads ); // do atomic addition
        if( i >= w->t->n )
            break;
#ifdef USE_STL
        // dirty ...
        w->t->func( w->t->data, i, (int)( w - &( w->t->w[ 0 ] ) ) );
#else
        w->t->func( w->t->data, i, w - w->t->w );
#endif
    }
    while( ( i = steal_work( w->t ) ) >= 0 )
#ifdef USE_STL
        w->t->func( w->t->data, i, (int)( w - &( w->t->w[ 0 ] ) ) );
#else
        w->t->func( w->t->data, i, w - w->t->w );

#endif
#ifdef USE_STL
    return 0;
#else
    pthread_exit( 0 );
#endif
} // function

void kt_for( int n_threads, void ( *func )( void*, long, int ), void* data, long n )
{
    if( n_threads > 1 )
    {
        int i;
        kt_for_t t;
#ifdef USE_STL
        std::vector<std::thread> tid;
#else
        pthread_t* tid;
#endif
        t.func = func, t.data = data, t.n_threads = n_threads, t.n = n;
#ifndef USE_STL
        t.w = (ktf_worker_t*)calloc( n_threads, sizeof( ktf_worker_t ) );
        tid = (pthread_t*)calloc( n_threads, sizeof( pthread_t ) );
#else
        t.w.resize( n_threads );
#endif
        for( i = 0; i < n_threads; ++i )
            t.w[ i ].t = &t, t.w[ i ].i = i;
        for( i = 0; i < n_threads; ++i )
#ifdef USE_STL
            tid.emplace_back( std::thread( ktf_worker, &t.w[ i ] ) );
#else
            pthread_create( &tid[ i ], 0, ktf_worker, &t.w[ i ] );
#endif
        for( i = 0; i < n_threads; ++i )
#ifdef USE_STL
            tid[ i ].join( );
#else
            pthread_join( tid[ i ], 0 );
#endif
#ifndef USE_STL
        free( tid );
        free( t.w );
#endif
    } // if
    else
    {
        long j;
        for( j = 0; j < n; ++j )
            func( data, j, 0 );
    } // else
} // function

/*****************
 * kt_pipeline() *
 *****************/

struct ktp_t;

typedef struct
{
    struct ktp_t* pl;
    int64_t index;
    int step;
    void* data;
} ktp_worker_t;

struct ktp_t
{
    void* shared;
    void* ( *func )( void*, int, void* );
    int64_t index;
    int n_workers, n_steps;

#ifdef USE_STL
    std::vector<ktp_worker_t> workers;
    std::mutex mutex;
    std::condition_variable cv;
#else
    ktp_worker_t* workers;
    pthread_mutex_t mutex;
    pthread_cond_t cv;
#endif
};

static void* ktp_worker( void* data )
{
    ktp_worker_t* w = (ktp_worker_t*)data;
    ktp_t* p = w->pl;
    while( w->step < p->n_steps )
    {
        // test whether we can kick off the job with this worker
        {
#ifdef USE_STL
            std::lock_guard<std::mutex> lk( p->mutex );
#else
            pthread_mutex_lock( &p->mutex );
#endif

            for( ;; )
            {
                int i;
                // test whether another worker is doing the same step
                for( i = 0; i < p->n_workers; ++i )
                {
                    if( w == &p->workers[ i ] )
                        continue; // ignore itself
                    if( p->workers[ i ].step <= w->step && p->workers[ i ].index < w->index )
                        break;
                }
                if( i == p->n_workers )
                    break; // no workers with smaller indices are doing w->step or the previous steps
#ifdef USE_STL
                {
                    std::unique_lock<std::mutex> ul( p->mutex );
                    ( p->cv ).wait( ul );
                } // scope unique_lock
#else
                pthread_cond_wait( &p->cv, &p->mutex );
#endif
            } // for
#ifndef USE_STL
            pthread_mutex_unlock( &p->mutex );
#endif
        } // scope lock_guard
        // working on w->step
        w->data = p->func( p->shared, w->step, w->step ? w->data : 0 ); // for the first step, input is NULL
        {
            // update step and let other workers know
#ifdef USE_STL
            std::lock_guard<std::mutex> lk( p->mutex );
#else
            pthread_mutex_lock( &p->mutex );
#endif
            w->step = w->step == p->n_steps - 1 || w->data ? ( w->step + 1 ) % p->n_steps : p->n_steps;
            if( w->step == 0 )
                w->index = p->index++;
#ifdef USE_STL
            ( p->cv ).notify_all( );
#else
            pthread_cond_broadcast( &p->cv );
            pthread_mutex_unlock( &p->mutex );
#endif
        } // scope lock_guard
    } // while
#ifdef USE_STL
    return 0;
#else
    pthread_exit( 0 );
#endif
} // function

void kt_pipeline( int n_threads, void* ( *func )(void*, int, void*), void* shared_data, int n_steps )
{
    ktp_t aux;
#ifdef USE_STL
    std::vector<std::thread> tid;
#else
    pthread_t* tid;
#endif
    int i;

    if( n_threads < 1 )
        n_threads = 1;
    aux.n_workers = n_threads;
    aux.n_steps = n_steps;
    aux.func = func;
    aux.shared = shared_data;
    aux.index = 0;

#ifndef USE_STL
    pthread_mutex_init( &aux.mutex, 0 ); // mutex of pipeline
    pthread_cond_init( &aux.cv, 0 ); // condition variable of pipeline
    aux.workers = (ktp_worker_t*)calloc( n_threads, sizeof( ktp_worker_t ) );
#else
    aux.workers.resize( n_threads );
#endif
    for( i = 0; i < n_threads; ++i )
    {
        ktp_worker_t* w = &aux.workers[ i ];
        w->step = 0;
        w->pl = &aux;
        w->data = 0;
        w->index = aux.index++;
    }
#ifndef USE_STL
    tid = (pthread_t*)calloc( n_threads, sizeof( pthread_t ) ); // allocate an array of n_threads
#endif
    for( i = 0; i < n_threads; ++i ) // create n_threads worker
#ifdef USE_STL
        tid.emplace_back( std::thread( ktp_worker, &aux.workers[ i ] ) );
#else
        pthread_create( &tid[ i ], 0, ktp_worker, &aux.workers[ i ] );
#endif

    for( i = 0; i < n_threads; ++i ) // join the threads
#ifdef USE_STL
        tid[ i ].join( );
#else
        pthread_join( tid[ i ], 0 );
#endif


#ifndef USE_STL
    free( tid );
    free( aux.workers );

    pthread_mutex_destroy( &aux.mutex );
    pthread_cond_destroy( &aux.cv );
#endif
}
