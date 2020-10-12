
#include "ma/container/nucSeq.h"
#include "ma/container/pack.h"
#include "ma/container/segment.h"
#include "ma/module/stripOfConsideration.h"
#include "ms/util/parameter.h"
#ifndef __cplusplus
#define __cplusplus TRUE
#endif
#ifdef WITH_ZLIB
#include "minimap.h"
#endif

#pragma once

namespace minimizer
{
#ifdef WITH_ZLIB

inline void mm_filter_none( mm128_t*, size_t&, void* )
{} // method

/**
 * @brief cpp wrapper for the minimap index
 * @details
 * at the moment multipart indices are not supported
 */
class Index : public libMS::Container
{
  private:
    mm_idx_t* pData;
    mm_idxopt_t xOptions;
    mm_mapopt_t xMapOpt;

    class IndexReader
    {
      public:
        mm_idx_reader_t* idx_rdr;

        IndexReader( std::string sIndexName, mm_idxopt_t* pOptions, const char* pcOut )
        {
            idx_rdr = mm_idx_reader_open( sIndexName.c_str( ), pOptions, pcOut );
        } // constructor

        ~IndexReader( )
        {
            if( idx_rdr != nullptr )
                mm_idx_reader_close( idx_rdr );
        } // destructor
    }; // class

    mm_tbuf_t* mm_tbuf_init( void )
    {
        mm_tbuf_t* b;
        b = (mm_tbuf_t*)calloc( 1, sizeof( mm_tbuf_t ) );
        if( !( mm_dbg_flag & 1 ) )
            b->km = km_init( );
        return b;
    }

    void mm_tbuf_destroy( mm_tbuf_t* b )
    {
        if( b == 0 )
            return;
        km_destroy( b->km );
        free( b );
    }

    void mm_mapopt_init( mm_mapopt_t* opt )
    {
        memset( opt, 0, sizeof( mm_mapopt_t ) );
        opt->seed = 11;
        opt->mid_occ_frac = 2e-4f;
        opt->sdust_thres = 0; // no SDUST masking

        opt->min_cnt = 3;
        opt->min_chain_score = 40;
        opt->bw = 500;
        opt->max_gap = 5000;
        opt->max_gap_ref = -1;
        opt->max_chain_skip = 25;

        // opt->mask_level = 0.5f;
        // opt->pri_ratio = 0.8f;
        // opt->best_n = 5;

        // opt->max_join_long = 20000;
        // opt->max_join_short = 2000;
        // opt->min_join_flank_sc = 1000;
        // opt->min_join_flank_ratio = 0.5f;

        opt->a = 2, opt->b = 4, opt->q = 4, opt->e = 2, opt->q2 = 24, opt->e2 = 1;
        opt->sc_ambi = 1;
        opt->zdrop = 400, opt->zdrop_inv = 200;
        opt->end_bonus = -1;
        opt->min_dp_max = opt->min_chain_score * opt->a;
        opt->min_ksw_len = 200;
        opt->anchor_ext_len = 20, opt->anchor_ext_shift = 6;
        opt->max_clip_ratio = 1.0f;
        opt->mini_batch_size = 500000000;

        opt->pe_ori = 0; // FF
        opt->pe_bonus = 33;
    }

    void mm_mapopt_update( mm_mapopt_t* opt, const mm_idx_t* mi )
    {
        if( ( opt->flag & MM_F_SPLICE_FOR ) || ( opt->flag & MM_F_SPLICE_REV ) )
            opt->flag |= MM_F_SPLICE;
        // if( opt->mid_occ <= 0 )
        //    opt->mid_occ = mm_idx_cal_max_occ( mi, opt->mid_occ_frac );
        // if( opt->mid_occ < opt->min_mid_occ )
        //    opt->mid_occ = opt->min_mid_occ;
        opt->mid_occ = 200;
        // if( mm_verbose >= 3 )
        //    fprintf( stderr, "[M::%s::%.3f*%.2f] mid_occ = %d\n", __func__, realtime( ) - mm_realtime0,
        //             cputime( ) / ( realtime( ) - mm_realtime0 ), opt->mid_occ );
    } // method

    void initOptions( )
    {
        mm_mapopt_init( &xMapOpt );
        xMapOpt.flag |= MM_F_CIGAR; // | MM_F_HEAP_SORT; // perform alignment
        // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        mm_mapopt_update( &xMapOpt, pData );
    } // method

  public:
    /**
     * @brief Opens in Index if index filename is given or creates index if fasta(q) file is given
     */
    Index( const ParameterSetManager& rParameters, std::string sIndexName )
        : xOptions{/* .k = */ rParameters.getSelected( )->xMinimizerK->get( ),
                   /* .w = */ rParameters.getSelected( )->xMinimizerW->get( ),
                   /* .flag = */ rParameters.getSelected( )->xMinimizerFlag->get( ),
                   /* .bucket_bits = */ rParameters.getSelected( )->xMinimizerBucketBits->get( ),
                   /* .mini_batch_size = */ rParameters.getSelected( )->xMinimizerMiniBatchSize->get( ),
                   /* .batch_size = */ rParameters.getSelected( )->xMinimizerBatchSize->get( )}
    {
        this->xOptions.k = rParameters.getSelected( )->xMinimizerK->get( );
        this->xOptions.w = rParameters.getSelected( )->xMinimizerW->get( );
        this->xOptions.flag = rParameters.getSelected( )->xMinimizerFlag->get( );
        this->xOptions.bucket_bits = rParameters.getSelected( )->xMinimizerBucketBits->get( );
        this->xOptions.mini_batch_size = rParameters.getSelected( )->xMinimizerMiniBatchSize->get( );
        this->xOptions.batch_size = rParameters.getSelected( )->xMinimizerBatchSize->get( );

        IndexReader xIndexReader( sIndexName, &xOptions, NULL );
        if( xIndexReader.idx_rdr == 0 )
            throw std::runtime_error( "failed to open file" + sIndexName );
        pData = mm_idx_reader_read( xIndexReader.idx_rdr, (int)rParameters.getNumThreads( ) );
        if( ( xOptions.flag & MM_F_CIGAR ) && ( pData->flag & MM_I_NO_SEQ ) )
            throw std::runtime_error( "the prebuilt index doesn't contain sequences" );
        initOptions( );
#if DEBUG_LEVEL > 0
        mm_idx_stat( pData );
#endif
        // make sure that there are no multipart indices
        if( mm_idx_reader_read( xIndexReader.idx_rdr, (int)rParameters.getNumThreads( ) ) != 0 )
            throw std::runtime_error( "index comes in multiple parts" );
    } // constructor

    Index( const ParameterSetManager& rParameters, std::vector<std::string> vContigs,
           std::vector<std::string> vContigsNames )
        : xOptions{/* .k = */ rParameters.getSelected( )->xMinimizerK->get( ),
                   /* .w = */ rParameters.getSelected( )->xMinimizerW->get( ),
                   /* .flag = */ rParameters.getSelected( )->xMinimizerFlag->get( ),
                   /* .bucket_bits = */ rParameters.getSelected( )->xMinimizerBucketBits->get( ),
                   /* .mini_batch_size = */ rParameters.getSelected( )->xMinimizerMiniBatchSize->get( ),
                   /* .batch_size = */ rParameters.getSelected( )->xMinimizerBatchSize->get( )}
    {
        xOptions.k = rParameters.getSelected( )->xMinimizerK->get( );
        xOptions.w = rParameters.getSelected( )->xMinimizerW->get( );
        xOptions.flag = rParameters.getSelected( )->xMinimizerFlag->get( );
        xOptions.bucket_bits = rParameters.getSelected( )->xMinimizerBucketBits->get( );
        xOptions.mini_batch_size = rParameters.getSelected( )->xMinimizerMiniBatchSize->get( );
        xOptions.batch_size = rParameters.getSelected( )->xMinimizerBatchSize->get( );


        std::vector<const char*> seq;
        std::vector<const char*> name;
        for( auto& sContig : vContigs )
            seq.push_back( sContig.data( ) );
        for( auto& sName : vContigsNames )
            name.push_back( sName.data( ) );

        pData = mm_idx_str( xOptions.w, xOptions.k, xOptions.flag & MM_I_HPC, xOptions.bucket_bits,
                            (int)vContigs.size( ), &seq[ 0 ], &name[ 0 ] );
        initOptions( );
#if DEBUG_LEVEL > 0
        mm_idx_stat( pData );
        std::cout << "max Occ:" << xMapOpt.mid_occ << std::endl;
#endif
    } // constructor

    ~Index( )
    {
        mm_idx_destroy( pData );
    } // destructor

    void dump( std::string sIndexName )
    {
        mm_idx_stat( pData );
        if( mm_idx_dump_name( sIndexName.c_str( ), pData ) != 0 )
            throw std::runtime_error( "failed to open file" + sIndexName );
    } // method


    std::shared_ptr<libMA::Seeds> seed_one( const char* sSeq, const int iSize, bool bRectangular,
                                            std::shared_ptr<libMA::Pack> pPack,
                                            void ( *mm_filter )( mm128_t*, size_t&, void* ), void* pFilterArg )
    {
        auto pRet = std::make_shared<libMA::Seeds>( );
        int64_t n_a = 0;

        mm_tbuf_t* tbuf = mm_tbuf_init( );
        // const char* seqs = (const char*)&pQuery->pxSequenceRef;
        mm128_t* a = collect_seeds( pData, 1, &iSize, &sSeq, tbuf, &xMapOpt, 0, &n_a, mm_filter, pFilterArg );
        if( a != NULL )
        {
            for( int64_t uiI = 0; uiI < n_a; uiI++ )
            {
                if( ( a[ uiI ].x >> 63 ) == 0 )
                {
                    // std::cout << "q_span: " << ( int32_t )( a[ uiI ].y >> 32 ) << std::endl;
                    // fprintf( stderr, "SD\t%d\t%c\t%d\t%d\t%d\n", (int32_t)a[ uiI ].x, "+-"[ a[ uiI ].x >> 63 ],
                    //         (int32_t)a[ uiI ].y, ( int32_t )( a[ uiI ].y >> 32 & 0xff ),
                    //         uiI == 0 ? 0
                    //                  : ( (int32_t)a[ uiI ].y - (int32_t)a[ uiI - 1 ].y ) -
                    //                        ( (int32_t)a[ uiI ].x - (int32_t)a[ uiI - 1 ].x ) );
                    uint64_t uiQSpan = ( int32_t )( a[ uiI ].y >> 32 & 0xff );
                    // assert( (int32_t)a[ uiI ].y + 1 >= (int32_t)uiQSpan );
                    // assert( (int32_t)a[ uiI ].x + 1 >= (int32_t)uiQSpan );
                    pRet->emplace_back(
                        /* minimap uses 1-based index i use 0-based indices...
                           a[ uiI ].y = last base of minimizer */
                        ( (int32_t)a[ uiI ].y ) + 1 - uiQSpan,
                        uiQSpan,
                        /* reference position is encoded in the lower 32 bits of x */
                        ( (int32_t)a[ uiI ].x ) + 1 - uiQSpan + pPack->startOfSequenceWithId( a[ uiI ].x << 1 >> 33 ),
                        true // minimap code : "+-"[a[i].x>>63]
                    );
                    assert( !pPack->bPositionIsOnReversStrand( pRet->back( ).start_ref( ) ) );
                    assert( !pPack->bPositionIsOnReversStrand( pRet->back( ).end_ref( ) - 1 ) );
                } // if
                else
                {
                    uint64_t uiQSpan = ( int32_t )( a[ uiI ].y >> 32 & 0xff );
                    pRet->emplace_back(
                        /* minimap uses 1-based index i use 0-based indices...
                           a[ uiI ].y = last base of minimizer */
                        iSize - 1 - ( (int32_t)a[ uiI ].y ),
                        uiQSpan,
                        /* reference position is encoded in the lower 32 bits of x */
                        (int32_t)a[ uiI ].x + pPack->startOfSequenceWithId( a[ uiI ].x << 1 >> 33 ),
                        false // minimap code : "+-"[a[i].x>>63]
                    );
                    assert( !pPack->bPositionIsOnReversStrand( pRet->back( ).start_ref( ) - 1 ) );
                    assert( !pPack->bPositionIsOnReversStrand( pRet->back( ).start_ref( ) - pRet->back( ).size( ) ) );
                } // else
            } // for
            kfree( tbuf->km, a );
        } // if
        else if( n_a > 0 )
            throw std::runtime_error( "minimizer vector is empty" );
        mm_tbuf_destroy( tbuf );

        // delta values need to be set in case we want to compute a SoC...
        // @todo this should be done by the SoC module but is done as well in the binary seeding module so i added it
        // here for now
        for( auto& rSeed : *pRet )
            libMA::ExtractSeeds::setDeltaOfSeed( rSeed, iSize, *pPack, !bRectangular );

        return pRet;
    } // method

    std::shared_ptr<libMA::Seeds> seed_one( std::string& sQuery, bool bRectangular, std::shared_ptr<libMA::Pack> pPack )
    {
        const char* sSeq = sQuery.c_str( );
        const int iSize = (int)sQuery.size( );
        return seed_one( sSeq, iSize, bRectangular, pPack, &mm_filter_none, nullptr );
    } // method

    std::shared_ptr<libMA::Seeds> seed_one( const char* sSeq, const int iSize, bool bRectangular,
                                            std::shared_ptr<libMA::Pack> pPack )
    {
        return seed_one( sSeq, iSize, bRectangular, pPack, &mm_filter_none, nullptr );
    } // method

#define VOID_MM UINT64_MAX
    static unsigned char seq_nt4_table[ 256 ];

    static inline uint64_t hash64( uint64_t key, uint64_t mask )
    {
        key = ( ~key + ( key << 21 ) ) & mask; // key = (key << 21) - key - 1;
        key = key ^ key >> 24;
        key = ( ( key + ( key << 3 ) ) + ( key << 8 ) ) & mask; // key * 265
        key = key ^ key >> 14;
        key = ( ( key + ( key << 2 ) ) + ( key << 4 ) ) & mask; // key * 21
        key = key ^ key >> 28;
        key = ( key + ( key << 31 ) ) & mask;
        return key;
    }
    /**
     * Adapted from the Minimap code
     * Find symmetric (w,k)-minimizers on a DNA sequence
     * Implements the algorithm proposed in:
     * "Minimap and miniasm: fast mapping and de novo assembly for noisy long sequences"
     *
     * @param km     thread-local memory pool; using NULL falls back to malloc()
     * @param str    DNA sequence
     * @param len    length of $str
     * @param w      find a minimizer for every $w consecutive k-mers (windows size)
     * @param k      k-mer size
     * @param rid    reference ID; will be copied to the output $p array
     * @param is_hpc homopolymer-compressed or not
     * @param p      minimizers
     *               p->a[i].x = kMer<<8 | kmerSpan
     *               p->a[i].y = rid<<32 | lastPos<<1 | strand
     *               where lastPos is the position of the last base of the i-th minimizer,
     *               and strand indicates whether the minimizer comes from the top or the bottom strand.
     *               Callers may want to set "p->n = 0"; otherwise results are appended to p
     */
    static std::vector<mm128_t> mm_sketch( const char* str, int len, int w, int k, uint32_t rid )
    {
        std::vector<mm128_t> p;
        uint64_t shift1 = 2 * ( k - 1 ), mask = ( 1ULL << 2 * k ) - 1, kmer[ 2 ] = {0, 0};
        int i, j, l, buf_pos, min_pos, kmer_span = 0;
        mm128_t buf[ 256 ]; // Circular Buffer of size w. stores kmer[ 0 ] and kmer[ 1 ] of each interation
        mm128_t min = {VOID_MM, VOID_MM}; // minimizer
        /* tiny_queue_t tq; */
        // std::cout << toBinString( mask ) << std::endl;
        assert( len > 0 && ( w > 0 && w < 256 ) &&
                ( k > 0 && k <= 28 ) ); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
        memset( buf, 0xff, w * 16 );
        /* memset( &tq, 0, sizeof( tiny_queue_t ) ); */
        // kv_resize( mm128_t, km, *p, p->n + len / w );
        p.reserve( p.size( ) + len / w );

        for( i = l = buf_pos = min_pos = 0; i < len; ++i )
        {
            int c = seq_nt4_table[ (uint8_t)str[ i ] ];
            mm128_t info = {VOID_MM, VOID_MM};
            if( c < 4 )
            { // not an ambiguous base
                int z;

                // if (is_hpc)
                // {
                //     int skip_len = 1;
                //     if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c)
                //     {
                //         for (skip_len = 2; i + skip_len < len; ++skip_len)
                //             if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
                //                 break;
                //         i += skip_len - 1; // put $i at the end of the current homopolymer run
                //     }
                //     tq_push(&tq, skip_len);
                //     kmer_span += skip_len;
                //     if (tq.count > k)
                //         kmer_span -= tq_shift(&tq);
                // }
                // else
                kmer_span = l + 1 < k ? l + 1 : k;
                // std::cout << "X" << kmer_span << " ";
                // I think her is a difference to the paper, where first comes hashing and than minimum
                kmer[ 0 ] = ( kmer[ 0 ] << 2 | c ) & mask; // forward k-mer
                kmer[ 1 ] = ( kmer[ 1 ] >> 2 ) | ( 3ULL ^ c ) << shift1; // reverse k-mer
                // std::cout << "0:" << toBinString<uint64_t, K_MER_SIZE * 2>( kmer[ 0 ] ) << std::endl;
                // std::cout << "1:" << toBinString<uint64_t, K_MER_SIZE * 2>( kmer[ 1 ] ) << std::endl;

                // Same k-mer on forward strand and reverse strand?
                if( kmer[ 0 ] == kmer[ 1 ] )
                    continue; // skip "symmetric k-mers" as we don't know it strand
                z = kmer[ 0 ] < kmer[ 1 ] ? 0 : 1; // strand
                // std::cout << "STRAND:" << z << std::endl;
                ++l;
                if( l >= k && kmer_span < 256 )
                {
                    // Hashing is for performance improvements merely
                    // See paper
                    info.x = hash64( kmer[ z ], mask ) << 8 | kmer_span;
                    info.y = (uint64_t)rid << 32 | (uint32_t)i << 1 | z;
                } // if
            } // if
            else
                l = 0, /* tq.count = tq.front = 0, */ kmer_span = 0;
            buf[ buf_pos ] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below

            if( l == w + k - 1 && min.x != VOID_MM ) // special case of the first window
            {
                // identical k-mers are not stored yet
                // std::cout << "FIRST WINDOW - l:" << l << std::endl;

                // Inspect the circular buffer except for buf_pos itself
                for( j = buf_pos + 1; j < w; ++j )
                    if( min.x == buf[ j ].x && buf[ j ].y != min.y )
                        // kv_push( mm128_t, km, *p, buf[ j ] );
                        p.push_back( buf[ j ] ); // store minimizer
                for( j = 0; j < buf_pos; ++j )
                    if( min.x == buf[ j ].x && buf[ j ].y != min.y )
                        // kv_push( mm128_t, km, *p, buf[ j ] );
                        p.push_back( buf[ j ] );
            } // if

            if( info.x <= min.x ) // a new minimum
            {
                // std::cout << "Case 2 - l:" << l << std::endl;

                if( l >= w + k && min.x != VOID_MM ) // record the old min
                    p.push_back( min ); // store minimizer
                min = info; // store new MM itself
                min_pos = buf_pos; // store position of new MM in buffer
            } // if

            else if( buf_pos == min_pos )
            { // old min has moved outside the window
                // std::cout << "Case 3 - l:" << l << std::endl;
                if( l >= w + k - 1 && min.x != VOID_MM )
                    // kv_push( mm128_t, km, *p, min );
                    p.push_back( min ); // store minimizer

                // the two loops are necessary when there are identical k-mers
                min.x = VOID_MM;
                // Search for minimal position in buffer
                for( j = buf_pos + 1; j < w; ++j )
                    if( min.x >= buf[ j ].x )
                        min = buf[ j ], min_pos = j; // >= is important s.t. min is always the closest k-mer

                for( j = 0; j <= buf_pos; ++j )
                    if( min.x >= buf[ j ].x )
                        min = buf[ j ], min_pos = j;

                if( l >= w + k - 1 && min.x != VOID_MM )
                { // write identical k-mers
                    for( j = buf_pos + 1; j < w; ++j ) // these two loops make sure the output is sorted
                        if( min.x == buf[ j ].x && min.y != buf[ j ].y )
                            // kv_push( mm128_t, km, *p, buf[ j ] );
                            p.push_back( buf[ j ] ); // store minimizer
                    for( j = 0; j <= buf_pos; ++j )
                        if( min.x == buf[ j ].x && min.y != buf[ j ].y )
                            // kv_push( mm128_t, km, *p, buf[ j ] );
                            p.push_back( buf[ j ] ); // store minimizer
                }
            }

            if( ++buf_pos == w )
                buf_pos = 0;
        }
        if( min.x != VOID_MM )
            // kv_push( mm128_t, km, *p, min );
            p.push_back( min ); // store minimizer
        return p;
    } // method

    static uint64_t _getHash( mm128_t kMerRepresentation )
    {
        return kMerRepresentation.x >> 8;
    } // method

    static std::vector<uint64_t> _getHash( std::vector<mm128_t> vKmers )
    {
        std::vector<uint64_t> vRet;
        vRet.reserve( vKmers.size( ) );
        for( auto kMer : vKmers )
            vRet.push_back( _getHash( kMer ) );
        return vRet;
    } // method


    static std::vector<uint64_t> _getHash( const char* sSeq, const int iSize, int uiK, int uiW )
    {
        return _getHash( mm_sketch( sSeq, iSize, uiW, uiK, 0 ) );
    } // method

    std::vector<uint64_t> getHash( const char* sSeq, const int iSize, int uiK, int uiW )
    {
        return _getHash( sSeq, iSize, xOptions.w, xOptions.k );
    } // method

    std::vector<std::shared_ptr<libMA::Seeds>>
    seed( std::vector<std::string> vQueries, bool bRectangular, std::shared_ptr<libMA::Pack> pPack )
    {
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;

        for( auto& sQuery : vQueries )
            vRet.push_back( seed_one( sQuery, bRectangular, pPack ) );

        return vRet;
    } // method

    void setMaxOcc( size_t uiNewMaxOcc )
    {
        xMapOpt.mid_occ = uiNewMaxOcc;
    }
}; // class

#endif
} // namespace minimizer


#ifdef WITH_PYTHON
void exportMinimizerIndex( libMS::SubmoduleOrganizer& xOrganizer );
#endif
