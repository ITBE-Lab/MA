
#include "container/container.h"
#include "container/nucSeq.h"
#include "container/pack.h"
#include "container/segment.h"
#include "util/parameter.h"
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
/**
 * @brief cpp wrapper for the minimap index
 * @details
 * at the moment multipart indices are not supported
 */
class Index : public libMA::Container
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

    // clang-format off
    Index( const ParameterSetManager& rParameters ) //
        : xOptions{ /* .k = */ (short)rParameters.getSelected( )->xMinimizerK->get( ),
                    /* .w = */ (short)rParameters.getSelected( )->xMinimizerW->get( ),
#if HIDE_MMI_PARAMETERS == 0
                    /* .flag = */ rParameters.getSelected( )->xMinimizerFlag->get( ),
                    /* .bucket_bits = */ rParameters.getSelected( )->xMinimizerBucketBits->get( ),
                    /* .mini_batch_size = */ rParameters.getSelected( )->xMinimizerMiniBatchSize->get( ),
                    /* .batch_size = */ rParameters.getSelected( )->xMinimizerBatchSize->get( )}
#else
                    /* .flag = */ 0,
                    /* .bucket_bits = */ 14,
                    /* .mini_batch_size = */ 50000000,
                    /* .batch_size = */ 4000000000ULL}
#endif
    {} // default constructor
    // clang-format on

  public:
    /**
     * @brief Opens in Index if index filename is given or creates index if fasta(q) file is given
     */
    Index( const ParameterSetManager& rParameters, std::string sIndexName ) : Index( rParameters )
    {
        IndexReader xIndexReader( sIndexName, &xOptions, NULL );
        if( xIndexReader.idx_rdr == 0 )
            throw std::runtime_error( "failed to open file" + sIndexName ); //@todo this does not seem to work...
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
        : Index( rParameters )
    {
        std::vector<const char*> seq;
        std::vector<const char*> name;
        for( auto& sContig : vContigs )
            seq.push_back( sContig.data( ) );
        for( auto& sName : vContigsNames )
            name.push_back( sName.data( ) );

        pData = mm_idx_str( xOptions.w, xOptions.k, xOptions.flag & MM_I_HPC, xOptions.bucket_bits,
                            (int)vContigs.size( ), &seq[ 0 ], &name[ 0 ] );
        initOptions( );
        mm_idx_stat( pData );
        std::cout << "max Occ:" << xMapOpt.mid_occ << std::endl;
    } // constructor

    ~Index( )
    {
        mm_idx_destroy( pData );
    } // destructor

    void dump( std::string sIndexName )
    {
        if( mm_idx_dump_c_str( sIndexName.c_str( ), pData ) != 0 )
            throw std::runtime_error( "failed to open file" + sIndexName );
    } // method


    void setMidOcc(int32_t iNewVal)
    {
        xMapOpt.mid_occ = iNewVal;
    }

    std::shared_ptr<libMA::Seeds> seed_one( std::string& sQuery, std::shared_ptr<libMA::Pack> pPack )
    {
        auto pRet = std::make_shared<libMA::Seeds>( );
        const char* sSeq = sQuery.c_str( );
        const int iSize = (int)sQuery.size( );
        int64_t n_a = 0;
        mm_tbuf_t* tbuf = mm_tbuf_init( );
        // const char* seqs = (const char*)&pQuery->pxSequenceRef;
        mm128_t* a = collect_seeds( pData, 1, &iSize, &sSeq, tbuf, &xMapOpt, 0, &n_a );
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
                        ( (int32_t)a[ uiI ].x ) + 1 - uiQSpan + pPack->startOfSequenceWithId( a[ uiI ].x << 1 >> 33 )
#if AMBIGUITY == ( 1 )
                            ,
                        1
#endif
                    );
                } // if
                else
                {
                    uint64_t uiQSpan = ( int32_t )( a[ uiI ].y >> 32 & 0xff );
                    pRet->emplace_back(
                        /* minimap uses 1-based index i use 0-based indices...
                            a[ uiI ].y = last base of minimizer */
                        sQuery.length( ) - 1 - ( (int32_t)a[ uiI ].y ),
                        uiQSpan,
                        /* reference position is encoded in the lower 32 bits of x */
                        pPack->uiPositionToReverseStrand( (int32_t)a[ uiI ].x +
                                                          pPack->startOfSequenceWithId( a[ uiI ].x << 1 >> 33 ) )
#if AMBIGUITY == ( 1 )
                            ,
                        1
#endif
                    );
                } // else
            } // for
            kfree( tbuf->km, a );
        } // if
        else if( n_a > 0 )
            throw std::runtime_error( "minimizer vector is empty" );
        mm_tbuf_destroy( tbuf );
        return pRet;
    } // method

    std::vector<std::shared_ptr<libMA::Seeds>>
    seed( std::vector<std::string> vQueries, std::shared_ptr<libMA::Pack> pPack )
    {
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;

        for( auto& sQuery : vQueries )
            vRet.push_back( seed_one( sQuery, pPack ) );

        return vRet;
    } // method
}; // class

#endif
} // namespace minimizer


#ifdef WITH_PYTHON
void exportMinimizerIndex( py::module& rxPyModuleId );
#endif
