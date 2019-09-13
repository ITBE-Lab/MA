#include "container/nucSeq.h"
#include "container/segment.h"
#include "util/parameter.h"
#ifndef __cplusplus
#define __cplusplus TRUE
#endif
#include "minimap.h"

#pragma once

namespace minimizer
{
/**
 * @brief cpp wrapper for the minimap index
 * @details
 * at the moment multipart indices are not supported
 */
class Index
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

        opt->mask_level = 0.5f;
        opt->pri_ratio = 0.8f;
        opt->best_n = 5;

        opt->max_join_long = 20000;
        opt->max_join_short = 2000;
        opt->min_join_flank_sc = 1000;
        opt->min_join_flank_ratio = 0.5f;

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

  public:
    /**
     * @brief Opens in Index if index filename is given or creates index if fasta(q) file is given
     */
    Index( const ParameterSetManager& rParameters, std::string sIndexName )
        : xOptions{.k = rParameters.getSelected( )->xMinimizerK->get( ),
                   .w = rParameters.getSelected( )->xMinimizerW->get( ),
                   .flag = rParameters.getSelected( )->xMinimizerFlag->get( ),
                   .bucket_bits = rParameters.getSelected( )->xMinimizerBucketBits->get( ),
                   .mini_batch_size = rParameters.getSelected( )->xMinimizerMiniBatchSize->get( ),
                   .batch_size = rParameters.getSelected( )->xMinimizerBatchSize->get( )}
    {
        IndexReader xIndexReader( sIndexName, &xOptions, NULL );
        if( xIndexReader.idx_rdr == 0 )
            throw std::runtime_error( "failed to open file" + sIndexName );
        pData = mm_idx_reader_read( xIndexReader.idx_rdr, rParameters.getNumThreads( ) );
        if( ( xOptions.flag & MM_F_CIGAR ) && ( pData->flag & MM_I_NO_SEQ ) )
            throw std::runtime_error( "the prebuilt index doesn't contain sequences" );
#if DEBUG_LEVEL > 0
        mm_idx_stat( pData );
#endif
        mm_mapopt_init( &xMapOpt );
        // make sure that there are no multipart indices
        assert( mm_idx_reader_read( xIndexReader.idx_rdr, rParameters.getNumThreads( ) ) == 0 );
    } // constructor

    Index( const ParameterSetManager& rParameters, std::vector<std::shared_ptr<libMA::NucSeq>> vContigs,
           std::vector<std::string> vContigsNames )
        : xOptions{.k = rParameters.getSelected( )->xMinimizerK->get( ),
                   .w = rParameters.getSelected( )->xMinimizerW->get( ),
                   .flag = rParameters.getSelected( )->xMinimizerFlag->get( ),
                   .bucket_bits = rParameters.getSelected( )->xMinimizerBucketBits->get( ),
                   .mini_batch_size = rParameters.getSelected( )->xMinimizerMiniBatchSize->get( ),
                   .batch_size = rParameters.getSelected( )->xMinimizerBatchSize->get( )}
    {
        std::vector<const char*> seq;
        std::vector<const char*> name;
        for( auto pContig : vContigs )
            seq.push_back( (const char*)pContig->pxSequenceRef );
        for( auto sName : vContigsNames )
            name.push_back( sName.c_str( ) );

        pData = mm_idx_str( xOptions.w, xOptions.k, xOptions.flag & MM_I_HPC, xOptions.bucket_bits, vContigs.size( ),
                            &seq[ 0 ], &name[ 0 ] );
    } // constructor

    ~Index( )
    {
        mm_idx_destroy( pData );
    } // destructor

    void dump( std::string sIndexName )
    {
        auto* fp_out = fopen( sIndexName.c_str( ), "wb" );
        if( fp_out == 0 )
            throw std::runtime_error( "failed to open file" + sIndexName );
        mm_idx_dump( fp_out, pData );
        fclose( fp_out );
    } // method


    std::vector<std::shared_ptr<libMA::Seeds>> seed( std::vector<std::shared_ptr<libMA::NucSeq>> vQueries )
    {
        std::vector<std::shared_ptr<libMA::Seeds>> vRet;
        int64_t n_a;

        mm_tbuf_t* tbuf = mm_tbuf_init( );
        for( auto pQuery : vQueries )
        {
            const char* seqs = (const char*)&pQuery->pxSequenceRef;
            mm128_t* a = collect_seeds( pData, 1, (const int*)&pQuery->uiSize, &seqs, tbuf, &xMapOpt, 0, &n_a );
            if( a != NULL )
            {
                vRet.push_back( std::make_shared<libMA::Seeds>( ) );
                for( int64_t uiI = 0; uiI < n_a; uiI++ )
                {
                    vRet.back( )->emplace_back( (int32_t)a[ uiI ].x, (libMA::nucSeqIndex)xOptions.k,
                                                (int32_t)a[ uiI ].y, ( ( a[ uiI ].x >> 63 ) == 0 ? true : false ) );
                    if( uiI + 1 < n_a && a[ uiI ].x << 1 >> 33 != a[ uiI + 1 ].x << 1 >> 33 )
                        vRet.push_back( std::make_shared<libMA::Seeds>( ) );
                }
                free( a );
            }
            else if( n_a > 0 )
                throw std::runtime_error( "minimizer vector is empty" );
        }
        mm_tbuf_destroy( tbuf );

        return vRet;
    }
}; // class

} // namespace minimizer

#ifdef WITH_PYTHON
void exportMinimizerIndex( py::module& rxPyModuleId );
#endif