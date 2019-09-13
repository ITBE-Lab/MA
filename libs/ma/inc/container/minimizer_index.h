#include "minimap.h"
#include "util/parameter.h"

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
            mm_idx_reader_close( idx_rdr );
        } // destructor
    }; // class
  public:
    /**
     * @brief Opens in Index if index filename is given or creates index if fasta(q) file is given
     */
    Index( std::shared_ptr<ParameterSetManager> pParameters, std::string sIndexName )
        : xOptions{.k = pParameters->getSelected( )->xMinimizerK->get( ),
                   .w = pParameters->getSelected( )->xMinimizerW->get( ),
                   .flag = pParameters->getSelected( )->xMinimizerFlag->get( ),
                   .bucket_bits = pParameters->getSelected( )->xMinimizerBucketBits->get( ),
                   .mini_batch_size = pParameters->getSelected( )->xMinimizerMiniBatchSize->get( ),
                   .batch_size = pParameters->getSelected( )->xMinimizerBatchSize->get( )}
    {
        IndexReader xIndexReader( sIndexName, &xOptions, NULL );
        if( xIndexReader.idx_rdr == 0 )
            throw std::runtime_error( "failed to open file" + sIndexName );
        pData = mm_idx_reader_read( xIndexReader.idx_rdr, pParameters->getNumThreads( ) );
        if( ( xOptions.flag & MM_F_CIGAR ) && ( pData->flag & MM_I_NO_SEQ ) )
            throw std::runtime_error( "the prebuilt index doesn't contain sequences" );
#if DEBUG_LEVEL > 0
        mm_idx_stat( pData );
#endif
        // make sure that there are no multipart indices
        assert( mm_idx_reader_read( xIndexReader.idx_rdr, pParameters->getNumThreads( ) ) == 0 );
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

}; // class

} // namespace minimizer

#ifdef WITH_PYTHON
void exportMinimizerIndex( py::module& rxPyModuleId );
#endif