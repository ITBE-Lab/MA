/**
 * @file pairedReads.h
 * @brief Picks alignment pairs for alignming paired reads.
 * @author Markus Schmidt
 */
#ifndef PAIRED_READS_H
#define PAIRED_READS_H

#include "ma/container/alignment.h"
#include "ms/module/module.h"
#include <cmath>

namespace libMA
{
/**
 * @brief Picks paired alignments for paired reads
 * @ingroup module
 * @details
 * Given two vector of alignments this module picks one alignment of each vector,
 * so that their combined score and penalty for the insert size is maximal.
 */
class PairedReads
    : public libMS::Module<libMS::ContainerVector<std::shared_ptr<Alignment>>, // out
                           false, // not volatile
                           NucSeq, NucSeq, // Query 1 and 2
                           libMS::ContainerVector<std::shared_ptr<Alignment>>, // alignments for first paired read
                           libMS::ContainerVector<std::shared_ptr<Alignment>>, // alignments for second paired read
                           Pack>
{
  public:
    /**
     * @brief Penalty for unpaired reads.
     * @details
     * the score penalty that is substracted from the SW scores of the alignment of two *
     * unpaired reads.
     * If the penalty is smaller than the penalty added by the insert size
     * between the paired seed we make an unpaired alignment.
     * The default value is taken from BWA-MEM as we use the same formula for pairing.
     */
    double u;
    ///@brief the mean of the insert size
    size_t mean;
    ///@brief the standard deviation of the insert size
    double std;

    PairedReads( const ParameterSetManager& rParameters )
        : u( rParameters.getSelected( )->xPairedBonus->get( ) ),
          mean( (size_t)rParameters.getSelected( )->xMeanPairedReadDistance->get( ) ),
          std( rParameters.getSelected( )->xStdPairedReadDistance->get( ) )
    {} // constructor

    /**
     * @brief The probability for a insert size >= d.
     * @details
     * Can be set to use normal or uniform distribution.
     * Note: both distributions are set using mean and standard deviation.
     */
    double DLL_PORT(MA) p( nucSeqIndex d ) const;

    std::shared_ptr<libMS::ContainerVector<std::shared_ptr<Alignment>>>
        DLL_PORT(MA) execute( std::shared_ptr<NucSeq> pQ1, std::shared_ptr<NucSeq> pQ2,
                          std::shared_ptr<libMS::ContainerVector<std::shared_ptr<Alignment>>> pA,
                          std::shared_ptr<libMS::ContainerVector<std::shared_ptr<Alignment>>> pB,
                          std::shared_ptr<Pack> pPack );
}; // class
} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the PairedReads and PairedAlignment @ref libMA::Module "module" to python.
 * @ingroup export
 */
#ifdef WITH_BOOST
void exportPairedReads( );
#else
void exportPairedReads( libMS::SubmoduleOrganizer& xOrganizer );
#endif
#endif


#endif