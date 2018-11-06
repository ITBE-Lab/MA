/**
 * @file stripOfConsideration.h
 * @brief Implements a Strip of Consideration.
 * @author Markus Schmidt
 */
#ifndef STRIP_OF_CONSIDERATION_H
#define STRIP_OF_CONSIDERATION_H

#include "container/segment.h"
#include "container/soc.h"
#include "module/module.h"
#include <cmath>

namespace libMA
{
/**
 * @brief Used to quickly find areas with high density of @ref Seed "seeds".
 * @ingroup module
 */
class StripOfConsideration : public Module<SoCPriorityQueue, false, SegmentVector, NucSeq, Pack, FMIndex>
{
  public:
    /// @brief Maximum ambiguity for a seed to be considered.
    size_t uiMaxAmbiguity = defaults::uiMaxAmbiguity;
    /// @brief Minimum seed length.
    size_t uiMinLen = defaults::uiMinLen;
    /**
     * @brief Minimal SoC score.
     * @details
     * Must be [0,-inf]!
     * The minimal score that should be allowed during SoC collection.
     * Is interpreted relative to the query length.
     */
    float fScoreMinimum = defaults::fSoCScoreMinimum;
    /**
     * @brief If the best SoC has seeds of accumulative length smaller than this, abort.
     * @details
     * Is multiplied by query length.
     * 0 = never abort.
     */
    const float fGiveUp;
    size_t uiCurrHarmScoreMin = defaults::uiCurrHarmScoreMin;

    /**
     * @brief disable fGiveUp and fRelMinSeedSizeAmount if genome is too short
     */
    nucSeqIndex uiMinGenomeSize = defaults::uiGenomeSizeDisable;

    size_t uiSoCWidth = defaults::uiSoCWidth;

    /**
     * @brief skip seeds with too much ambiguity
     * @details
     * True: skip all seeds with to much ambiguity
     * False: use max_hits instances of the seeds with more ambiguity
     */
    bool bSkipLongBWTIntervals = defaults::bSkipLongBWTIntervals;

    inline static nucSeqIndex getPositionForBucketing( nucSeqIndex uiQueryLength, const Seed& xS )
    {
        return xS.start_ref( ) + ( uiQueryLength - xS.start( ) );
    } // function

    inline nucSeqIndex getStripSize( nucSeqIndex uiQueryLength, int iMatch, int iExtend, int iGap ) const
    {
        if( uiSoCWidth != 0 )
            return uiSoCWidth;

        return ( iMatch * uiQueryLength - iGap ) / iExtend - ( int64_t )( fScoreMinimum * uiQueryLength );
    } // function

    template <class FUNCTOR>
    inline void forEachNonBridgingSeed( SegmentVector& rSegmentVector, FMIndex& rxFM_index, Pack& rxRefSequence,
                                        FUNCTOR&& fDo, // std::function<void(const Seed &rxS)> fDo,
                                        const nucSeqIndex uiQLen, nucSeqIndex addSize = 0 )
    {
        rSegmentVector.forEachSeed( rxFM_index, uiQLen, uiMaxAmbiguity, bSkipLongBWTIntervals,
                                    [&rxRefSequence, &rxFM_index, &fDo, &addSize]( Seed&& xS ) {
                                        // check if the match is bridging the forward/reverse strand
                                        // or bridging between two chromosomes
                                        if( !rxRefSequence.bridgingSubsection(
                                                // prevent negative index
                                                xS.start_ref( ) > addSize ? xS.start_ref( ) - addSize : 0, // from
                                                // prevent index larger than reference
                                                xS.end_ref( ) + addSize <= rxFM_index.getRefSeqLength( )
                                                    ? xS.size( ) - 1 + addSize
                                                    : rxFM_index.getRefSeqLength( ) - xS.start_ref( ) - 1 // to
                                                ) )
                                        {
                                            // if non-bridging use this seed
                                            fDo( xS );
                                        } // if
                                        // returning true since we want to continue extracting seeds
                                        return true;
                                    } // lambda
        );
    } // method


    inline void emplaceAllNonBridgingSeed( SegmentVector& rSegmentVector, FMIndex& rxFM_index, Pack& rxRefSequence,
                                           Seeds& rvSeedVector, const nucSeqIndex uiQLen )
    {
        rSegmentVector.emplaceAllEachSeeds(
            rxFM_index, uiQLen, uiMaxAmbiguity, uiMinLen, rvSeedVector,
            [&rxRefSequence, &rxFM_index, &rvSeedVector, &uiQLen]( ) {
                /*
                 * @note this bridging check is not required since we check weather a SoC
                 * is brigding in general.
                 * If any of the seeds within a SoC are bridging then the SoC is bridging.
                 * Could turn this into a debug assertion...
                 */
                rvSeedVector.back( ).uiDelta = getPositionForBucketing( uiQLen, rvSeedVector.back( ) );

                size_t uiContig = rxRefSequence.uiSequenceIdForPosition( rvSeedVector.back( ).start_ref( ) );
                // move the delta position to the right by the contig index
                // this way we ensure seeds are ordered by contigs first
                // and then their delta positions within the contig.
                rvSeedVector.back( ).uiDelta += ( uiQLen + 10 ) * uiContig;
                // returning true since we want to continue extracting seeds
                return true;
            } // lambda
        );
    } // method

  public:
    StripOfConsideration( ) : fGiveUp( defaults::fGiveUp )
    {} // default constructor

    virtual std::shared_ptr<SoCPriorityQueue> EXPORTED execute( std::shared_ptr<SegmentVector> pSegments,
                                                                std::shared_ptr<NucSeq>
                                                                    pQuerySeq,
                                                                std::shared_ptr<Pack>
                                                                    pRefSeq,
                                                                std::shared_ptr<FMIndex>
                                                                    pFM_index );
}; // class


} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the bucketing @ref Module "module" to python.
 * @ingroup export
 */
void exportStripOfConsideration( );
#endif

#endif