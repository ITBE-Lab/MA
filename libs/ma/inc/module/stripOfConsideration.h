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
#include "module/harmonization.h"
#include <cmath>

namespace libMA
{
class StripOfConsiderationBase
{
  public:
    /// @brief Maximum ambiguity for a seed to be considered.
    const size_t uiMaxAmbiguity;
    /// @brief Minimum seed length.
    const size_t uiMinLen;
    /**
     * @brief If the best SoC has seeds of accumulative length smaller than this, abort.
     * @details
     * Is multiplied by query length.
     * 0 = never abort.
     */
    const double fGiveUp;
    const size_t uiCurrHarmScoreMin;

    /**
     * @brief disable fGiveUp and fRelMinSeedSizeAmount if genome is too short
     */
    const nucSeqIndex uiMinGenomeSize;

    const size_t uiSoCWidth;

    /**
     * @brief skip seeds with too much ambiguity
     * @details
     * True: skip all seeds with to much ambiguity
     * False: use max_hits instances of the seeds with more ambiguity
     */
    const bool bSkipLongBWTIntervals;

    inline static nucSeqIndex getPositionForBucketing( nucSeqIndex uiQueryLength, const Seed& xS )
    {
        return xS.start_ref( ) + ( uiQueryLength - xS.start( ) );
    } // function

    inline nucSeqIndex getStripSize( nucSeqIndex uiQueryLength, int iMatch, int iExtend, int iGap ) const
    {
        if( uiSoCWidth != 0 )
            return uiSoCWidth;

        return ( iMatch * uiQueryLength - iGap ) / iExtend;
    } // function

    StripOfConsiderationBase( const ParameterSetManager& rParameters )
        : uiMaxAmbiguity( rParameters.getSelected( )->xMaximalSeedAmbiguity->get( ) ),
          uiMinLen( rParameters.getSelected( )->xMinSeedLength->get( ) ),
          fGiveUp( rParameters.getSelected( )->xHarmScoreMinRel->get( ) ),
          uiCurrHarmScoreMin( rParameters.getSelected( )->xHarmScoreMin->get( ) ),
          uiMinGenomeSize( rParameters.getSelected( )->xGenomeSizeDisable->get( ) ),
          uiSoCWidth( rParameters.getSelected( )->xSoCWidth->get( ) ),
          bSkipLongBWTIntervals( rParameters.getSelected( )->xSkipAmbiguousSeeds->get( ) )
    {} // default constructor
}; // class

/**
 * @brief Used to quickly find areas with high density of @ref Seed "seeds".
 * @ingroup module
 */
class StripOfConsideration : public Module<SoCPriorityQueue, false, SegmentVector, NucSeq, Pack, FMIndex>,
                             StripOfConsiderationBase
{
  public:
    template <class FUNCTOR>
    inline void forEachNonBridgingSeed( SegmentVector& rSegmentVector, FMIndex& rxFM_index, Pack& rxRefSequence,
                                        FUNCTOR&& fDo, // std::function<void(const Seed &rxS)> fDo,
                                        nucSeqIndex addSize = 0 )
    {
        rSegmentVector.forEachSeed( rxFM_index, uiMaxAmbiguity, bSkipLongBWTIntervals,
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
        rSegmentVector.emplaceAllEachSeeds( rxFM_index, uiMaxAmbiguity, uiMinLen, bSkipLongBWTIntervals, rvSeedVector,
                                            [&rxRefSequence, &rxFM_index, &rvSeedVector, &uiQLen]( ) {
        /*
         * @note this bridging check is not required since we check weather a SoC
         * is brigding in general.
         * If any of the seeds within a SoC are bridging then the SoC is bridging.
         * Could turn this into a debug assertion...
         */
#if 0 // enable / disable the bridging check for all seeds
                    constexpr const nucSeqIndex addSize = 0;
                    auto& rS = rvSeedVector.back();
                    // check if the match is bridging the forward/reverse strand 
                    // or bridging between two chromosomes
                    if ( rxRefSequence.bridgingSubsection(
                            //prevent negative index
                            rS.start_ref() > addSize ? rS.start_ref() - addSize : 0,//from
                            //prevent index larger than reference
                            rS.end_ref() + addSize <= rxFM_index.getRefSeqLength() ?
                                rS.size() - 1 + addSize :
                                rxFM_index.getRefSeqLength() - rS.start_ref() - 1//to
                            ) 
                        )
                    {
                        //if bridging remove this seed
                        rvSeedVector.pop_back();
                    }//if
#if DELTA_CACHE == ( 1 )
                    else
                    {
                        // set the delta cache
                        rS.uiDelta = getPositionForBucketing( uiQLen, rS );
                    }// else
#endif
#elif DELTA_CACHE == ( 1 )
                                                rvSeedVector.back( ).uiDelta =
                                                    getPositionForBucketing( uiQLen, rvSeedVector.back( ) );
#endif
#if CONTIG_ID_CACHE == ( 1 )
                                                rvSeedVector.back( ).uiContigId = rxRefSequence.uiSequenceIdForPosition(
                                                    rvSeedVector.back( ).start_ref( ) );
#endif
                                                // returning true since we want to continue extracting seeds
                                                return true;
                                            } // lambda
        );
    } // method

  public:
    StripOfConsideration( const ParameterSetManager& rParameters ) : StripOfConsiderationBase( rParameters )
    {} // default constructor

    virtual std::shared_ptr<SoCPriorityQueue> EXPORTED execute( std::shared_ptr<SegmentVector> pSegments,
                                                                std::shared_ptr<NucSeq>
                                                                    pQuerySeq,
                                                                std::shared_ptr<Pack>
                                                                    pRefSeq,
                                                                std::shared_ptr<FMIndex>
                                                                    pFM_index );
}; // class

class MinimizerStripOfConsideration : public Module<SoCPriorityQueue, false, Seeds, NucSeq, Pack>,
                                      StripOfConsiderationBase
{
  public:
    SeedLumping<true> xLumper;
    SeedExtender xExtender;

    MinimizerStripOfConsideration( const ParameterSetManager& rParameters )
        : StripOfConsiderationBase( rParameters ), xLumper( rParameters ), xExtender( rParameters )
    {} // constructor

    virtual std::shared_ptr<SoCPriorityQueue> EXPORTED execute( std::shared_ptr<Seeds> pSeeds,
                                                                std::shared_ptr<NucSeq> pQuerySeq,
                                                                std::shared_ptr<Pack> pPack );
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief export the bucketing @ref libMA::Module "module" to python.
 * @ingroup export
 */
#ifdef WITH_BOOST
void exportStripOfConsideration( );
#else
void exportStripOfConsideration( py::module& rxPyModuleId );
#endif
#endif

#endif