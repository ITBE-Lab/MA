/**
 * @file needlemanWunsch.h
 * @brief Implements NMW
 * @author Markus Schmidt
 */
#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H

#define OLD_KSW ( 0 )

#include "container/alignment.h"
#include "kswcpp.h"
#include "module/module.h"


namespace libMA
{
// small wrapper that takes care of deallocation
class Wrapper_ksw_extz_t
{
  public:
#if OLD_KSW == 1
    ksw_extz_t* ez;
#else
    kswcpp_extz_t* ez;
#endif

    Wrapper_ksw_extz_t( )
    {
#if OLD_KSW == 1
        ez = new ksw_extz_t{}; // {} forces zero initialization
#else
        ez = new kswcpp_extz_t{}; // {} forces zero initialization
#endif

    } // default constructor

    ~Wrapper_ksw_extz_t( )
    {
        free( ez->cigar ); // malloced in c code
        delete ez; // allocated by new in cpp code
    } // default constructor
}; // class

/**
 * @brief implements NMW
 * @details
 * Returns a finished alignment if given a sound selection of seeds.
 * @ingroup module
 */
class NeedlemanWunsch : public Module<ContainerVector<std::shared_ptr<Alignment>>, false,
                                      ContainerVector<std::shared_ptr<Seeds>>, NucSeq, Pack>
{
    std::vector<std::vector<std::vector<int>>> s;
    std::vector<std::vector<std::vector<char>>> dir;

    void EXPORTED dynPrg( const std::shared_ptr<NucSeq> pQuery, const std::shared_ptr<NucSeq> pRef,
                          const nucSeqIndex fromQuery, const nucSeqIndex toQuery, const nucSeqIndex fromRef,
                          const nucSeqIndex toRef,
                          std::shared_ptr<Alignment> pAlignment, // in & output
                          AlignedMemoryManager& rMemoryManager, const bool bLocalBeginning, const bool bLocalEnd );

    void ksw_dual_ext( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<NucSeq> pRef, nucSeqIndex fromQuery,
                       nucSeqIndex toQuery, nucSeqIndex fromRef, nucSeqIndex toRef,
                       std::shared_ptr<Alignment> pAlignment, AlignedMemoryManager& rMemoryManager );

    void ksw( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<NucSeq> pRef, nucSeqIndex fromQuery, nucSeqIndex toQuery,
              nucSeqIndex fromRef, nucSeqIndex toRef, std::shared_ptr<Alignment> pAlignment,
              AlignedMemoryManager& rMemoryManager );

    const KswCppParam<5> xKswParameters;
    const nucSeqIndex uiMaxGapArea;
    const nucSeqIndex uiPadding;
    const nucSeqIndex uiMissMatch;
    const size_t uiZDrop;
    /*
     * @details
     * Minimal bandwith when filling the gap between seeds.
     * This bandwidth gets increased if the seeds are not on a diagonal.
     */
    const int iMinBandwidthGapFilling;
    // bandwidth when extending the edge of an alignment
    const int iBandwidthDPExtension;

  public:
    bool bLocal = false;

    NeedlemanWunsch( const ParameterSetManager& rParameters )
        : xKswParameters( rParameters.getSelected( )->xMatch->get( ),
                          rParameters.getSelected( )->xMisMatch->get( ),
                          rParameters.getSelected( )->xGap->get( ),
                          rParameters.getSelected( )->xExtend->get( ),
                          rParameters.getSelected( )->xGap2->get( ),
                          rParameters.getSelected( )->xExtend2->get( ) ),
          uiMaxGapArea( rParameters.getSelected( )->xMaxGapArea->get( ) ),
          uiPadding( rParameters.getSelected( )->xPadding->get( ) ),
          uiMissMatch( rParameters.getSelected( )->xMisMatch->get( ) ),
          uiZDrop( rParameters.getSelected( )->xZDrop->get( ) ),
          iMinBandwidthGapFilling( rParameters.getSelected( )->xMinBandwidthGapFilling->get( ) ),
          iBandwidthDPExtension( rParameters.getSelected( )->xBandwidthDPExtension->get( ) ){}; // default constructor


    std::shared_ptr<Alignment> EXPORTED execute_one( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuery,
                                                     std::shared_ptr<Pack> pRefPack,
                                                     AlignedMemoryManager& rMemoryManager );

    // overload
    virtual std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
    execute( std::shared_ptr<ContainerVector<std::shared_ptr<Seeds>>> pSeedSets, std::shared_ptr<NucSeq> pQuery,
             std::shared_ptr<Pack> pRefPack )
    {
        AlignedMemoryManager xMemoryManager;
        auto pRet = std::make_shared<ContainerVector<std::shared_ptr<Alignment>>>( );
        for( auto pSeeds : *pSeedSets )
        {
            // if we have a seed set with an inversion we need to generate two alignments...
            //if( !pSeeds->mainStrandIsForward( ) )
            //    pSeeds->mirror( pRefPack->uiStartOfReverseStrand( ), pQuery->length( ) );
            pRet->push_back( execute_one( pSeeds, pQuery, pRefPack, xMemoryManager ) );

            // auto pRevSeeds = pSeeds->splitOnStrands( pRefPack->uiStartOfReverseStrand( ), pQuery->length() );
            // if( !pRevSeeds->empty( ) )
            //    pRet->push_back( execute_one( pRevSeeds, pQuery, pRefPack, xMemoryManager ) );
            // if( !pSeeds->empty( ) )
            //    pRet->push_back( execute_one( pSeeds, pQuery, pRefPack, xMemoryManager ) );
        } // for
        // we need to move the best alignment to the first spot
        std::sort( pRet->begin( ), pRet->end( ),
                   []( std::shared_ptr<Alignment>& pA, std::shared_ptr<Alignment>& pB ) { return pA->larger( pB ); } );
        return pRet;
    } // function

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief Exposes the Alignment container to boost python.
 * @ingroup export
 */
#ifdef WITH_BOOST
void exportNeedlemanWunsch( );
#else
void exportNeedlemanWunsch( py::module& rxPyModuleId );
#endif
#endif

#endif