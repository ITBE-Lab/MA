/**
 * @file needlemanWunsch.h
 * @brief Implements NMW
 * @author Markus Schmidt
 */
#pragma once

#include "container/alignment.h"
#include "kswcpp.h"
#include "module/module.h"
#include "module/needlemanWunsch.h"
#include <limits>

namespace libMA
{
/**
 * @brief implements NMW
 * @details
 * Returns a finished alignment if given a sound selection of seeds.
 * @ingroup module
 */
class SmallInversions : public Module<ContainerVector<std::shared_ptr<Alignment>>, false,
                                      ContainerVector<std::shared_ptr<Alignment>>, NucSeq, Pack>
{

    const KswCppParam<5> xKswParameters;
    const nucSeqIndex uiMaxGapArea;
    const size_t uiZDrop;
    const size_t uiZDropInversion;
    const size_t uiBandwidth;
    const int iMatch;
    const int iMissMatch;
    const int iExtend;
    const int iGap;
    const bool bDisableHeuristics;
    const int iHarmScoreMin;
    const double dHarmScoreMinRel;

  public:
    bool bLocal = false;

    SmallInversions( const ParameterSetManager& rParameters )
        : xKswParameters( rParameters.getSelected( )->xMatch->get( ),
                          rParameters.getSelected( )->xMisMatch->get( ),
                          rParameters.getSelected( )->xGap->get( ),
                          rParameters.getSelected( )->xExtend->get( ),
                          rParameters.getSelected( )->xGap2->get( ),
                          rParameters.getSelected( )->xExtend2->get( ) ),
          uiMaxGapArea( rParameters.getSelected( )->xMaxGapArea->get( ) ),
          uiZDrop( rParameters.getSelected( )->xZDrop->get( ) ),
          uiZDropInversion( rParameters.getSelected( )->xZDropInversion->get( ) ),
          uiBandwidth( rParameters.getSelected( )->xBandwidthDPExtension->get( ) ),
          iMatch( rParameters.getSelected( )->xMatch->get( ) ),
          iMissMatch( rParameters.getSelected( )->xMisMatch->get( ) ),
          iExtend( rParameters.getSelected( )->xExtend->get( ) ),
          iGap( rParameters.getSelected( )->xGap->get( ) ),
          bDisableHeuristics( rParameters.getSelected( )->xDisableHeuristics->get( ) ),
          iHarmScoreMin( rParameters.getSelected( )->xHarmScoreMin->get( ) ),
          dHarmScoreMinRel( rParameters.getSelected( )->xHarmScoreMinRel->get( ) ){}; // default constructor

    inline void forAllDropPos( std::function<void( nucSeqIndex, nucSeqIndex )> fDo,
                               std::shared_ptr<Alignment> pAlignment )
    {
        nucSeqIndex uiMaxScorePosQ = pAlignment->uiBeginOnQuery;
        nucSeqIndex uiPosQ = uiMaxScorePosQ;
        nucSeqIndex uiMaxScorePosR = pAlignment->uiBeginOnRef;
        nucSeqIndex uiPosR = uiMaxScorePosR;
        int iMaxScore = std::numeric_limits<int>::min( );
        int iCurrScore = 0;
        for( std::pair<MatchType, nucSeqIndex> section : pAlignment->data )
        {
            switch( section.first )
            {
                case MatchType::seed:
                case MatchType::match:
                    iCurrScore += iMatch * section.second;
                    uiPosQ += section.second;
                    uiPosR += section.second;
                    break;
                case MatchType::missmatch:
                    iCurrScore += iMissMatch * section.second;
                    uiPosQ += section.second;
                    uiPosR += section.second;
                    break;
                case MatchType::insertion:
                    iCurrScore -= iGap + iExtend * section.second;
                    uiPosQ += section.second;
                    break;
                case MatchType::deletion:
                    iCurrScore -= iGap + iExtend * section.second;
                    uiPosR += section.second;
                    break;
                default:
                    std::cerr << "WARNING invalid cigar symbol" << std::endl;
                    break;
            } // switch
            if( iCurrScore >= iMaxScore )
            {
                iMaxScore = iCurrScore;
                uiMaxScorePosQ = uiPosQ;
                uiMaxScorePosR = uiPosR;
            } // if
            size_t uiDiff = std::max( uiPosQ - uiMaxScorePosQ, uiPosR - uiMaxScorePosR );
            int iZDropCurr = iMaxScore - iCurrScore - uiDiff * iExtend;
            if( iZDropCurr >= (int)uiZDropInversion )
                fDo( uiPosQ, uiPosR );
        } // for
    } // method

    inline std::shared_ptr<Alignment> tryInversionExtension( nucSeqIndex uiQFrom, std::shared_ptr<NucSeq> pQuery,
                                                             std::shared_ptr<NucSeq> pRef,
                                                             AlignedMemoryManager& rMemoryManager )
    {
        Wrapper_ksw_extz_t ez;
        kswcpp_dispatch( pQuery->length( ) - uiQFrom, pQuery->pxSequenceRef + uiQFrom, pRef->length( ),
                         pRef->pxSequenceRef, xKswParameters, uiBandwidth, (int)uiZDrop, KSW_EZ_EXTZ_ONLY, ez.ez,
                         rMemoryManager );
        auto pRet = std::make_shared<Alignment>( );

        nucSeqIndex qPos = uiQFrom;
        nucSeqIndex rPos = 0;
        for( int i = 0; i < ez.ez->n_cigar; ++i )
        {
            uint32_t uiSymbol = ez.ez->cigar[ i ] & 0xf;
            uint32_t uiAmount = ez.ez->cigar[ i ] >> 4;
            switch( uiSymbol )
            {
                case 0:
                    for( uint32_t uiPos = 0; uiPos < uiAmount; uiPos++ )
                    {
                        if( ( *pQuery )[ uiPos + qPos ] == ( *pRef )[ uiPos + rPos ] )
                            pRet->append( MatchType::match );
                        else
                            pRet->append( MatchType::missmatch );
                    } // for
                    qPos += uiAmount;
                    rPos += uiAmount;
                    break;
                case 1:
                    pRet->append( MatchType::insertion, uiAmount );
                    qPos += uiAmount;
                    break;
                case 2:
                    pRet->append( MatchType::deletion, uiAmount );
                    rPos += uiAmount;
                    break;
                default:
                    std::cerr << "obtained wierd symbol from ksw: " << uiSymbol << std::endl;
                    assert( false );
                    break;
            } // switch
        } // for

        return pRet;
    } // method

    // overload
    virtual std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>>
    execute( std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> pAlignments, std::shared_ptr<NucSeq> pQuery,
             std::shared_ptr<Pack> pRefPack )
    {
        AlignedMemoryManager xMemoryManager;
        // holds the additionally computed alignments
        ContainerVector<std::shared_ptr<Alignment>> vExtra;
        auto pRevQuery = std::make_shared<NucSeq>( );
        pRevQuery->vAppend( pQuery->pxSequenceRef, pQuery->length( ) );

        for( std::shared_ptr<Alignment> pAlignment : *pAlignments )
            forAllDropPos( //
                [&]( nucSeqIndex uiPosQ, nucSeqIndex uiPosR ) //
                {
                    bool bRevStrand = pRefPack->bPositionIsOnReversStrand( pAlignment->uiBeginOnRef );
                    auto uiRPos = pRefPack->uiPositionToReverseStrand( uiPosR );
                    auto pInvAlignment = tryInversionExtension( pQuery->length( ) - uiPosQ,
                                                                bRevStrand ? pQuery : pRevQuery,
                                                                pRefPack->vExtract( uiRPos, uiRPos + uiPosQ ),
                                                                xMemoryManager );
                    pInvAlignment->uiBeginOnQuery += uiPosQ;
                    pInvAlignment->uiEndOnQuery += uiPosQ;
                    pInvAlignment->uiBeginOnRef += uiRPos;
                    pInvAlignment->uiEndOnRef += uiRPos;
                    pInvAlignment->bSupplementary = true;
                    pInvAlignment->xStats = pAlignment->xStats;
                    pInvAlignment->fMappingQuality = 0;
                    if( bDisableHeuristics || pInvAlignment->score( ) > iHarmScoreMin * iMatch )
                        vExtra.push_back( pInvAlignment );
                }, // lambda
                pAlignment );
        for( std::shared_ptr<Alignment> pAlignment : vExtra )
            pAlignments->push_back( pAlignment );
        return pAlignments;
    } // function

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief Exposes the SmallInversions Module to python.
 * @ingroup export
 */
void exportSmallInversions( py::module& rxPyModuleId );
#endif
