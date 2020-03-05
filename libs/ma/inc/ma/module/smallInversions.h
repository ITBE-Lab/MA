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
class SmallInversions : public libMS::Module<libMS::ContainerVector<std::shared_ptr<Alignment>>, false,
                                             libMS::ContainerVector<std::shared_ptr<Alignment>>, NucSeq, Pack>
{

    const KswCppParam<5> xKswParameters;
    const nucSeqIndex uiMaxGapArea;
    const size_t uiZDrop;
    const size_t uiZDropInversion;
    const size_t uiBandwidth;
    const bool bDisableHeuristics;
    const int iHarmScoreMin;
    const double dHarmScoreMinRel;

  public:
    bool bLocal = false;

    SmallInversions( const ParameterSetManager& rParameters )
        : xKswParameters( pGlobalParams->iMatch->get( ),
                          pGlobalParams->iMissMatch->get( ),
                          pGlobalParams->iGap->get( ),
                          pGlobalParams->iExtend->get( ),
                          pGlobalParams->iGap2->get( ),
                          pGlobalParams->iExtend2->get( ) ),
          uiMaxGapArea( rParameters.getSelected( )->xMaxGapArea->get( ) ),
          uiZDrop( rParameters.getSelected( )->xZDrop->get( ) ),
          uiZDropInversion( rParameters.getSelected( )->xZDropInversion->get( ) ),
          uiBandwidth( rParameters.getSelected( )->xBandwidthDPExtension->get( ) ),
          bDisableHeuristics( rParameters.getSelected( )->xDisableHeuristics->get( ) ),
          iHarmScoreMin( rParameters.getSelected( )->xHarmScoreMin->get( ) ),
          dHarmScoreMinRel( rParameters.getSelected( )->xHarmScoreMinRel->get( ) ){}; // default constructor

    inline void forAllDropPos( std::function<void( nucSeqIndex, nucSeqIndex, nucSeqIndex, nucSeqIndex )> fDo,
                               std::shared_ptr<Alignment>
                                   pAlignment )
    {
        nucSeqIndex uiMaxScorePosQ = pAlignment->uiBeginOnQuery;
        nucSeqIndex uiPosQ = uiMaxScorePosQ;
        nucSeqIndex uiStartQ = uiMaxScorePosQ;
        nucSeqIndex uiMaxScorePosR = pAlignment->uiBeginOnRef;
        nucSeqIndex uiPosR = uiMaxScorePosR;
        nucSeqIndex uiStartR = uiMaxScorePosR;
        int iMaxScore = std::numeric_limits<int>::min( );
        int iCurrScore = 0;
        int iMaxDrop = 0;
        for( std::pair<MatchType, nucSeqIndex> section : pAlignment->data )
        {
            switch( section.first )
            {
                case MatchType::seed:
                    if( iMaxDrop >= (int)uiZDropInversion )
                    {
                        // std::cout << "INVERSION: " << uiStartQ << " - " << uiPosQ << std::endl;
                        fDo( uiStartQ, uiStartR, uiPosQ, uiPosR );
                    } // if
                    uiStartQ = section.second + uiPosQ;
                    uiStartR = section.second + uiPosR;
                    iMaxDrop = 0;
                    iCurrScore = 0;
                    iMaxScore = std::numeric_limits<int>::min( );
                case MatchType::match:
                    iCurrScore += pGlobalParams->iMatch->get( ) * (int)section.second;
                    uiPosQ += section.second;
                    uiPosR += section.second;
                    break;
                case MatchType::missmatch:
                    iCurrScore -= pGlobalParams->iMissMatch->get( ) * (int)section.second;
                    uiPosQ += section.second;
                    uiPosR += section.second;
                    break;
                case MatchType::insertion:
                    iCurrScore -= pGlobalParams->iGap->get( ) + pGlobalParams->iExtend->get( ) * (int)section.second;
                    uiPosQ += section.second;
                    break;
                case MatchType::deletion:
                    iCurrScore -= pGlobalParams->iGap->get( ) + pGlobalParams->iExtend->get( ) * (int)section.second;
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
            else
            {
                int uiDiff = (int)std::max( uiPosQ - uiMaxScorePosQ, uiPosR - uiMaxScorePosR );
                int iZDropCurr = iMaxScore - iCurrScore - uiDiff * pGlobalParams->iExtend->get( );
                iMaxDrop = std::max( iMaxDrop, iZDropCurr );
            } // else
        } // for
    } // method

    inline std::shared_ptr<Alignment> tryInversionExtension( nucSeqIndex uiQFrom,
                                                             nucSeqIndex uiQTo,
                                                             std::shared_ptr<NucSeq>
                                                                 pQuery,
                                                             std::shared_ptr<NucSeq>
                                                                 pRef,
                                                             AlignedMemoryManager& rMemoryManager )
    {
        Wrapper_ksw_extz_t ez;
        kswcpp_dispatch( (int)uiQTo - (int)uiQFrom, pQuery->pGetSequenceRef( ) + uiQFrom, (int)pRef->length( ),
                         pRef->pGetSequenceRef( ), xKswParameters, (int)uiBandwidth, (int)uiZDrop, 0, ez.ez,
                         rMemoryManager );
        auto pRet = std::make_shared<Alignment>( );

        // std::cout << pQuery->length( ) - uiQFrom << " | " << pRef->length( ) << " | " << uiBandwidth << std::endl;

        // std::cout << "Query: " << uiQFrom << " " << pQuery->toString() << std::endl;
        // std::cout << "Ref: " << pRef->toString() << std::endl;
        // std::cout << "len: " << ez.ez->n_cigar << std::endl;
        // std::cout << "score: " << ez.ez->score << std::endl;
        // std::cout << "zdropped: " << ez.ez->zdropped << std::endl;

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
    virtual std::shared_ptr<libMS::ContainerVector<std::shared_ptr<Alignment>>>
    execute( std::shared_ptr<libMS::ContainerVector<std::shared_ptr<Alignment>>> pAlignments,
             std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRefPack )
    {
        AlignedMemoryManager xMemoryManager;
        // the container vector of alignments that is returned.
        auto pRet = std::make_shared<libMS::ContainerVector<std::shared_ptr<Alignment>>>( );
        for( std::shared_ptr<Alignment> pAlignment : *pAlignments )
        {
            pRet->push_back( pAlignment );
            forAllDropPos( //
                [&]( nucSeqIndex uiStartQ, nucSeqIndex uiStartR, nucSeqIndex uiEndQ, nucSeqIndex uiEndR ) {
                    auto uiStartRRevStr = pRefPack->uiPositionToReverseStrand( uiEndR );
                    auto uiEndRRevStr = pRefPack->uiPositionToReverseStrand( uiStartR );
                    auto pRef = pRefPack->vExtract( uiStartRRevStr, uiEndRRevStr );
                    auto pInvAlignment = tryInversionExtension( uiStartQ, uiEndQ, pQuery, pRef, xMemoryManager );
                    if( bDisableHeuristics || pInvAlignment->score( ) > iHarmScoreMin * pGlobalParams->iMatch->get( ) )
                    {
                        pInvAlignment->uiBeginOnQuery += uiStartQ;
                        pInvAlignment->uiEndOnQuery += uiStartQ;
                        pInvAlignment->uiBeginOnRef += uiStartRRevStr;
                        pInvAlignment->uiEndOnRef += uiStartRRevStr;
                        pInvAlignment->bSupplementary = true;
                        pInvAlignment->xStats = pAlignment->xStats;
                        pInvAlignment->fMappingQuality = 0;
                        pRet->push_back( pInvAlignment );
                    } // if
                }, // lambda
                pAlignment );
        } // for
        return pRet;
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
