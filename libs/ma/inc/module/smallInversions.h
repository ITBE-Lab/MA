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

    inline void forAllDropPos( std::function<std::pair<nucSeqIndex, nucSeqIndex>( nucSeqIndex, nucSeqIndex )> fDo,
                               std::shared_ptr<Alignment>
                                   pAlignment )
    {
        nucSeqIndex uiMaxScorePosQ = pAlignment->uiBeginOnQuery;
        nucSeqIndex uiPosQ = uiMaxScorePosQ;
        nucSeqIndex uiMaxScorePosR = pAlignment->uiBeginOnRef;
        nucSeqIndex uiPosR = uiMaxScorePosR;
        nucSeqIndex uiInvEndQ = 0;
        nucSeqIndex uiInvEndR = 0;
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
                    iCurrScore -= iMissMatch * section.second;
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
            else if( uiPosQ >= uiInvEndQ && uiPosR >= uiInvEndR )
            {
                size_t uiDiff = std::max( uiPosQ - uiMaxScorePosQ, uiPosR - uiMaxScorePosR );
                int iZDropCurr = iMaxScore - iCurrScore - uiDiff * iExtend;
                if( iZDropCurr >= (int)uiZDropInversion )
                {
                    std::cout << "INVERSION: " << uiPosQ << ", " << uiPosR << std::endl;
                    auto xInvEnd = fDo( uiPosQ, uiPosR );
                    uiInvEndQ = xInvEnd.first;
                    uiInvEndR = xInvEnd.second;
                } // if
            } // else
        } // for
    } // method

    inline std::shared_ptr<Alignment> tryInversionExtension( nucSeqIndex uiQFrom, std::shared_ptr<NucSeq> pQuery,
                                                             std::shared_ptr<NucSeq> pRef,
                                                             AlignedMemoryManager& rMemoryManager )
    {
        Wrapper_ksw_extz_t ez;
        kswcpp_dispatch( pQuery->length( ) - uiQFrom, pQuery->pGetSequenceRef( ) + uiQFrom, pRef->length( ),
                         pRef->pGetSequenceRef( ), xKswParameters, -1, -1, KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR, ez.ez,
                         rMemoryManager );
        auto pRet = std::make_shared<Alignment>( );

        std::cout << pQuery->length( ) - uiQFrom << " | " << pRef->length( ) << " | " << uiBandwidth << std::endl;

        // std::cout << "Query: " << uiQFrom << " " << pQuery->toString() << std::endl;
        // std::cout << "Ref: " << pRef->toString() << std::endl;
        std::cout << "len: " << ez.ez->n_cigar << std::endl;
        std::cout << "score: " << ez.ez->score << std::endl;
        std::cout << "zdropped: " << ez.ez->zdropped << std::endl;

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
        // the container vector of alignments that is returned.
        auto pRet = std::make_shared<ContainerVector<std::shared_ptr<Alignment>>>( );
        // auto pRevQuery = std::make_shared<NucSeq>( );
        // pRevQuery->vAppend( pQuery->pxSequenceRef, pQuery->length( ) );
        // pRevQuery->vReverse();

        for( std::shared_ptr<Alignment> pAlignment : *pAlignments )
        {
            pRet->push_back( pAlignment );
            forAllDropPos( //
                [&]( nucSeqIndex uiPosQ, nucSeqIndex uiPosR ) -> std::pair<nucSeqIndex, nucSeqIndex> {
                    auto uiREnd = pRefPack->uiPositionToReverseStrand( uiPosR );
                    auto pRef = pRefPack->vExtract( uiREnd - ( pQuery->length( ) - uiPosQ ) , uiREnd );
                    pRef->vReverse();
                    auto pInvAlignment =
                        tryInversionExtension( uiPosQ,
                                               pQuery,
                                               pRef,
                                               xMemoryManager );
                    pInvAlignment->uiBeginOnQuery += uiPosQ;
                    pInvAlignment->uiEndOnQuery += uiPosQ;
                    auto uiRStart = pRefPack->uiPositionToReverseStrand( uiPosR + pInvAlignment->uiEndOnRef );
                    pInvAlignment->uiBeginOnRef = uiRStart;
                    pInvAlignment->uiEndOnRef = uiREnd;
                    pInvAlignment->bSupplementary = true;
                    pInvAlignment->xStats = pAlignment->xStats;
                    pInvAlignment->fMappingQuality = 0;
                    if( bDisableHeuristics || pInvAlignment->score( ) > iHarmScoreMin * iMatch )
                        pRet->push_back( pInvAlignment );
                    return std::make_pair( pInvAlignment->uiEndOnQuery, pInvAlignment->uiBeginOnRef );
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
