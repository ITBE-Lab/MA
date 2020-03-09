/**
 * @file connectorPatternFilter.h
 */
#pragma once

#include "ms/container/sv_db/pool_container.h"
#include "ma/module/needlemanWunsch.h"
#include "msv/module/sweepSvJumps.h"

namespace libMSV
{

/**
 * @brief filters SV calls out that occur due to ambiguities on the reference
 * @details
 * Uses DP to check if the reference before and after the breakpoints matches better to itself than the reads to it.
 * If so it discards the SV.
 */
template <typename DBCon>
class ConnectorPatternFilter : public Module<CompleteBipartiteSubgraphClusterVector, false,
                                             CompleteBipartiteSubgraphClusterVector, Pack, libMS::PoolContainer<DBCon>>
{
  public:
    nucSeqIndex uiMaxExtensionSize = 100;
    const KswCppParam<5> xKswParameters;
    const size_t uiZDrop;

    ConnectorPatternFilter( const ParameterSetManager& rParameters )
        : uiZDrop( rParameters.getSelected( )->xZDrop->get( ) )
    {} // constructor

    std::shared_ptr<CompleteBipartiteSubgraphClusterVector>
    execute( std::shared_ptr<CompleteBipartiteSubgraphClusterVector> pCalls, std::shared_ptr<Pack> pRef,
             std::shared_ptr<libMS::PoolContainer<DBCon>> pPool )
    {
        return pPool->xPool.run(
            [this]( auto pConnection, std::shared_ptr<CompleteBipartiteSubgraphClusterVector> pCalls,
                    std::shared_ptr<Pack> pRef ) {
                SQLQuery<DBCon, std::shared_ptr<CompressedNucSeq>> xGetRead(
                    pConnection, "SELECT sequence FROM read_table WHERE id = ? " );
                AlignedMemoryManager xMemoryManager;

                auto pRet = std::make_shared<CompleteBipartiteSubgraphClusterVector>( );

                for( auto pCall : pCalls->vContent )
                {
                    int64_t iReferenceScore = 0;
                    int64_t iReadScore = 0;

                    // compute the reference score before the jump
                    uint64_t iBegin =
                        pCall->xXAxis.start( ) > uiMaxExtensionSize ? pCall->xXAxis.start( ) - uiMaxExtensionSize : 0;
                    uint64_t iSize =
                        pCall->xXAxis.start( ) > uiMaxExtensionSize ? uiMaxExtensionSize : pCall->xXAxis.start( );
                    if( pRef->bridgingSubsection( iBegin, iSize ) )
                        pRef->unBridgeSubsection( iBegin, iSize );
                    auto pNucSeqLeft = pRef->vExtract( iBegin, iBegin + iSize );
                    pNucSeqLeft->vReverse( );

                    iBegin = pCall->xXAxis.end( );
                    iSize = uiMaxExtensionSize;
                    if( pRef->bridgingSubsection( iBegin, iSize ) )
                        pRef->unBridgeSubsection( iBegin, iSize );
                    auto pNucSeqRight = pRef->vExtract( iBegin, iBegin + iSize );

                    iBegin =
                        pCall->xYAxis.start( ) > uiMaxExtensionSize ? pCall->xYAxis.start( ) - uiMaxExtensionSize : 0;
                    iSize = pCall->xYAxis.start( ) > uiMaxExtensionSize ? uiMaxExtensionSize : pCall->xYAxis.start( );
                    if( pRef->bridgingSubsection( iBegin, iSize ) )
                        pRef->unBridgeSubsection( iBegin, iSize );
                    auto pNucSeqDown = pRef->vExtract( iBegin, iBegin + iSize );
                    pNucSeqDown->vReverse( );

                    iBegin = pCall->xYAxis.end( );
                    iSize = uiMaxExtensionSize;
                    if( pRef->bridgingSubsection( iBegin, iSize ) )
                        pRef->unBridgeSubsection( iBegin, iSize );
                    auto pNucSeqUp = pRef->vExtract( iBegin, iBegin + iSize );

                    if( pCall->bSwitchStrand )
                    {
                        pNucSeqUp.swap( pNucSeqDown );
                        pNucSeqUp->vSwitchAllBasePairsToComplement( );
                        pNucSeqDown->vSwitchAllBasePairsToComplement( );
                    } // if

                    {
                        Wrapper_ksw_extz_t ez;
                        kswcpp_dispatch( (int)pNucSeqLeft->length( ), pNucSeqLeft->pxSequenceRef,
                                         (int)pNucSeqDown->length( ), pNucSeqDown->pxSequenceRef, xKswParameters, 100,
                                         (int)uiZDrop, KSW_EZ_EXTZ_ONLY, ez.ez, xMemoryManager );
                        if( ez.ez->score > 0 )
                            iReferenceScore += ez.ez->score;
                    } // scope for ez
                    {
                        Wrapper_ksw_extz_t ez;
                        kswcpp_dispatch( (int)pNucSeqRight->length( ), pNucSeqRight->pxSequenceRef,
                                         (int)pNucSeqUp->length( ), pNucSeqUp->pxSequenceRef, xKswParameters, 100,
                                         (int)uiZDrop, KSW_EZ_EXTZ_ONLY, ez.ez, xMemoryManager );
                        if( ez.ez->score > 0 )
                            iReferenceScore += ez.ez->score;
                    } // scope for ez


                    // compute the read scores
                    for( auto pJump : pCall->vSupportingJumps )
                    {
                        // this db select puts quite the strain on the system...
                        auto pRead = xGetRead.scalar( pJump->iReadId )->pUncomNucSeq;
                        assert( pRead != nullptr );
                        pRead->iId = pJump->iReadId;

                        int64_t iBegin =
                            pJump->uiQueryFrom > uiMaxExtensionSize ? pJump->uiQueryFrom - uiMaxExtensionSize : 0;
                        int64_t iSize =
                            pJump->uiQueryFrom > uiMaxExtensionSize ? uiMaxExtensionSize : pJump->uiQueryFrom;
                        auto pNucSeqLeft = std::make_shared<NucSeq>( pRead->fromTo( iBegin, iBegin + iSize ) );
                        pNucSeqLeft->vReverse( );

                        iBegin = pJump->uiQueryTo;
                        iSize = pJump->uiQueryTo + uiMaxExtensionSize < pRead->length( )
                                    ? uiMaxExtensionSize
                                    : pRead->length( ) - pJump->uiQueryTo;
                        auto pNucSeqRight = std::make_shared<NucSeq>( pRead->fromTo( iBegin, iBegin + iSize ) );

                        {
                            Wrapper_ksw_extz_t ez;
                            kswcpp_dispatch( (int)pNucSeqLeft->length( ), pNucSeqLeft->pxSequenceRef,
                                             (int)pNucSeqUp->length( ), pNucSeqUp->pxSequenceRef, xKswParameters, 100,
                                             (int)uiZDrop, KSW_EZ_EXTZ_ONLY, ez.ez, xMemoryManager );
                            if( ez.ez->score > 0 )
                                iReadScore += ez.ez->score;
                        } // scope for ez
                        {
                            Wrapper_ksw_extz_t ez;
                            kswcpp_dispatch( (int)pNucSeqRight->length( ), pNucSeqRight->pxSequenceRef,
                                             (int)pNucSeqDown->length( ), pNucSeqDown->pxSequenceRef, xKswParameters,
                                             100, (int)uiZDrop, KSW_EZ_EXTZ_ONLY, ez.ez, xMemoryManager );
                            if( ez.ez->score > 0 )
                                iReadScore += ez.ez->score;
                        } // scope for ez

                        if( iReadScore / (int64_t)pCall->vSupportingJumps.size( ) > iReferenceScore )
                            break;
                    } // for

                    // if the score for the read is better than the one purely on the reference, keep the call
                    if( iReadScore / (int64_t)pCall->vSupportingJumps.size( ) > iReferenceScore )
                        pRet->vContent.push_back( pCall );
                } // for

                return pRet;
            },
            pCalls, pRef );
    } // method
}; // class


} // namespace libMSV

#ifdef WITH_PYTHON
void exportConnectorPatternFilter( libMS::SubmoduleOrganizer& xOrganizer );
#endif
