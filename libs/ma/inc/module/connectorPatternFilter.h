/**
 * @file connectorPatternFilter.h
 */
#pragma once

#include "module/needlemanWunsch.h"
#include "module/sweepSvJumps.h"

namespace libMA
{
/**
 * @brief filters SV calls out that occur due to ambiguities on the reference
 * @details
 * Uses DP to check if the reference before and after the breakpoints matches better to itself than the reads to it.
 * If so it discards the SV.
 */
class ConnectorPatternFilter
    : public Module<CompleteBipartiteSubgraphClusterVector, false, CompleteBipartiteSubgraphClusterVector, Pack>
{
  public:
    nucSeqIndex uiMaxExtensionSize = 100;
    const KswCppParam<5> xKswParameters;
    const size_t uiZDrop;
    SQLQuery<DBCon, NucSeqSql> xGetRead;

    ConnectorPatternFilter( const ParameterSetManager& rParameters, std::shared_ptr<SV_Schema<DBCon>> pDBEnv )
        : uiZDrop( rParameters.getSelected( )->xZDrop->get( ) ),
          xGetRead( pDBEnv->pDatabase, "SELECT sequence FROM read_table WHERE id = ? " )
    {} // constructor
    std::shared_ptr<CompleteBipartiteSubgraphClusterVector>
    execute( std::shared_ptr<CompleteBipartiteSubgraphClusterVector> pCalls, std::shared_ptr<Pack> pRef )
    {
        AlignedMemoryManager xMemoryManager;

        auto pRet = std::make_shared<CompleteBipartiteSubgraphClusterVector>( );

        for( auto pCall : pCalls->vContent )
        {
            int64_t iReferenceScore = 0;
            int64_t iReadScore = 0;

            // compute the reference score before the jump
            uint64_t iBegin =
                pCall->xXAxis.start( ) > uiMaxExtensionSize ? pCall->xXAxis.start( ) - uiMaxExtensionSize : 0;
            uint64_t iSize = pCall->xXAxis.start( ) > uiMaxExtensionSize ? uiMaxExtensionSize : pCall->xXAxis.start( );
            if( pRef->bridgingSubsection( iBegin, iSize ) )
                pRef->unBridgeSubsection( iBegin, iSize );
            auto pNucSeqLeft = pRef->vExtract( iBegin, iBegin + iSize );
            pNucSeqLeft->vReverse( );

            iBegin = pCall->xXAxis.end();
            iSize = uiMaxExtensionSize;
            if( pRef->bridgingSubsection( iBegin, iSize ) )
                pRef->unBridgeSubsection( iBegin, iSize );
            auto pNucSeqRight = pRef->vExtract( iBegin, iBegin + iSize );

            iBegin = pCall->xYAxis.start( ) > uiMaxExtensionSize ? pCall->xYAxis.start( ) - uiMaxExtensionSize : 0;
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
                kswcpp_dispatch( (int)pNucSeqLeft->length( ), pNucSeqLeft->pxSequenceRef, (int)pNucSeqDown->length( ),
                                 pNucSeqDown->pxSequenceRef, xKswParameters, 100, (int)uiZDrop, KSW_EZ_EXTZ_ONLY, ez.ez,
                                 xMemoryManager );
                if( ez.ez->score > 0 )
                    iReferenceScore += ez.ez->score;
            } // scope for ez
            {
                Wrapper_ksw_extz_t ez;
                kswcpp_dispatch( (int)pNucSeqRight->length( ), pNucSeqRight->pxSequenceRef, (int)pNucSeqUp->length( ),
                                 pNucSeqUp->pxSequenceRef, xKswParameters, 100, (int)uiZDrop, KSW_EZ_EXTZ_ONLY, ez.ez,
                                 xMemoryManager );
                if( ez.ez->score > 0 )
                    iReferenceScore += ez.ez->score;
            } // scope for ez


            // compute the read scores
            for( auto pJump : pCall->vSupportingJumps )
            {
                // this db read puts quite the strain on the system...
                auto pRead = xGetRead.scalar( pJump->iReadId ).pNucSeq;
                assert( pRead != nullptr );
                pRead->iId = pJump->iReadId;

                int64_t iBegin = pJump->uiQueryFrom > uiMaxExtensionSize ? pJump->uiQueryFrom - uiMaxExtensionSize : 0;
                int64_t iSize = pJump->uiQueryFrom > uiMaxExtensionSize ? uiMaxExtensionSize : pJump->uiQueryFrom;
                auto pNucSeqLeft = std::make_shared<NucSeq>( pRead->fromTo( iBegin, iBegin + iSize ) );
                pNucSeqLeft->vReverse( );

                iBegin = pJump->uiQueryTo;
                iSize = pJump->uiQueryTo + uiMaxExtensionSize < pRead->length( ) ? uiMaxExtensionSize
                                                                                 : pRead->length( ) - pJump->uiQueryTo;
                auto pNucSeqRight = std::make_shared<NucSeq>( pRead->fromTo( iBegin, iBegin + iSize ) );

                {
                    Wrapper_ksw_extz_t ez;
                    kswcpp_dispatch( (int)pNucSeqLeft->length( ), pNucSeqLeft->pxSequenceRef, (int)pNucSeqUp->length( ),
                                     pNucSeqUp->pxSequenceRef, xKswParameters, 100, (int)uiZDrop, KSW_EZ_EXTZ_ONLY,
                                     ez.ez, xMemoryManager );
                    if( ez.ez->score > 0 )
                        iReadScore += ez.ez->score;
                } // scope for ez
                {
                    Wrapper_ksw_extz_t ez;
                    kswcpp_dispatch( (int)pNucSeqRight->length( ), pNucSeqRight->pxSequenceRef,
                                     (int)pNucSeqDown->length( ), pNucSeqDown->pxSequenceRef, xKswParameters, 100,
                                     (int)uiZDrop, KSW_EZ_EXTZ_ONLY, ez.ez, xMemoryManager );
                    if( ez.ez->score > 0 )
                        iReadScore += ez.ez->score;
                } // scope for ez

                // @todo optimization: this loop does not need to be executed all the way...
            } // for
            iReadScore /= pCall->vSupportingJumps.size( );

            // if the score for the read is better than the one purely on the reference, keep the call
            if( iReadScore > iReferenceScore )
                pRet->vContent.push_back( pCall );
        } // for

        return pRet;
    } // method
}; // class


} // namespace libMA

#ifdef WITH_PYTHON
void exportConnectorPatternFilter( py::module& rxPyModuleId );
#endif
