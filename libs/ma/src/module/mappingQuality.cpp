/**
 * @file mappingQuality.cpp
 * @author Markus Schmidt
 */
#include "module/mappingQuality.h"

using namespace libMA;
using namespace libMS;


std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> MappingQuality::execute(
    std::shared_ptr<NucSeq> pQuery, std::shared_ptr<ContainerVector<std::shared_ptr<Alignment>>> pAlignments )
{
    std::sort(
        pAlignments->begin( ),
        pAlignments->end( ),
        []( std::shared_ptr<Alignment>& pA, std::shared_ptr<Alignment>& pB ) { return pA->score( ) > pB->score( ); } );

    auto pSupplementaries = std::make_shared<ContainerVector<std::shared_ptr<Alignment>>>( );

    // if no alignment was found we cannot set any quality...
    if( pAlignments->size( ) == 0 )
        return std::make_shared<ContainerVector<std::shared_ptr<Alignment>>>( );

    // compute the mapping quality for the best alignment
    std::shared_ptr<Alignment> pFirst = std::dynamic_pointer_cast<Alignment>( ( *pAlignments )[ 0 ] ); // dc
    pFirst->bSecondary = false;

    // set the mapping quality for all alignments to zero
    // mapping qual of best one will be overridden
    size_t uiSupplementaries = 0;
    bool bFirst = true;
    for( auto pAlign : *pAlignments )
    {
        if( bFirst )
        {
            bFirst = false;
            continue;
        } // if
        std::shared_ptr<Alignment> pCasted = std::dynamic_pointer_cast<Alignment>( pAlign ); // dc
        pCasted->fMappingQuality = 0.0;
        if( uiSupplementaries < uiMaxSupplementaryPerPrim && pCasted->overlap( *pFirst ) < dMaxOverlapSupplementary )
        {
            pCasted->bSupplementary = true;
            pCasted->bSecondary = false;
            uiSupplementaries++;
        } // if
        else
        {
            pCasted->bSupplementary = false;
            pCasted->bSecondary = true;
        } // else
    } // for


    // mapping quality based on scores
    if( pAlignments->size( ) - uiSupplementaries >= 2 )
    {
        size_t uiI = 1;
        std::shared_ptr<Alignment> pSecond;
        while( pSecond == nullptr || pSecond->bSupplementary )
        {
            pSecond = std::dynamic_pointer_cast<Alignment>( // dc
                ( *pAlignments )[ uiI ] );
            uiI++;
        } // while
        // this formula is given in the paper and is very similar to Heng li's approach in BWA-SW
        assert( pFirst->score( ) >= pSecond->score( ) );
        if( pFirst->score( ) == 0 )
            pFirst->fMappingQuality = 0;
        else
        {
            pFirst->fMappingQuality =
                static_cast<double>( pFirst->score( ) - pSecond->score( ) ) / static_cast<double>( pFirst->score( ) );
        } // else
    } // if
    else
        // the score of the second best alignment is 0 if we do not even find one...
        pFirst->fMappingQuality = pFirst->score( ) / (double)( pGlobalParams->iMatch->get() * pQuery->length( ) );

    if(pFirst->getNumSeeds() <= 1)
        pFirst->fMappingQuality /= 2;

    if(pFirst->score() >= pGlobalParams->iMatch->get() * pQuery->length( ) * 0.8 && pAlignments->size() >= 3)
        pFirst->fMappingQuality *= 2;
    if(pFirst->fMappingQuality > 1)
        pFirst->fMappingQuality = 1;

    assert( !pFirst->xStats.bSetMappingQualityToZero );
    if( pFirst->xStats.bSetMappingQualityToZero )
        pFirst->fMappingQuality = 0;

    if( uiSupplementaries > 0 )
    {
        // set supp mapping quality
        for( size_t uiI = 1; uiI < pAlignments->size( ); uiI++ )
        {
            std::shared_ptr<Alignment> pCasted = std::dynamic_pointer_cast<Alignment>( // dc
                ( *pAlignments )[ uiI ] );
            if( pCasted->bSupplementary )
                pCasted->fMappingQuality = pFirst->fMappingQuality;
        } // for
        // move supplementary alignments forward
        std::sort( pAlignments->begin( ), pAlignments->end( ),
                   []( std::shared_ptr<Alignment> a, std::shared_ptr<Alignment> b ) { return a->larger( b ); } // lambda
        ); // sort function call
    } // if

    // factors
    // penalty for too little seeds
    // (this improves the mapping quality estimation quite significantly)
    // double dA = std::max(std::min(5* pFirst->numBySeeds() / (double)pQuery->length(), 1.0),
    // 0.01); pFirst->fMappingQuality *= dA;

    /// maybe this should be moved into it's own module but whatever...
    auto pRet = std::make_shared<ContainerVector<std::shared_ptr<Alignment>>>( pAlignments );
    if( uiReportNBest != 0 && pRet->size( ) > uiReportNBest + uiSupplementaries )
    {
        // remove the smallest elements
        pRet->erase( pRet->begin( ) + uiReportNBest + uiSupplementaries, pRet->end( ) );
        assert( pRet->size( ) == uiReportNBest + uiSupplementaries );
    } // if

    // remove alignments below the minimal score
    pRet->erase( std::remove_if(
                     pRet->begin( ), pRet->end( ),
                     [&]( const std::shared_ptr<Alignment>& pA ) { return pA->score( ) < (long)uiMinAlignmentScore; } ),
                 pRet->end( ) );

    return pRet;
} // function

#ifdef WITH_PYTHON

void exportMappingQuality( py::module& rxPyModuleId )
{
    // export the MappingQuality class
    exportModule<MappingQuality>( rxPyModuleId, "MappingQuality" );
} // function
#endif
