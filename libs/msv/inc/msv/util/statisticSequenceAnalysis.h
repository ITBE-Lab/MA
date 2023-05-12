/**
 * @file statisticSequenceAnalysis.h
 * @brief Implements statistical ways to analyze repetitiveness in sequences.
 * @author Markus Schmidt
 */
#include "ma/container/nucSeq.h"
#include "ma/container/pack.h"
#include "util/geom.h"

#pragma once

using namespace libMA;
namespace libMSV
{

/**
 * @brief Determine the appropriate k-mers size for a "rectangle"
 * @details The formula used over here is:
 * 1 - t <= (1 - 1/4^k)^( (w-k+1)*(h-k+1) )
 * where (1 - 1/4^k) is the probability that two k-sized nucleotide sequences do not match.
 *	     (w-k+1)*(h-k+1) ) is the number of possible K-mer combinations within the rectangle.
 */
nucSeqIndex DLL_PORT( MSV ) getKMerSizeForRectangle( geom::Rectangle<nucSeqIndex>& rRect, double t );

/**
 * @brief returns the size at which all k-mer's on the reference interval of xArea are unique.
 * @details
 * Currently implemented inefficiently.
 */
nucSeqIndex DLL_PORT( MSV ) sampleKMerSize( NucSeq& rSequenceA, NucSeq& rSequenceB, double t );
inline nucSeqIndex sampleKMerSize( NucSeq& rSequence, double t )
{
    NucSeq xEmpty;
    return sampleKMerSize( rSequence, xEmpty, t );
} // function

/**
 * @brief returns the number of non random nucleotide chains matching among the sequences
 * @details
 * Computes and lumps k-mers, so that all chains have lengths that have a likelyhood <= t to occur randomly.
 * Compares both sequences with themselves and each other, so the minimal returned value should be
 * the length of both sequences added together.
 */
nucSeqIndex DLL_PORT( MSV ) sampleSequenceAmbiguity( NucSeq& rSequenceA, NucSeq& rSequenceB, double t );
inline nucSeqIndex sampleSequenceAmbiguity( NucSeq& rSequence, double t )
{
    NucSeq xEmpty;
    return sampleSequenceAmbiguity( rSequence, xEmpty, t );
} // function


inline nucSeqIndex DLL_PORT( MSV ) sampleAmbiguity( std::shared_ptr<NucSeq> pSeqA, std::shared_ptr<NucSeq> pSeqB )
{
    return std::max( 1, (int)sampleSequenceAmbiguity( *pSeqA, *pSeqB, 0.001 ) - (int)pSeqA->length( ) -
                            (int)pSeqB->length( ) );
} // method


inline std::shared_ptr<NucSeq> DLL_PORT( MSV )
    getRegion( nucSeqIndex uiPos, bool bLeftDirection, std::shared_ptr<Pack> pPack, nucSeqIndex uiDistance )
{
    // due to their fuzziness calls can reach past the end of the genome
    if( uiPos >= pPack->uiUnpackedSizeForwardStrand )
        uiPos = pPack->uiUnpackedSizeForwardStrand - 1;

    auto uiSeqId = pPack->uiSequenceIdForPosition( uiPos );
    if( bLeftDirection )
    {
        nucSeqIndex iStartOfContig = pPack->startOfSequenceWithId( uiSeqId );
        nucSeqIndex uiStart = uiPos > iStartOfContig + uiDistance ? uiPos - uiDistance : iStartOfContig;
        nucSeqIndex uiSize = uiPos - uiStart;
        // return empty sequence for size = 0 cause pack throws exception otherwise
        if( uiSize == 0 )
            return std::make_shared<NucSeq>( );
        if( pPack->bridgingSubsection( uiStart, uiSize ) )
            pPack->unBridgeSubsection( uiStart, uiSize );
        return pPack->vExtract( uiStart, uiStart + uiSize );
    } // if
    else
    {
        nucSeqIndex iEndOfContig = pPack->endOfSequenceWithId( uiSeqId );
        nucSeqIndex uiEnd = uiPos + uiDistance < iEndOfContig ? uiPos + uiDistance : iEndOfContig;
        nucSeqIndex uiSize = uiEnd - uiPos;
        // return empty sequence for size = 0 cause pack throws exception otherwise
        if( uiSize == 0 )
            return std::make_shared<NucSeq>( );
        if( pPack->bridgingSubsection( uiPos, uiSize ) )
            pPack->unBridgeSubsection( uiPos, uiSize );
        return pPack->vExtract( uiPos, uiPos + uiSize );
    } // else
} // method

inline nucSeqIndex DLL_PORT( MSV )
    sampleSequenceAmbiguity( nucSeqIndex uiFrom, nucSeqIndex uiTo, bool bFromForward, bool bToForward,
                             std::shared_ptr<Pack> pPack, nucSeqIndex uiDistanceMax, nucSeqIndex uiDistanceMin )
{
    auto uiCallSize = uiFrom >= uiTo ? uiFrom - uiTo : uiTo - uiFrom;

    if( uiCallSize > uiDistanceMin || bFromForward != bToForward )
    {
        auto uiDistance = std::min( uiCallSize, uiDistanceMax );
        auto pLeftFrom = getRegion( uiFrom, true, pPack, uiDistance );
        auto pRightFrom = getRegion( uiFrom, false, pPack, uiDistance );
        auto pLeftTo = getRegion( uiTo, true, pPack, uiDistance );
        auto pRightTo = getRegion( uiTo, false, pPack, uiDistance );

        // if we switch strand we have to compare forward and reverse strands
        if( bFromForward != bToForward )
        {
            pLeftTo->vReverseAll( );
            pLeftTo->vSwitchAllBasePairsToComplement( );
            pRightTo->vReverseAll( );
            pRightTo->vSwitchAllBasePairsToComplement( );
        } // if

        auto a = sampleAmbiguity( pLeftFrom, ( bFromForward != bToForward ) ? pRightTo : pLeftTo );
        auto b = sampleAmbiguity( pRightFrom, ( bFromForward != bToForward ) ? pLeftTo : pRightTo );

        return std::max( a, b );
    } // if
    else
        // @todo how to evaluate such calls?
        return 1;
}

} // namespace libMSV