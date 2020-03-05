/**
 * @file statisticSequenceAnalysis.h
 * @brief Implements statistical ways to analyze repetitiveness in sequences.
 * @author Markus Schmidt
 */
#include "geom.h"
#include "container/nucSeq.h"

#pragma once

namespace libMA
{

/**
 * @brief Determine the appropriate k-mers size for a "rectangle"
 * @details The formula used over here is:
 * 1 - t <= (1 - 1/4^k)^( (w-k+1)*(h-k+1) )
 * where (1 - 1/4^k) is the probability that two k-sized nucleotide sequences do not match.
 *	     (w-k+1)*(h-k+1) ) is the number of possible K-mer combinations within the rectangle.
 */
nucSeqIndex getKMerSizeForRectangle( geom::Rectangle<nucSeqIndex>& rRect, double t );

/**
 * @brief returns the size at which all k-mer's on the reference interval of xArea are unique.
 * @details
 * Currently implemented inefficiently.
 */
nucSeqIndex sampleKMerSize( NucSeq& rSequenceA, NucSeq& rSequenceB, double t );
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
nucSeqIndex sampleSequenceAmbiguity( NucSeq& rSequenceA, NucSeq& rSequenceB, double t );
inline nucSeqIndex sampleSequenceAmbiguity( NucSeq& rSequence, double t )
{
    NucSeq xEmpty;
    return sampleSequenceAmbiguity( rSequence, xEmpty, t );
} // function

} // namespace libMA