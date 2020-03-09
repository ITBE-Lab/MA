/**
 * @file otherSeeding.cpp
 * @brief implements two futher seeding techniques for texting purpose.
 * @author Markus Schmidt
 */
#include "ma/module/otherSeeding.h"

using namespace libMA;
using namespace libMS;

#define INCLUDE_SELF_TESTS ( 1 )

#include <atomic>
#include <chrono>
#include <memory>
#include <vector>

#define complement( x ) ( uint8_t ) NucSeq::nucleotideComplement( x )


void OtherSeeding::bowtieExtension( std::shared_ptr<FMIndex> pFM_index,
                                    std::shared_ptr<NucSeq>
                                        pQuerySeq,
                                    std::shared_ptr<SegmentVector>
                                        pSegmentVector )
{
    unsigned int iSize = 16;
    unsigned int iStep = 1;

    const uint8_t* q = pQuerySeq->pGetSequenceRef( );
    for( nucSeqIndex i = 0; i < pQuerySeq->length( ) - iSize; i += iStep )
    {
        SAInterval ik( pFM_index->L2[ complement( q[ i ] ) ] + 1, pFM_index->L2[ (int)q[ i ] ] + 1,
                       pFM_index->L2[ (int)q[ i ] + 1 ] - pFM_index->L2[ (int)q[ i ] ] );
        for( nucSeqIndex i2 = 1; i2 <= iSize; i2++ )
        {
            // this is the extension
            // actually forward since we extend backwards on the reverse complement...
            ik = pFM_index->extend_backward( ik, complement( q[ i2 + i ] ) );
            if( ik.size( ) == 0 )
                break;
        } // for
        if( ik.size( ) == 0 )
            continue;
        // found a seed
        pSegmentVector->emplace_back( i, iSize, ik.revComp( ) );
    } // for
} // function

/**
 * This is what blasr does:
 *
 * for each position in the query:
 *  extend maximally backwards
 *  save a seed that is one shorter than the maximal extension
 *  if the seed is longer than K = 12 adn maxambiguity = ?
 */
void OtherSeeding::doBlasrExtension( std::shared_ptr<FMIndex> pFM_index,
                                     std::shared_ptr<NucSeq>
                                         pQuerySeq,
                                     std::shared_ptr<SegmentVector>
                                         pSegmentVector )
{
    const uint8_t* q = pQuerySeq->pGetSequenceRef( );
    for( nucSeqIndex i = 0; i < pQuerySeq->length( ); i += 1 )
    {
        SAInterval ik( pFM_index->L2[ (int)q[ i ] ] + 1, pFM_index->L2[ complement( q[ i ] ) ] + 1,
                       pFM_index->L2[ (int)q[ i ] + 1 ] - pFM_index->L2[ (int)q[ i ] ] );
        SAInterval lk( 0, 0, 0 ); // will hold the maximal extension
        SAInterval llk( 0, 0, 0 ); // will hold the interval one shorter than the maximal extension
        nucSeqIndex i2 = 0;
        while( i2 <= i )
        {
            llk = lk;
            lk = ik;
            ik = pFM_index->extend_backward( ik, q[ i - i2 ] );
            if( ik.size( ) == 0 )
                break;
            i2++;
        } // for
        if( i2 <= 12 )
            continue;
        // found a seed
        pSegmentVector->emplace_back( i - i2 + 1, i2 - 1, llk );
    } // for
} // function

std::shared_ptr<SegmentVector> OtherSeeding::execute( std::shared_ptr<FMIndex> pFM_index,
                                                      std::shared_ptr<NucSeq> pQuerySeq )
{

    std::shared_ptr<SegmentVector> pSegmentVector( new SegmentVector( ) );
    if( pQuerySeq == nullptr )
        return pSegmentVector;

    DEBUG_2( std::cout << pQuerySeq->fastaq( ) << std::endl; )

    if( bBowtie )
        bowtieExtension( pFM_index, pQuerySeq, pSegmentVector );
    else
        doBlasrExtension( pFM_index, pQuerySeq, pSegmentVector );

    return pSegmentVector;
} // function


#ifdef WITH_PYTHON

void exportOtherSeeding( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the OtherSeeding class
    exportModule<OtherSeeding, bool>( xOrganizer, "OtherSeeding" );
} // function
#endif
