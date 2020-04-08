#include "msv/util/statisticSequenceAnalysis.h"
#include "ma/module/hashMapSeeding.h"


nucSeqIndex libMSV::getKMerSizeForRectangle( geom::Rectangle<nucSeqIndex>& rRect, double t )
{
    auto w = rRect.xXAxis.size( ); // width of rectangle
    auto h = rRect.xYAxis.size( ); // height of rectangle

    int64_t iDenominator = 4 * 4;
    for( nucSeqIndex uiK = 2; uiK <= std::min( w, h ); uiK++ )
    {
        // p is the probability of NOT having random match in the rectange (wh, w).
        // So, if p is small, then the probability that we have a random match is high.
        // If p is large, the probabilty of having a random match is low.
        double p = std::pow( 1.0 - ( 1.0 / ( (double)iDenominator ) ), ( w - uiK + 1 ) * ( h - uiK + 1 ) );
        if( p >= 1 - t )
            // The probabilty for a match of length uiK is now below 1 - t.
            return uiK;
        iDenominator *= 4;
    } // for

    // Give up (we cannot reach the required probabilty) and return an impossible k-mer size with respect to matching.
    return std::min( w, h ) + 1;
} // function

nucSeqIndex libMSV::sampleKMerSize( NucSeq& rSequenceA, NucSeq& rSequenceB, double t )
{
    geom::Rectangle<nucSeqIndex> xRect( 0, 0, rSequenceA.length( ) + rSequenceB.length( ),
                                        rSequenceA.length( ) + rSequenceB.length( ) );
    nucSeqIndex uiStaticsticalSize = getKMerSizeForRectangle( xRect, t );
    nucSeqIndex uiSeedSize = uiStaticsticalSize;
    for( ; uiSeedSize < rSequenceA.length( ) || uiSeedSize < rSequenceB.length( ); uiSeedSize++ )
    {
        std::set<std::string> xKMerSet;
        bool bAllUnique = true;
        for( auto& rSequence : {rSequenceA, rSequenceB} )
            for( nucSeqIndex uiPos = 0; uiPos + uiSeedSize <= rSequence.length( ); uiPos++ )
            {
                std::string sKMer = rSequence.fromTo( uiPos, uiPos + uiSeedSize );
                if( xKMerSet.find( sKMer ) != xKMerSet.end( ) )
                {
                    bAllUnique = false;
                    break;
                } // if
                xKMerSet.insert( sKMer );
            } // for
        if( bAllUnique )
            break;
    } // while

    return uiSeedSize - uiStaticsticalSize;
} // function

nucSeqIndex libMSV::sampleSequenceAmbiguity( NucSeq& rSequenceA, NucSeq& rSequenceB, double t )
{
    HashMapSeeding xSeeder;

    geom::Rectangle<nucSeqIndex> xRect( 0, 0, rSequenceA.length( ) + rSequenceB.length( ),
                                        rSequenceA.length( ) + rSequenceB.length( ) );
    xSeeder.uiSeedSize = getKMerSizeForRectangle( xRect, t );

    auto pSeeds = xSeeder.execute( rSequenceA, rSequenceB );
    pSeeds->append( xSeeder.execute( rSequenceA, rSequenceA ) );
    pSeeds->append( xSeeder.execute( rSequenceB, rSequenceB ) );

    auto pLumped = SeedLumping( ).execute( pSeeds, rSequenceA, rSequenceB );

    nucSeqIndex uiSum = 0;
    for( auto xSeed : *pLumped )
        uiSum += xSeed.size( );

    return uiSum;
} // function