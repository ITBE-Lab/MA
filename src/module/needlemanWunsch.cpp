/**
 * @file needlemanWunsch.cpp
 * @author Markus Schmidt
 */

#define OLD_KSW ( 0 )

#include "module/needlemanWunsch.h"

#if OLD_KSW == 1
#include "ksw/ksw2.h"
#endif

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

using namespace libMA;


inline void ksw_ext( int qlen,
                     const uint8_t* query,
                     int tlen,
                     const uint8_t* target,
                     const KswCppParam<5>& rKswParameters,
                     int w,
                     size_t uiZDrop,
#if OLD_KSW == 1
                     ksw_extz_t* ez,
#else
                     kswcpp_extz_t* ez,
#endif
                     AlignedMemoryManager& rMemoryManager,
                     bool bRef )
{
#if OLD_KSW == 1
    if( bRef )
        ksw_extd2_sse_( nullptr, qlen, query, tlen, target, 5, &rKswParameters.mat[ 0 ], rKswParameters.q,
                       rKswParameters.e, rKswParameters.q2, rKswParameters.e2, w, uiZDrop, -1,
                       KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR, ez );
    else
        ksw_extd2_sse_( nullptr, qlen, query, tlen, target, 5, &rKswParameters.mat[ 0 ], rKswParameters.q,
                       rKswParameters.e, rKswParameters.q2, rKswParameters.e2, w, uiZDrop, -1, KSW_EZ_EXTZ_ONLY, ez );
#else
    if( bRef )
        kswcpp_dispatch( qlen, query, tlen, target, rKswParameters, w, (int)uiZDrop,
                         KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR, ez, rMemoryManager );
    else
        kswcpp_dispatch( qlen, query, tlen, target, rKswParameters, w, (int)uiZDrop, KSW_EZ_EXTZ_ONLY, ez, rMemoryManager );
#endif
} // function

inline void ksw_simplified( int qlen, const uint8_t* query, int tlen, const uint8_t* target,
                            const KswCppParam<5>& rKswParameters, int w,
#if OLD_KSW == 1
                            ksw_extz_t* ez,
#else
                            kswcpp_extz_t* ez,
#endif
                            AlignedMemoryManager& rMemoryManager )
{
    int minAddBandwidth = 10; // must be >= 0 otherwise ksw will not align till the end
    /*
     * Adjust the bandwith according to the delta distance of the seeds creating this gap
     * the add minAddBandwidth so that the alignment can go a little further out.
     */
    if( std::abs( tlen - qlen ) + minAddBandwidth > w )
        w = std::abs( tlen - qlen ) + minAddBandwidth;

#if OLD_KSW == 1
    ksw_extd2_sse_( nullptr, qlen, query, tlen, target, 5, &rKswParameters.mat[ 0 ], rKswParameters.q, rKswParameters.e,
                   rKswParameters.q2, rKswParameters.e2, w, -1, -1, 0, ez );
#else
    kswcpp_dispatch( qlen, query, tlen, target, rKswParameters, w, -1, 0, ez, rMemoryManager );
#endif
} // function


// small wrapper that takes care of deallocation
class Wrapper_ksw_extz_t
{
  public:
#if OLD_KSW == 1
    ksw_extz_t* ez;
#else
    kswcpp_extz_t* ez;
#endif

    Wrapper_ksw_extz_t( )
    {
#if OLD_KSW == 1
        ez = new ksw_extz_t{}; // {} forces zero initialization
#else
        ez = new kswcpp_extz_t{}; // {} forces zero initialization
#endif

    } // default constructor

    ~Wrapper_ksw_extz_t( )
    {
        free( ez->cigar ); // malloced in c code
        delete ez; // allocated by new in cpp code
    } // default constructor
}; // class

// banded global NW
void NeedlemanWunsch::ksw( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<NucSeq> pRef, nucSeqIndex fromQuery,
                           nucSeqIndex toQuery, nucSeqIndex fromRef, nucSeqIndex toRef,
                           std::shared_ptr<Alignment> pAlignment, AlignedMemoryManager& rMemoryManager )
{
    // sanity checks
    if( toRef <= fromRef )
        if( toQuery <= fromQuery )
            return;
    if( toQuery <= fromQuery )
    {
        pAlignment->append( MatchType::deletion, toRef - fromRef );
        return;
    } // if
    if( toRef <= fromRef )
    {
        pAlignment->append( MatchType::insertion, toQuery - fromQuery );
        return;
    } // if

    Wrapper_ksw_extz_t ez;

    assert( toQuery < pQuery->length( ) );
    assert( toRef < pRef->length( ) );
    ksw_simplified( (int)(toQuery - fromQuery), pQuery->pGetSequenceRef( ) + fromQuery, (int) (toRef - fromRef),
                    pRef->pGetSequenceRef( ) + fromRef, xKswParameters, iMinBandwidthGapFilling,
                    ez.ez, // return value
                    rMemoryManager );

    nucSeqIndex qPos = fromQuery;
    nucSeqIndex rPos = fromRef;
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
                        pAlignment->append( MatchType::match );
                    else
                        pAlignment->append( MatchType::missmatch );
                } // for
                qPos += uiAmount;
                rPos += uiAmount;
                break;
            case 1:
                pAlignment->append( MatchType::insertion, uiAmount );
                qPos += uiAmount;
                break;
            case 2:
                pAlignment->append( MatchType::deletion, uiAmount );
                rPos += uiAmount;
                break;
            default:
                std::cerr << "obtained wierd symbol from ksw: " << uiSymbol << std::endl;
                assert( false );
                break;
        } // switch
    } // for
#if DEBUG_LEVEL >= 1
    const char vMIDN[] = {'M', 'I', 'D', 'N'};
    if( qPos != (uint32_t)toQuery && rPos != (uint32_t)toRef )
    {
        std::cerr << "ksw did neither extend till end of query nor ref " << (int)toQuery - qPos << "; "
                  << (int)toRef - rPos << " bandwidth: " << iMinBandwidthGapFilling << std::endl;
        std::cout << pQuery->fromTo( fromQuery, toQuery ) << std::endl;
        std::cout << pRef->fromTo( fromRef, toRef ) << std::endl;
        std::cout << "CIGAR:";
        for( int i = 0; i < ez.ez->n_cigar; ++i )
        {
            uint32_t uiSymb = ez.ez->cigar[ i ] & 0xf;
            uint32_t uiAmount = ez.ez->cigar[ i ] >> 4;
            std::cout << uiAmount << vMIDN[ uiSymb ] << " ";
        } // for
        std::cout << std::endl;
        assert( false );
    } // if
#endif
    // ensure we do not append something negative.. (underflow)
    assert( toQuery >= qPos );
    assert( toRef >= rPos );
    // Append the remaining insertion or deletion.
    // Remaining lengths of zero are taken care of in append
    pAlignment->append( MatchType::deletion, toQuery - qPos );
    pAlignment->append( MatchType::insertion, toRef - rPos );
} // function

class MyPrinterMemory
{
  public:
    std::shared_ptr<NucSeq> pQuery;
    std::shared_ptr<NucSeq> pRef;
    std::shared_ptr<Alignment> pAlignment;
    nucSeqIndex uiQPos, uiRPos;
    DEBUG( nucSeqIndex uiQueryExtensionSize = 0; nucSeqIndex uiRefExtensionSize = 0; )

    MyPrinterMemory( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<NucSeq> pRef,
                     std::shared_ptr<Alignment> pAlignment, nucSeqIndex uiQueryStart, nucSeqIndex uiRefStart )
        : pQuery( pQuery ), pRef( pRef ), pAlignment( pAlignment ), uiQPos( uiQueryStart ), uiRPos( uiRefStart )
    {} // constructor
}; // class

////// DEPRECATED
//// char vTrans[ 9 ] = {'?', 'A', 'C', '?', 'G', '?', '?', '?', 'T'};
//// int printer( void* pVoid, uint64_t uiLen, char c )
//// {
////     MyPrinterMemory* pM = ( MyPrinterMemory* )pVoid;
////     assert( pM != nullptr );
////     assert( pM->pAlignment != nullptr );
////     // std::cout << pM->pAlignment->length() << std::endl;
////     // std::cout << pM->pAlignment->data.size() << std::endl;
////     DEBUG_3( std::cout << c << uiLen << " "; )
////     switch ( c )
////     {
////         case 'M':
////             for ( size_t i = 0; i < uiLen; i++ )
////             {
////                 if ( ( *pM->pQuery )[ pM->uiQPos + i ] == ( *pM->pRef )[ pM->uiRPos + i ] )
////                     pM->pAlignment->append( MatchType::match );
////                 else
////                     pM->pAlignment->append( MatchType::missmatch );
////             } // for
////             pM->uiQPos += uiLen;
////             pM->uiRPos += uiLen;
////             DEBUG( pM->uiQueryExtensionSize += uiLen; pM->uiRefExtensionSize += uiLen; )
////             break;
////         case 'I':
////             pM->pAlignment->append( MatchType::insertion, uiLen );
////             pM->uiQPos += uiLen;
////             DEBUG( pM->uiQueryExtensionSize += uiLen; )
////             break;
////         case 'D':
////             pM->pAlignment->append( MatchType::deletion, uiLen );
////             pM->uiRPos += uiLen;
////             DEBUG( pM->uiRefExtensionSize += uiLen; )
////             break;
////         default:
////             std::cout << "GABA cigar contains unknown symbol: " << c << std::endl;
////             assert( false );
////     } // switch
////     return 0;
//// } // function


// banded global NW
void NeedlemanWunsch::ksw_dual_ext( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<NucSeq> pRef, nucSeqIndex fromQuery,
                                    nucSeqIndex toQuery, nucSeqIndex fromRef, nucSeqIndex toRef,
                                    std::shared_ptr<Alignment> pAlignment, AlignedMemoryManager& rMemoryManager )
{
    Wrapper_ksw_extz_t ez_left;
    Wrapper_ksw_extz_t ez_right;

    ksw_ext( (int)(toQuery - fromQuery), pQuery->pGetSequenceRef( ) + fromQuery, (int)(toRef - fromRef),
             pRef->pGetSequenceRef( ) + fromRef, xKswParameters, iBandwidthDPExtension, uiZDrop,
             ez_left.ez, // return value
             rMemoryManager, false );

    pQuery->vReverse( fromQuery, toQuery );
    pRef->vReverse( fromRef, toRef );
    ksw_ext( (int)(toQuery - fromQuery), pQuery->pGetSequenceRef( ) + fromQuery, (int)(toRef - fromRef),
             pRef->pGetSequenceRef( ) + fromRef, xKswParameters, iBandwidthDPExtension, uiZDrop,
             ez_right.ez, // return value
             rMemoryManager, true );
    pQuery->vReverse( fromQuery, toQuery );
    pRef->vReverse( fromRef, toRef );

    // center between both extensions
    nucSeqIndex qCenter = ( fromQuery + ez_left.ez->max_q + ( toQuery - ez_right.ez->max_q - 1 ) ) / 2;
    // in case one of the extensions reaches to the other end but the other drops immediately
    qCenter = std::max( fromQuery, std::min( toQuery, qCenter ) );
    nucSeqIndex rCenter = ( fromRef + ez_left.ez->max_t + ( toRef - ez_right.ez->max_t - 1 ) ) / 2;
    // in case one of the extensions reaches to the other end but the other drops immediately
    rCenter = std::max( fromRef, std::min( toRef, rCenter ) );

    nucSeqIndex qPos = fromQuery;
    nucSeqIndex rPos = fromRef;

    if( rPos != rCenter && qPos != qCenter )
        for( int i = 0; i < ez_left.ez->n_cigar; ++i )
        {
            assert( qPos < qCenter );
            assert( rPos < rCenter );
            uint32_t uiSymbol = ez_left.ez->cigar[ i ] & 0xf;
            uint32_t uiAmount = ez_left.ez->cigar[ i ] >> 4;
            switch( uiSymbol )
            {
                case 0:
                    // dont go over the center
                    if( qPos + uiAmount > qCenter )
                    {
                        assert( qCenter >= qPos );
                        uiAmount = (uint32_t)(qCenter - qPos);
                    } // if
                    // dont go over the center
                    if( rPos + uiAmount > rCenter )
                    {
                        assert( rCenter >= rPos );
                        uiAmount = (uint32_t)(rCenter - rPos);
                    } // if
                    for( uint32_t uiPos = 0; uiPos < uiAmount; uiPos++ )
                    {
                        if( ( *pQuery )[ uiPos + qPos ] == ( *pRef )[ uiPos + rPos ] )
                            pAlignment->append( MatchType::match );
                        else
                            pAlignment->append( MatchType::missmatch );
                    } // for
                    qPos += uiAmount;
                    rPos += uiAmount;
                    break;
                case 1:
                    // dont go over the center
                    if( qPos + uiAmount > qCenter )
						uiAmount = (uint32_t)(qCenter - qPos);
                    pAlignment->append( MatchType::insertion, uiAmount );
                    qPos += uiAmount;
                    break;
                case 2:
                    // dont go over the center
                    if( rPos + uiAmount > rCenter )
						uiAmount = (uint32_t)(rCenter - rPos);
                    pAlignment->append( MatchType::deletion, uiAmount );
                    rPos += uiAmount;
                    break;
                default:
                    std::cerr << "obtained wierd symbol from ksw: " << uiSymbol << std::endl;
                    assert( false );
                    break;
            } // switch
            assert( rPos <= rCenter );
            assert( qPos <= qCenter );
            if( rPos == rCenter )
                break;
            if( qPos == qCenter )
                break;
        } // for
    assert( rPos <= rCenter );
    assert( qPos <= qCenter );
    assert( qPos == qCenter || rPos == rCenter || ez_left.ez->zdropped == 1 );
    nucSeqIndex rPosRight = toRef - ez_right.ez->max_t - 1;
    nucSeqIndex qPosRight = toQuery - ez_right.ez->max_q - 1;
    uint32_t uiAmountNotUnrolled = 0;
    // should be overwritten or unused anyways...
    MatchType xTypeLastUnrolledCigar = MatchType::seed;
    int i = 0;
    // unroll the cigar until both query and refernce positions are past the center
    // it might be necessary to unroll one cigar operation partially...
    for( ; i < ez_right.ez->n_cigar; ++i )
    {
        // if we are past both centerlines stop unrolling
        if( rPosRight >= rCenter && qPosRight >= qCenter )
            break;
        // uiAmountNotUnrolled should only be set in the very last iteration of this loop
        assert( uiAmountNotUnrolled == 0 );
        uint32_t uiSymbol = ez_right.ez->cigar[ i ] & 0xf;
        uint32_t uiAmount = ez_right.ez->cigar[ i ] >> 4;
        switch( uiSymbol )
        {
            case 0:
                if( rPosRight + uiAmount >= rCenter && qPosRight + uiAmount >= qCenter )
                {
                    if( rPosRight < rCenter && ( qPosRight >= qCenter || rCenter - rPosRight > qCenter - qPosRight ) )
                    {
                        assert( rCenter > rPosRight );
                        assert( uiAmount >= ( rCenter - rPosRight ) );
                        uiAmountNotUnrolled = uiAmount - (uint32_t)( rCenter - rPosRight );
                        uiAmount = (uint32_t)(rCenter - rPosRight);
                    } // if
                    else
                    {
                        assert( qCenter > qPosRight );
                        assert( uiAmount >= ( qCenter - qPosRight ) );
                        uiAmountNotUnrolled = uiAmount - (uint32_t)( qCenter - qPosRight );
                        uiAmount = (uint32_t)(qCenter - qPosRight);
                    } // else
                } // if
                qPosRight += uiAmount;
                rPosRight += uiAmount;
                xTypeLastUnrolledCigar = MatchType::match;
                break;
            case 1:
                if( qPosRight + uiAmount > qCenter && rPosRight >= rCenter )
                {
                    assert( qCenter > qPosRight );
                    assert( uiAmount >= ( qCenter - qPosRight ) );
                    uiAmountNotUnrolled = uiAmount - (uint32_t)( qCenter - qPosRight );
					uiAmount = (uint32_t)(qCenter - qPosRight);
                } // if
                qPosRight += uiAmount;
                xTypeLastUnrolledCigar = MatchType::insertion;
                break;
            case 2:
                if( rPosRight + uiAmount > rCenter && qPosRight >= qCenter )
                {
                    assert( rCenter > rPosRight );
                    assert( uiAmount >= ( rCenter - rPosRight ) );
                    uiAmountNotUnrolled = uiAmount - (uint32_t)( rCenter - rPosRight );
					uiAmount = (uint32_t)(rCenter - rPosRight);
                } // if
                rPosRight += uiAmount;
                xTypeLastUnrolledCigar = MatchType::deletion;
                break;
            default:
                std::cerr << "obtained wierd symbol from ksw: " << uiSymbol << std::endl;
                assert( false );
                break;
        } // switch
    } // for

    // fill in the gap between the left and right extension
    assert( rPosRight >= rPos );
    pAlignment->append( MatchType::deletion, rPosRight - rPos );
    assert( qPosRight >= qPos );
    pAlignment->append( MatchType::insertion, qPosRight - qPos );

    // add the last cigar operation (it might have been partially unrolled...)
    if( xTypeLastUnrolledCigar == MatchType::match )
        for( uint32_t uiPos = 0; uiPos < uiAmountNotUnrolled; uiPos++ )
        {
            if( ( *pQuery )[ uiPos + qPosRight ] == ( *pRef )[ uiPos + rPosRight ] )
                pAlignment->append( MatchType::match );
            else
                pAlignment->append( MatchType::missmatch );
        } // for
    else
        pAlignment->append( xTypeLastUnrolledCigar, uiAmountNotUnrolled );

    // move q & r pos forward
    switch( xTypeLastUnrolledCigar )
    {
        case MatchType::match:
            qPosRight += uiAmountNotUnrolled;
            rPosRight += uiAmountNotUnrolled;
            break;
        case MatchType::insertion:
            qPosRight += uiAmountNotUnrolled;
            break;
        case MatchType::deletion:
            rPosRight += uiAmountNotUnrolled;
            break;
        default: // if xTypeLastUnrolledCigar is unset: do nothing
            break;
    } // switch

    // add all remaining cigar operations
    for( ; i < ez_right.ez->n_cigar; ++i )
    {
        uint32_t uiSymbol = ez_right.ez->cigar[ i ] & 0xf;
        uint32_t uiAmount = ez_right.ez->cigar[ i ] >> 4;
        switch( uiSymbol )
        {
            case 0:
                for( uint32_t uiPos = 0; uiPos < uiAmount; uiPos++ )
                {
                    if( ( *pQuery )[ uiPos + qPosRight ] == ( *pRef )[ uiPos + rPosRight ] )
                        pAlignment->append( MatchType::match );
                    else
                        pAlignment->append( MatchType::missmatch );
                } // for
                qPosRight += uiAmount;
                rPosRight += uiAmount;
                break;
            case 1:
                pAlignment->append( MatchType::insertion, uiAmount );
                qPosRight += uiAmount;
                break;
            case 2:
                pAlignment->append( MatchType::deletion, uiAmount );
                rPosRight += uiAmount;
                break;
            default:
                std::cerr << "obtained wierd symbol from ksw: " << uiSymbol << std::endl;
                assert( false );
                break;
        } // switch
    } // for
    assert( qPosRight == toQuery );
    assert( rPosRight == toRef );
} // function

void NeedlemanWunsch::dynPrg( const std::shared_ptr<NucSeq> pQuery, const std::shared_ptr<NucSeq> pRef,
                              const nucSeqIndex fromQuery, const nucSeqIndex toQuery, const nucSeqIndex fromRef,
                              const nucSeqIndex toRef,
                              std::shared_ptr<Alignment> pAlignment, // in & output
                              AlignedMemoryManager& rMemoryManager, const bool bLocalBeginning, const bool bLocalEnd )
{
    // do some checking for empty sequences
    if( toRef <= fromRef )
        if( toQuery <= fromQuery )
        {
            DEBUG_3( std::cout << "dynProg end" << std::endl; )
            return;
        } // if
    if( toQuery <= fromQuery )
    {
        pAlignment->append( MatchType::deletion, toRef - fromRef );
        DEBUG_3( std::cout << "dynProg end" << std::endl; )
        return;
    } // if
    if( toRef <= fromRef )
    {
        pAlignment->append( MatchType::insertion, toQuery - fromQuery );
        DEBUG_3( std::cout << "dynProg end" << std::endl; )
        return;
    } // if

    // do not actually compute through gaps that are larger than a set maximum
#if 1
    if( toQuery - fromQuery > uiMaxGapArea || toRef - fromRef > uiMaxGapArea )
    {
        ksw_dual_ext( pQuery, pRef, fromQuery, toQuery, fromRef, toRef, pAlignment, rMemoryManager );
        return;
    } // if
#endif

    DEBUG_3( std::cout << "dynProg begin" << std::endl; )
    if( !bLocalBeginning && !bLocalEnd )
    {
        ksw( pQuery, pRef, fromQuery, toQuery, fromRef, toRef, pAlignment, rMemoryManager );
        DEBUG_3( std::cout << "dynProg end" << std::endl; )
        return;
    } // if

    DEBUG_3( std::cout << "sw1" << std::endl; )

    assert( !( bLocalBeginning && bLocalEnd ) );
    assert( bLocalBeginning || bLocalEnd );

    const bool bReverse = bLocalBeginning;

    // if we reached this point we actually have to align something
    DEBUG_2( std::cout << pQuery->toString( ) << std::endl; std::cout << pRef->toString( ) << std::endl; )

    DEBUG_3( std::cout << "sw2" << std::endl; )
    /*
     * do the NW alignment
     */

    Wrapper_ksw_extz_t ez;

    assert( toQuery < pQuery->length( ) );
    assert( toRef < pRef->length( ) );
    if( bReverse )
    {
        pQuery->vReverse( fromQuery, toQuery );
        pRef->vReverse( fromRef, toRef );
    } // if
    ksw_ext( (int)(toQuery - fromQuery), pQuery->pGetSequenceRef( ) + fromQuery, (int)(toRef - fromRef),
             pRef->pGetSequenceRef( ) + fromRef, xKswParameters, iBandwidthDPExtension, uiZDrop,
             ez.ez, // return value
             rMemoryManager, bReverse );
    if( bReverse )
    {
        pQuery->vReverse( fromQuery, toQuery );
        pRef->vReverse( fromRef, toRef );
    } // if

    nucSeqIndex qPos = fromQuery;
    nucSeqIndex rPos = fromRef;
    if( bReverse )
    {
        rPos += toRef - ez.ez->max_t - 1;
        qPos += toQuery - ez.ez->max_q - 1;
    } // if

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
                        pAlignment->append( MatchType::match );
                    else
                        pAlignment->append( MatchType::missmatch );
                } // for
                qPos += uiAmount;
                rPos += uiAmount;
                break;
            case 1:
                pAlignment->append( MatchType::insertion, uiAmount );
                qPos += uiAmount;
                break;
            case 2:
                pAlignment->append( MatchType::deletion, uiAmount );
                rPos += uiAmount;
                break;
            default:
                std::cerr << "obtained wierd symbol from ksw: " << uiSymbol << std::endl;
                assert( false );
                break;
        } // switch
    } // for
    /*
     * Warning:
     * Order is important: the shifting needs to be done after the cigar extraction
     */
    if( bReverse )
    {
        pAlignment->shiftOnRef( toRef - ez.ez->max_t - 1 );
        pAlignment->shiftOnQuery( toQuery - ez.ez->max_q - 1 );
    } // if

} // function


std::shared_ptr<Alignment> NeedlemanWunsch::execute_one( std::shared_ptr<Seeds> pSeeds, std::shared_ptr<NucSeq> pQuery,
                                                         std::shared_ptr<Pack> pRefPack,
                                                         AlignedMemoryManager& rMemoryManager )
{

    if( pSeeds == nullptr )
        return std::shared_ptr<Alignment>( new Alignment( ) );

    // no seeds => no spot found at all...
    if( pSeeds->empty( ) )
    {
        DEBUG( std::cerr << "WARNING: no seeds found for query: " + pQuery->sName << std::endl; )
        std::shared_ptr<Alignment> pRet( new Alignment( ) );
        pRet->xStats = pSeeds->xStats;
        pRet->xStats.sName = pQuery->sName;
        return pRet;
    } // if
#if DEBUG_LEVEL >= 2
    std::cout << "seedlist: (start_ref, end_ref; start_query, end_query)" << std::endl;
    for( Seed& rSeed : *pSeeds )
    {
        std::cout << rSeed.start_ref( ) << ", " << rSeed.end_ref( ) << "; " << rSeed.start( ) << ", " << rSeed.end( )
                  << std::endl;
    } // for
#endif

    // Determine the query and reverence coverage of the seeds

    nucSeqIndex beginRef = pSeeds->front( ).start_ref( );
    nucSeqIndex endRef = pSeeds->back( ).end_ref( );
    // seeds are sorted by ther startpos so we
    // actually need to check all seeds to get the proper end
    nucSeqIndex endQuery = pSeeds->back( ).end( );
    nucSeqIndex beginQuery = pSeeds->front( ).start( );
    for( auto xSeed : *pSeeds )
    {
        if( endRef < xSeed.end_ref( ) )
            endRef = xSeed.end_ref( );
        if( beginRef > xSeed.start_ref( ) )
            beginRef = xSeed.start_ref( );
        if( endQuery < xSeed.end( ) )
            endQuery = xSeed.end( );
        if( beginQuery > xSeed.end( ) )
            beginQuery = xSeed.start( );
        assert( xSeed.start( ) <= xSeed.end( ) );
    } // for
    DEBUG_2( std::cout << beginRef << ", " << endRef << "; " << beginQuery << ", " << endQuery << std::endl; ) // DEEBUG

    if( beginRef >= endRef || pRefPack->bridgingSubsection( beginRef, endRef - beginRef ) )
    {
#if 0
        // sometimes we can save the situation by making the last seed smaller...
        int64_t iContig = pRefPack->uiSequenceIdForPositionOrRev(pSeeds->back().start_ref());
        if( iContig == pRefPack->uiSequenceIdForPositionOrRev(beginRef) )
        {
            pSeeds->back().size( 
                pRefPack->endOfSequenceWithIdOrReverse(iContig) - pSeeds->back().start_ref() - 1
             );
            endRef = pSeeds->back().end_ref();
            for (auto xSeed : *pSeeds)
                if(xSeed.end_ref() > endRef)
                {
                    xSeed.size(endRef - xSeed.start_ref());
                    assert(endRef >= xSeed.end_ref());
                }// if
        }// if
        else
#endif
        {
#if CONTIG_ID_CACHE == ( 1 )
            DEBUG( std::cerr << "WARNING: computed bridging alignment:\n";
                   std::cerr << beginRef << " - " << endRef << std::endl;
                   std::cerr << pRefPack->nameOfSequenceForPosition( beginRef ) << " - "
                             << pRefPack->nameOfSequenceForPosition( endRef ) << std::endl;
                   std::cerr << pRefPack->iAbsolutePosition( beginRef ) << " - "
                             << pRefPack->iAbsolutePosition( endRef ) << std::endl;
                   auto names = pRefPack->contigNames( );
                   auto starts = pRefPack->contigStarts( );
                   auto lengths = pRefPack->contigLengths( );
                   for( size_t i = 0; i < names.size( ); i++ ) {
                       if( starts[ i ] + lengths[ i ] >= (nucSeqIndex)pRefPack->iAbsolutePosition( beginRef ) ||
                           starts[ i ] + lengths[ i ] >= (nucSeqIndex)pRefPack->iAbsolutePosition( endRef ) )
                           std::cerr << names[ i ] << ": [" << starts[ i ] << "-" << starts[ i ] + lengths[ i ]
                                     << "] revComp: [" << pRefPack->uiPositionToReverseStrand( starts[ i ] ) << "-"
                                     << pRefPack->uiPositionToReverseStrand( starts[ i ] + lengths[ i ] ) << "]"
                                     << std::endl;
                       if( starts[ i ] >= (nucSeqIndex)pRefPack->iAbsolutePosition( beginRef ) &&
                           starts[ i ] >= (nucSeqIndex)pRefPack->iAbsolutePosition( endRef ) )
                           break;
                   } // for
                   ) // DEBUG
#endif
            std::shared_ptr<Alignment> pRet( new Alignment( ) );
            pRet->xStats = pSeeds->xStats;
            pRet->xStats.sName = pQuery->sName;
            return pRet;
        } // else
    } // if

    // here we have enough query coverage to attemt to fill in the gaps merely
    DEBUG_2( std::cout << "filling in gaps" << std::endl; )
    assert( !pRefPack->bridgingSubsection( beginRef, endRef - beginRef ) );

    std::shared_ptr<Alignment> pRet;

#define PADDING_W_Q ( 0 )

    if( !bLocal )
    {
        int64_t iOldContig = pRefPack->uiSequenceIdForPositionOrRev( beginRef );
#if PADDING_W_Q == 1
        beginRef -= uiPadding + beginQuery;
#else
        beginRef -= uiPadding;
#endif
        if( beginRef > endRef ) // check for underflow
            beginRef = 0;
#if PADDING_W_Q == 1
        endRef += uiPadding + ( pQuery->length( ) - endQuery );
#else
        endRef += uiPadding;
#endif
        if( beginRef > endRef ) // check for overflow
            endRef = pRefPack->uiUnpackedSizeForwardPlusReverse( );
        endQuery = pQuery->length( );
        beginQuery = 0;
        if( pRefPack->uiSequenceIdForPositionOrRev( beginRef ) != iOldContig )
            beginRef = pRefPack->startOfSequenceWithIdOrReverse( iOldContig );
        if( pRefPack->uiSequenceIdForPositionOrRev( endRef ) != iOldContig )
            endRef = pRefPack->endOfSequenceWithIdOrReverse( iOldContig );

        DEBUG( if( beginRef >= endRef || pRefPack->bridgingSubsection( beginRef, endRef - beginRef ) ) {
            std::cerr << "ERROR: produced bridging alignment:\n";
            std::cerr << beginRef << " - " << endRef << std::endl;
            std::cerr << pRefPack->nameOfSequenceForPosition( beginRef ) << " - "
                      << pRefPack->nameOfSequenceForPosition( endRef ) << std::endl;
            std::cerr << pRefPack->iAbsolutePosition( beginRef ) << " - " << pRefPack->iAbsolutePosition( endRef )
                      << std::endl;
            auto names = pRefPack->contigNames( );
            auto starts = pRefPack->contigStarts( );
            auto lengths = pRefPack->contigLengths( );
            for( size_t i = 0; i < names.size( ); i++ )
            {
                if( starts[ i ] + lengths[ i ] >= (nucSeqIndex)pRefPack->iAbsolutePosition( beginRef ) ||
                    starts[ i ] + lengths[ i ] >= (nucSeqIndex)pRefPack->iAbsolutePosition( endRef ) )
                    std::cerr << names[ i ] << ": [" << starts[ i ] << "-" << starts[ i ] + lengths[ i ]
                              << "] revComp: [" << pRefPack->uiPositionToReverseStrand( starts[ i ] ) << "-"
                              << pRefPack->uiPositionToReverseStrand( starts[ i ] + lengths[ i ] ) << "]" << std::endl;
                if( starts[ i ] >= (nucSeqIndex)pRefPack->iAbsolutePosition( beginRef ) &&
                    starts[ i ] >= (nucSeqIndex)pRefPack->iAbsolutePosition( endRef ) )
                    break;
            } // for
        } // if
        )
        assert( !pRefPack->bridgingSubsection( beginRef, endRef - beginRef ) );
        assert( beginRef <= pSeeds->front( ).start_ref( ) );
    } // if
    assert( endQuery <= pQuery->length( ) );
    pRet = std::shared_ptr<Alignment>( new Alignment( beginRef, beginQuery ) );

    // save the strip of consideration stats in the alignment
    pRet->xStats = pSeeds->xStats;
    pRet->xStats.sName = pQuery->sName;

    DEBUG_2( std::cout << beginRef << " " << endRef << std::endl; )
    std::shared_ptr<NucSeq> pRef = pRefPack->vExtract( beginRef, endRef );

    // create the actual alignment

    if( !bLocal )
    {
        dynPrg( pQuery, pRef, 0, pSeeds->front( ).start( ), 0, pSeeds->front( ).start_ref( ) - beginRef, pRet,
                rMemoryManager, true, false );
    } // if

    nucSeqIndex endOfLastSeedQuery = pSeeds->front( ).end( );
    nucSeqIndex endOfLastSeedReference = pSeeds->front( ).end_ref( ) - beginRef;

    DEBUG( if( pRet->uiEndOnQuery != pSeeds->front( ).start( ) ) {
        std::cout << pRet->uiEndOnQuery << " ?= " << pSeeds->front( ).start( ) << std::endl;
        std::cout << pRet->uiEndOnRef << " ?= " << pSeeds->front( ).start_ref( ) << std::endl;
        assert( false );
    } // if
           if( pRet->uiEndOnRef != pSeeds->front( ).start_ref( ) ) {
               std::cout << pRet->uiEndOnQuery << " ?= " << pSeeds->front( ).start( ) << std::endl;
               std::cout << pRet->uiEndOnRef << " ?= " << pSeeds->front( ).start_ref( ) << std::endl;
               assert( false );
           } // if
           ) // DEBUG

    pRet->append( MatchType::seed, pSeeds->front( ).size( ) );
    bool bSkip = true;
    for( Seed& rSeed : *pSeeds )
    {
        // skip the first seed
        // we do this since the seed has already been appended before the loop
        // this makes the loop structure easier since this way
        // we can always first compute the NW and then append a seed
        if( bSkip )
        {
            bSkip = false;
            continue;
        } // if
        if( rSeed.size( ) == 0 )
            continue;
        nucSeqIndex ovQ = endOfLastSeedQuery - rSeed.start( );
        if( rSeed.start( ) > endOfLastSeedQuery )
            ovQ = 0;
        nucSeqIndex ovR = endOfLastSeedReference - ( rSeed.start_ref( ) - beginRef );
        if( rSeed.start_ref( ) > endOfLastSeedReference + beginRef )
            ovR = 0;
        nucSeqIndex len = rSeed.size( );
        nucSeqIndex overlap = std::max( ovQ, ovR );
        DEBUG_2( std::cout << "overlap: " << overlap << std::endl; ) // DEBUG
        if( len > overlap )
        {
            dynPrg( pQuery, pRef, endOfLastSeedQuery, rSeed.start( ), endOfLastSeedReference,
                    rSeed.start_ref( ) - beginRef, pRet, rMemoryManager, false, false );
            DEBUG(
                // std::cout << pRet->vGapsScatter.size() << std::endl;
                pRet->vGapsScatter.push_back(
                    std::make_pair( rSeed.start( ) - endOfLastSeedQuery,
                                    rSeed.start_ref( ) - beginRef - endOfLastSeedReference ) ); ) // DEBUG
            if( ovQ > ovR )
                pRet->append( MatchType::deletion, ovQ - ovR );
            DEBUG_2( for( nucSeqIndex i = ovR; i < ovQ; i++ ) std::cout << "d"; )
            if( ovR > ovQ )
                pRet->append( MatchType::insertion, ovR - ovQ );
            DEBUG_2( for( nucSeqIndex i = ovQ; i < ovR; i++ ) std::cout << "i"; std::cout << std::endl; ) // DEBUG
            pRet->append( MatchType::seed, len - overlap );
            DEBUG_2( std::cout << len - overlap << std::endl; ) // DEBUG_2
            DEBUG_2( for( nucSeqIndex i = overlap; i < len; i++ ) std::cout << pQuery->charAt( i + rSeed.start( ) );
                     std::cout << std::endl;
                     for( nucSeqIndex i = overlap; i < len; i++ ) std::cout
                     << pRef->charAt( i + rSeed.start_ref( ) - beginRef );
                     std::cout << std::endl; ) // DEBUG
            DEBUG_2( for( nucSeqIndex i = 0; i < len - overlap; i++ ) std::cout << "m"; ) // DEBUG_2
            if( rSeed.end( ) > endOfLastSeedQuery )
                endOfLastSeedQuery = rSeed.end( );
            if( rSeed.end_ref( ) > endOfLastSeedReference + beginRef )
                endOfLastSeedReference = rSeed.end_ref( ) - beginRef;
        } // if
    } // for

    if( bLocal )
        assert( std::get<0>( pRet->data.front( ) ) == MatchType::seed );
    assert( std::get<0>( pRet->data.back( ) ) == MatchType::seed );

    DEBUG_2( std::cout << std::endl; )
    if( bLocal )
        pRet->makeLocal( );
    else
    {
        dynPrg( pQuery, pRef, endOfLastSeedQuery, endQuery - 1, endOfLastSeedReference, endRef - beginRef - 1, pRet,
                rMemoryManager, false, true );
        // there should never be dangeling deletions with libGaba
        pRet->removeDangeling( );
    } // else
    return pRet;
} // function


/* Random nucleotide sequence of length uiLen, represented as codes.
 */
std::vector<char> randomNucSeq( const size_t uiLen )
{
    static const char nucleotides[] = {0, 1, 2, 3};

    std::vector<char> vNucSeq( uiLen );
    for( size_t i = 0; i < uiLen; ++i )
    {
        vNucSeq[ i ] = nucleotides[ rand( ) % ( sizeof( nucleotides ) - 1 ) ];
    } // for

    return vNucSeq;
} // function

#if 0
void testKsw( )
{
    /* Seed the random number generator */
    auto uiSeed = (unsigned)time( NULL );
    //// std::cout << "SEED:" << uiSeed << std::endl;
    srand( uiSeed );

    size_t uiRefSize = 1;
    size_t uiQuerySize = 10000;

    int8_t iMatch = 1;
    int8_t iMissMatch = 4;
    int8_t iGap = 10;
    int8_t iExtend = 1;
    int iBandwidth = 100;

#if 0 // adaptive bandwidth
    int minAddBandwidth = 10; // must be >= 0 otherwise ksw will not align till the end
    if( std::abs(uiRefSize - uiQuerySize) > iBandwidth + minAddBandwidth)
        iBandwidth = std::abs(uiRefSize - uiQuerySize) + minAddBandwidth;
#endif

    Wrapper_ksw_extz_t ez;

    for( int uiCounter = 0; uiCounter < 100; ++uiCounter )
    {
        std::cout << "score=" << ez.ez->score << " ciglen= " << ez.ez->n_cigar << " cigar=" << std::endl;
        // create match/missmatch matrix
        int8_t mat[ 25 ];
        ksw_gen_simple_mat( 5, mat, iMatch, iMissMatch );

        auto vRefSeq = randomNucSeq( uiRefSize ); // nucleotide sequence (codes)
        auto vQuerySeq = randomNucSeq( uiQuerySize ); // nucleotide sequence (codes)

        const char vMIDN[] = {'M', 'I', 'D', 'N'};

        std::cout << "dp" << std::endl;
#if 0
		ksw_extz2_sse(
            nullptr, // used in kmalloc so we do not need this
            vQuerySeq.size(), 
            (const uint8_t*)&vQuerySeq[0],
            vRefSeq.size(), 
            (const uint8_t*)&vRefSeq[0], 
            5, //number of input symbols
            mat, // 5x5 matrix match missmatch
            iGap,
            iExtend,
            iBandwidth,
            -1, // break alignment if score drops to fast (-1 == disabled)
            -1, // bonus(?) for reaching the end of the query (-1 == disabled)
            0, // flags (0 == no flags set)
            ez.ez // return value
        );
#else
        ksw_simplified( vQuerySeq.size( ), (const uint8_t*)&vQuerySeq[ 0 ], vRefSeq.size( ),
                        (const uint8_t*)&vRefSeq[ 0 ], iGap, iExtend, iGap2, iExtend2, iBandwidth,
                        ez.ez // return value
        );
#endif

        std::cout << "score=" << ez.ez->score << " ciglen= " << ez.ez->n_cigar << " cigar=" << std::endl;
        uint32_t qPos = 0;
        uint32_t rPos = 0;
        for( int i = 0; i < ez.ez->n_cigar; ++i )
        {
            uint32_t uiSymb = ez.ez->cigar[ i ] & 0xf;
            uint32_t uiAmount = ez.ez->cigar[ i ] >> 4;
            std::cout << uiAmount << vMIDN[ uiSymb ] << " ";
            if( uiSymb == 0 || uiSymb == 1 )
                qPos += uiAmount;
            if( uiSymb == 0 || uiSymb == 2 )
                rPos += uiAmount;
        } // for
        std::cout << std::endl;
        std::cout << "qPos=" << qPos << " rPos= " << rPos << std::endl;
        assert( qPos == (uint32_t)uiQuerySize );
        assert( rPos == (uint32_t)uiRefSize );
    } // for
} // function

std::string run_ksw( std::string sA, std::string sB, int8_t iM, int8_t iMm, int8_t iO, int8_t iO2, int8_t iE,
                     int8_t iE2, int iW )
{
    // make matrix
    int8_t mat[ 25 ];
    ksw_gen_simple_mat( 5, mat, iM, iMm );

    Wrapper_ksw_extz_t ez;
    std::vector<uint8_t> vA;
    for( size_t i = 0; i < sA.size( ); i++ )
        vA.push_back( sA[ i ] );
    std::vector<uint8_t> vB;
    for( size_t i = 0; i < sB.size( ); i++ )
        vB.push_back( sB[ i ] );
    ksw_simplified( sA.size( ), &vA[ 0 ], sB.size( ), &vB[ 0 ], iO, iE, iO2, iE2, iW,
                    ez.ez, // return value
                    mat // match mismatch matrix
    );

    const char vMIDN[] = {'M', 'I', 'D', 'X'};
    std::string sRet = "";
    for( int i = 0; i < ez.ez->n_cigar; ++i )
    {
        uint32_t uiSymb = ez.ez->cigar[ i ] & 0xf;
        uint32_t uiAmount = ez.ez->cigar[ i ] >> 4;
        sRet += std::to_string( uiAmount );
        sRet += vMIDN[ uiSymb ];
        sRet += ",";
    } // for
    return sRet;
} // function
#endif

#ifdef WITH_PYTHON
void exportNeedlemanWunsch( )
{
    // test ksw function
    // DEBUG( boost::python::def( "testKsw", &testKsw ); ) // DEBUG
    // boost::python::def( "run_ksw", &run_ksw );

    // export the NeedlemanWunsch class
    exportModule<NeedlemanWunsch>( "NeedlemanWunsch" );
} // function
#endif