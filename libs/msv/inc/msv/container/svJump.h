#pragma once

#include "ma/container/nucSeq.h"
#include "ma/container/seed.h"
#include "ms/util/parameter.h"
#include "msv/util/statisticSequenceAnalysis.h"
#include "util/geom.h"
#include <cmath>
#include <limits>

using namespace libMA;
namespace libMSV
{

class SvJump : public libMS::Container
{
    static inline nucSeqIndex dist( nucSeqIndex uiA, nucSeqIndex uiB )
    {
        return uiA < uiB ? uiB - uiA : uiA - uiB;
    } // method

    nucSeqIndex sDFActivate( ) const
    {
        return (nucSeqIndex)pGlobalParams->xSeedDirFuzziness->get( ) * 2;
    } // method
  public:
    static bool validJump( const Seed& rA, const Seed& rB, const bool bFromSeedStart )
    {
        // do not create edges between seeds that are overlapping more than 5 nt on the query.
        if( rA.end( ) > rB.start( ) + 5 && rB.end( ) > rA.start( ) + 5 )
            return false;
        if( bFromSeedStart != rB.bOnForwStrand ) // cases (0,2) (0,3) (3,0) (3,1)
            return true;
        if( !rA.bOnForwStrand && bFromSeedStart && rB.bOnForwStrand ) // case (1,2)
            return true;
        if( rA.bOnForwStrand && !bFromSeedStart && !rB.bOnForwStrand ) // case (2,1)
            return true;
        return false;
    } // method

    /*const*/ bool bWasMirrored; // this should be called was_mirrored
    /*const*/ nucSeqIndex uiFrom; // inclusive
    /*const*/ nucSeqIndex uiTo; // inclusive
    /*const*/ nucSeqIndex uiQueryFrom; // inclusive
    /*const*/ nucSeqIndex uiQueryTo; // inclusive
    /*const*/ bool bFromForward;
    /*const*/ bool bToForward;
    nucSeqIndex uiNumSupportingNt;
    int64_t iId;
    int64_t iReadId;


#if DEBUG_LEVEL > 0
    size_t uiSeedAId = 0, uiSeedBId = 0;
#endif

    SvJump( const SvJump& rOther )
        : bWasMirrored( rOther.bWasMirrored ),
          uiFrom( rOther.uiFrom ),
          uiTo( rOther.uiTo ),
          uiQueryFrom( rOther.uiQueryFrom ),
          uiQueryTo( rOther.uiQueryTo ),
          bFromForward( rOther.bFromForward ),
          bToForward( rOther.bToForward ),
          uiNumSupportingNt( rOther.uiNumSupportingNt ),
          iId( rOther.iId ),
          iReadId( rOther.iReadId )
    {}

    /**
     * @brief creates a jump from completely given data
     * @details
     * Used for constructing jump objects with data from the DB
     */
    SvJump( const nucSeqIndex uiFrom,
            const nucSeqIndex uiTo,
            const nucSeqIndex uiQueryFrom,
            const nucSeqIndex uiQueryTo,
            const bool bFromForward,
            const bool bToForward,
            const bool bWasMirrored,
            const nucSeqIndex uiNumSupportingNt,
            int64_t iId = -1, /* -1 == no id obtained */
            int64_t iReadId = -1 /* -1 == no id obtained */ )
        : bWasMirrored( bWasMirrored ),
          uiFrom( uiFrom ),
          uiTo( uiTo ),
          uiQueryFrom( uiQueryFrom ),
          uiQueryTo( uiQueryTo ),
          bFromForward( bFromForward ),
          bToForward( bToForward ),
          uiNumSupportingNt( uiNumSupportingNt ),
          iId( iId ),
          iReadId( iReadId )
    {
        assert( uiQueryFrom <= uiQueryTo );
        // necessary for mapping switch strand jumps rightwards
        assert( uiFrom * 8 + 1000 < static_cast<nucSeqIndex>( std::numeric_limits<int64_t>::max( ) ) );
    } // constructor

    SvJump( const nucSeqIndex uiFrom_,
            const nucSeqIndex uiTo_,
            const nucSeqIndex uiQueryFrom,
            const nucSeqIndex uiQueryTo,
            const bool bFromForward,
            const bool bToForward,
            const nucSeqIndex uiNumSupportingNt,
            int64_t iId,
            int64_t iReadId )
        : bWasMirrored( ( uiTo_ < uiFrom_ || ( uiTo_ == uiFrom_ && !bFromForward && bToForward ) ) &&
                        uiFrom_ != std::numeric_limits<uint32_t>::max( ) ),
          uiFrom( bWasMirrored ? uiTo_ : uiFrom_ ),
          uiTo( bWasMirrored ? uiFrom_ : uiTo_ ),
          uiQueryFrom( uiQueryFrom ),
          uiQueryTo( uiQueryTo ),
          bFromForward( bWasMirrored ? !bToForward : bFromForward ),
          bToForward( bWasMirrored ? !bFromForward : bToForward ),
          uiNumSupportingNt( uiNumSupportingNt ),
          iId( iId ),
          iReadId( iReadId )
    {
        assert( uiQueryFrom <= uiQueryTo );
        // necessary for mapping switch strand jumps rightwards
        assert( uiFrom * 8 + 1000 < static_cast<nucSeqIndex>( std::numeric_limits<int64_t>::max( ) ) );

        assert( ( (uint32_t)uiFrom ) != std::numeric_limits<uint32_t>::max( ) ||
                ( (uint32_t)uiTo ) != std::numeric_limits<uint32_t>::max( ) );
    } // constructor


    /**
     * @brief creates a jump between the given two seeds
     * @details
     * rB must occur after rA on the query
     */
    SvJump( const Seed& rA, const Seed& rB, int64_t iReadId )
        : SvJump(
              /* uiFrom = */
              rA.bOnForwStrand ? rA.end_ref( ) - 1
                               // @note rA's direction is mirrored on reference if rA is on rev comp strand
                               : rA.start_ref( ) - rA.size( ) + 1,
              /* uiTo = */
              rB.bOnForwStrand ? rB.start_ref( ) : rB.start_ref( ),
              /* uiQueryFrom = */
              std::min( rA.end( ) - 1, rB.start( ) ),
              /* uiQueryTo = */
              std::max( rA.end( ) - 1, rB.start( ) ),
              /* bFromForward = */ rA.bOnForwStrand,
              /* bToForward = */ rB.bOnForwStrand,
              /* uiNumSupportingNt = */ rA.size( ) + rB.size( ),
              /* iID */ -1,
              /* iReadId */ iReadId )
    {
#if 0
        std::string vS[] = {"false", "true"};
        std::cout << "jump -- " << rA.start( ) << ", " << rA.start_ref( ) << ", " << rA.size( ) << " -> " << rB.start( )
                  << ", " << rB.start_ref( ) << ", " << rB.size( ) << ": " << uiFrom << " - " << uiTo << ", "
                  << vS[ bFromForward ] << ", " << vS[ bToForward ] << std::endl;
#endif
#if DEBUG_LEVEL > 0
        uiSeedAId = rA.uiId;
        uiSeedBId = rB.uiId;
#endif
    } // constructor

    /**
     * @brief creates a dummy jump
     * @details
     */
    SvJump( const Seed& rA, const nucSeqIndex qLen, const bool bFirstSeed, int64_t iReadId, nucSeqIndex uiMaxJumpLen )
        : SvJump( /* uiFrom = */
                  bFirstSeed == rA.bOnForwStrand
                      // if we jump to the start of the first seed we don't know where we are coming from
                      ? std::numeric_limits<uint32_t>::max( )
                      : ( rA.bOnForwStrand ? rA.end_ref( ) - 1
                                           // @note rA's direction is mirrored on reference if rA is on rev comp strand
                                           : 1 + rA.start_ref( ) - rA.size( ) ),
                  /* uiTo = */
                  bFirstSeed != rA.bOnForwStrand
                      // if we jump from the end of the last seed we don't know where we are going to
                      ? std::numeric_limits<uint32_t>::max( )
                      : rA.start_ref( ),
                  /* uiQueryFrom = */
                  bFirstSeed ? ( rA.start( ) > uiMaxJumpLen ? rA.start( ) - uiMaxJumpLen : 0 ) : rA.end( ) - 1,
                  /* uiQueryTo = */
                  !bFirstSeed ? ( rA.end( ) + uiMaxJumpLen < qLen ? rA.end( ) + uiMaxJumpLen : qLen - 1 ) : rA.start( ),
                  /* bFromForward = */ true,
                  /* bToForward = */ true,
                  /* uiNumSupportingNt = */ rA.size( ),
                  /* iID */ -1,
                  /* iReadId */ iReadId )
    {
#if DEBUG_LEVEL > 0
        uiSeedAId = rA.uiId;
        uiSeedBId = rA.uiId;
#endif

        assert( ( (uint32_t)uiFrom ) != std::numeric_limits<uint32_t>::max( ) ||
                ( (uint32_t)uiTo ) != std::numeric_limits<uint32_t>::max( ) );
    } // constructor

    bool does_switch_strand( ) const
    {
        return bFromForward != bToForward;
    } // method

    bool from_known( ) const
    {
        return uiFrom != std::numeric_limits<uint32_t>::max( );
    } // method

    bool to_known( ) const
    {
        return uiTo != std::numeric_limits<uint32_t>::max( );
    } // method

    bool switch_strand_known( ) const
    {
        return from_known( ) && to_known( );
    } // method

    bool is_dummy( ) const
    {
        return !from_known( ) || !to_known( );
    } // method

    bool from_fuzziness_is_rightwards( ) const
    {
        if( !from_known( ) )
            return false;
        if( !to_known( ) )
            return true;
        return bFromForward;
    } // method

    nucSeqIndex fuzziness( ) const
    {
        double x = (double)std::max( dist( uiFrom, uiTo ), uiQueryTo - uiQueryFrom );
        assert( x >= 0 );
        double h_min = 0;
        return (nucSeqIndex)std::min(
            pGlobalParams->xJumpH->get( ),
            h_min + std::max( 0.0, x - ( uiTo >= uiFrom || uiQueryTo - uiQueryFrom >= uiFrom - uiTo
                                             ? pGlobalParams->xJumpS->get( )
                                             : pGlobalParams->xJumpSNeg->get( ) ) ) *
                        pGlobalParams->xJumpM->get( ) );
    } // method

    // down == left
    bool to_fuzziness_is_downwards( ) const
    {
        if( !from_known( ) )
            return true;
        if( !to_known( ) )
            return false;
        return bToForward;
    } // method

    nucSeqIndex query_distance( ) const
    {
        return uiQueryTo - uiQueryFrom;
    } // method

    int64_t getSeedDirFuzziness( ) const
    {
        if( !from_known( ) || !to_known( ) )
            return query_distance( ) > sDFActivate( ) ? (int64_t)pGlobalParams->xSeedDirFuzziness->get( ) : 0;
        return fuzziness( ) > sDFActivate( ) ? (int64_t)pGlobalParams->xSeedDirFuzziness->get( ) : 0;
    }

    int64_t from_start_same_strand( ) const
    {
        if( !from_known( ) )
            return std::max( (int64_t)0, ( (int64_t)uiTo ) - (int64_t)query_distance( ) + getSeedDirFuzziness( ) );
        if( !to_known( ) )
            return std::max( (int64_t)0, ( (int64_t)uiFrom ) - getSeedDirFuzziness( ) );

        if( from_fuzziness_is_rightwards( ) )
            return std::max( (int64_t)0, ( (int64_t)uiFrom ) - getSeedDirFuzziness( ) );
        return std::max( (int64_t)0, ( (int64_t)uiFrom ) - ( (int64_t)fuzziness( ) ) );
    } // method

    static const int64_t FROM_POS_NUM_SECTIONS = 8;
    static const int64_t FROM_POS_NUM_USED_SECTIONS = 8;
    int64_t from_start( ) const
    {
        auto iRet = from_start_same_strand( );
        // separate the 5 different types on the x axis
        // for that we separate the dataspace into 8 sections:
        // | 0               | 1               | 2               | 3               | 4     | 5      | 6      | 7      |
        // | forward-forward | forward-reverse | reverse-forward | reverse-reverse | dummy | unused | unused | unused |
        if( is_dummy( ) )
            return iRet + std::numeric_limits<int64_t>::max( ) / ( FROM_POS_NUM_SECTIONS / 4 );
        if( !bFromForward )
            iRet += std::numeric_limits<int64_t>::max( ) / ( FROM_POS_NUM_SECTIONS / 2 );
        if( !bToForward )
            iRet += std::numeric_limits<int64_t>::max( ) / FROM_POS_NUM_SECTIONS;
        return iRet;
    } // method

    nucSeqIndex from_size( ) const
    {
        if( !to_known( ) || !from_known( ) )
            return query_distance( ) + getSeedDirFuzziness( );

        return fuzziness( ) + getSeedDirFuzziness( );
    } // method

    int64_t from_end( ) const
    {
        return from_start( ) + (int64_t)from_size( );
    } // method

    int64_t from_end_same_strand( ) const
    {
        return from_start_same_strand( ) + (int64_t)from_size( );
    } // method

    int64_t to_start( ) const
    {
        if( !from_known( ) )
            return std::max( (int64_t)0, ( (int64_t)uiTo ) - (int64_t)query_distance( ) + getSeedDirFuzziness( ) ) + 1;
        if( !to_known( ) )
            return std::max( (int64_t)0, ( (int64_t)uiFrom ) - getSeedDirFuzziness( ) ) + 1;

        if( !to_fuzziness_is_downwards( ) )
            return std::max( (int64_t)0, ( (int64_t)uiTo ) - (int64_t)getSeedDirFuzziness( ) );
        return std::max( (int64_t)0, ( (int64_t)uiTo ) - (int64_t)fuzziness( ) );
    } // method

    nucSeqIndex to_size( ) const
    {
        if( !to_known( ) || !from_known( ) )
            return 0;

        return fuzziness( ) + getSeedDirFuzziness( );
    } // method

    int64_t to_end( ) const
    {
        return to_start( ) + to_size( );
    } // method

    int64_t sweep_end( ) const
    {
        return this->switch_strand_known( ) ? to_end( ) : to_start( ) + from_size( );
    } // method

    nucSeqIndex ref_distance( ) const
    {
        return uiTo < uiFrom ? uiFrom - uiTo : uiTo - uiFrom;
    } // method

    nucSeqIndex size( ) const
    {
        if( !from_known( ) || !to_known( ) )
            return std::numeric_limits<nucSeqIndex>::max( ) / (nucSeqIndex)4;
        return std::max( query_distance( ), ref_distance( ) );
    } // method

    nucSeqIndex numSupportingNt( ) const
    {
        return uiNumSupportingNt;
    } // method

    int64_t insert_ratio( ) const
    {
        if( !switch_strand_known( ) )
            return std::numeric_limits<int64_t>::max( ) / 4;
        return (int64_t)query_distance( ) - (int64_t)ref_distance( );
    } // method

    nucSeqIndex referenceAmbiguity( nucSeqIndex uiDistanceMin, nucSeqIndex uiDistanceMax, std::shared_ptr<Pack> pPack )
    {
        return sampleSequenceAmbiguity( uiFrom, uiTo, bFromForward, bToForward, pPack, uiDistanceMin, uiDistanceMax );
    }
}; // class

inline std::ostream& operator<<( std::ostream& os, const SvJump& rJ )
{
    os << "bWasMirrored= " << rJ.bWasMirrored << " uiFrom= " << rJ.uiFrom << " uiTo= " << rJ.uiTo
       << " uiQueryFrom= " << rJ.uiQueryFrom << " uiQueryTo= " << rJ.uiQueryTo << " bFromForward= " << rJ.bFromForward
       << " bToForward= " << rJ.bToForward << " uiNumSupportingNt= " << rJ.uiNumSupportingNt << " iId= " << rJ.iId
       << " iReadId= " << rJ.iReadId;
    return os;
}

// @todo move this to it's own file
class SvCall : public libMS::Container, public geom::Rectangle<nucSeqIndex>
{
  public:
#if 0
    class Regex
    {
      public:
        std::string sRegex;
        uint32_t uiState;
        int64_t iId;

        std::string unrollSimplest( size_t uiS, size_t uiE )
        {
            std::string sRet = "";
            for( size_t uiI = uiS; uiI < uiE; uiI++ )
            {
                if( sRegex[ uiI ] == '+' )
                    continue; // ignore 1 to many loop symbols (e.i. print the symbol once...)
                if( sRegex[ uiI ] == '|' )
                    break; // if we are in a choice only output the first element
                std::string sApp;
                if( sRegex[ uiI ] == '(' )
                {
                    // search the end of the brackets:
                    size_t uiC = 1;
                    size_t uiJ = uiI + 1;
                    while( uiC > 0 )
                    {
                        assert( uiJ < uiE );
                        if( sRegex[ uiJ ] == '(' )
                            uiC += 1;
                        if( sRegex[ uiJ ] == ')' )
                            uiC -= 1;
                        uiJ++;
                    } // while
                    sApp = unrollSimplest( uiI + 1, uiJ );
                    uiI = uiJ;
                } // if
                else
                    sApp.push_back( sRegex[ uiI ] );
                if( uiI + 1 < uiE && ( sRegex[ uiI + 1 ] == '*' || sRegex[ uiI + 1 ] == '?' ) )
                    continue; // delete symbols of 0 to many loops and zero or one occurrences
                sRet.append( sApp );
            } // for
            return sRet;
        } // method

        std::string unrollSimplest( )
        {
            return unrollSimplest( 0, sRegex.size( ) );
        } // method

        Regex( std::string sRegex, uint32_t uiState, int64_t iId = -1 /* -1 == no id obtained */ )
            : sRegex( sRegex ), uiState( uiState ), iId( iId )
        {}

        Regex( int64_t iCallId ) : Regex( std::to_string( iCallId ).append( "+" ), 0 )
        {}

    }; // class
#endif

    bool bFromForward;
    bool bToForward;
    nucSeqIndex uiNumSuppReads;
    nucSeqIndex uiReferenceAmbiguity;
    std::vector<int64_t> vSupportingJumpIds;
    int64_t iId;
    int64_t iOrderID = -1;
    bool bMirrored = false;
    bool bDummy = false;
#if 0
    Regex xRegex;
#endif
    size_t uiOpenEdges = 0;
    /**
     * @brief used for statistical position estimation of SV calls
     * @details
     * Contains the horizontal component for two-sided calls
     * and the minimum compunent for one-sided calls.
     */
    std::vector<nucSeqIndex> vHorizontal;
    /**
     * @brief used for statistical position estimation of SV calls
     * @details
     * Contains the vertical component for two-sided calls
     * and the maximum compunent for one-sided calls.
     */
    std::vector<nucSeqIndex> vVertical;

    // these can be empty
    std::shared_ptr<NucSeq> pInsertedSequence;
    std::vector<std::shared_ptr<SvJump>> vSupportingJumps;

#if 0
    SvCall( const SvCall& rOther )
        : geom::Rectangle<nucSeqIndex>(rOther)
          bFromForward( rOther.bFromForward ),
          bToForward( rOther.bToForward ),
          uiNumSuppReads( rOther.uiNumSuppReads ),
          uiReferenceAmbiguity( rOther.uiReferenceAmbiguity ),
          vSupportingJumpIds( rOther.vSupportingJumpIds ),
          iId( rOther.iId ),
          iOrderID( rOther.iOrderID ),
          bMirrored( rOther.bMirrored ),
          bDummy( rOther.bDummy ),
          uiOpenEdges( rOther.uiOpenEdges ),
          vHorizontal( rOther.vHorizontal ),
          vVertical( rOther.vVertical ),
          pInsertedSequence(
              rOther.pInsertedSequence == nullptr ? nullptr : std::make_shared<NucSeq>( *rOther.pInsertedSequence ) )
    {
        for( auto pJump : rOther.vSupportingJumps )
            vSupportingJumps.push_back( std::make_shared<SvJump>( *pJump ) );
    }
#endif

    SvCall( nucSeqIndex uiFromStart,
            nucSeqIndex uiToStart,
            nucSeqIndex uiFromSize,
            nucSeqIndex uiToSize,
            bool bFromForward,
            bool bToForward,
            nucSeqIndex uiNumSuppReads,
            std::vector<int64_t> vSupportingJumpIds = { },
            int64_t iId = -1, /* -1 == no id obtained */
            bool bMirrored = false,
            bool bDummy = false
#if 0
            , Regex xRegex = Regex( "", 0 )
#endif
            )
        : Rectangle( uiFromStart, uiToStart, uiFromSize, uiToSize ),
          bFromForward( bFromForward ),
          bToForward( bToForward ),
          uiNumSuppReads( uiNumSuppReads ),
          uiReferenceAmbiguity( 1 ),
          vSupportingJumpIds( vSupportingJumpIds ),
          iId( iId ),
          bMirrored( bMirrored ),
          bDummy( bDummy )
#if 0
          ,
          xRegex( xRegex )
#endif
    {} // constructor

    SvCall( nucSeqIndex uiFromStart,
            nucSeqIndex uiToStart,
            nucSeqIndex uiFromSize,
            nucSeqIndex uiToSize,
            bool bFromForward,
            bool bToForward,
            nucSeqIndex uiNumSuppReads,
            uint32_t uiReferenceAmbiguity )
        : SvCall( uiFromStart, uiToStart, uiFromSize, uiToSize, bFromForward, bToForward, uiNumSuppReads )
    {
        this->uiReferenceAmbiguity = uiReferenceAmbiguity;
    } // constructor

    SvCall( std::shared_ptr<SvJump> pJump, bool bRememberJump = true )
        : SvCall( pJump->from_start_same_strand( ),
                  pJump->to_start( ),
                  pJump->from_size( ),
                  pJump->to_size( ),
                  pJump->bFromForward,
                  pJump->bToForward,
                  1,
                  std::vector<int64_t>{ pJump->iId },
                  -1,
                  pJump->bWasMirrored,
                  !pJump->from_known( ) || !pJump->to_known( ) )
    {
        if( bRememberJump )
            vSupportingJumps.push_back( pJump );
        addJumpToEstimateClusterSize( pJump );
        uiOpenEdges = 1;
    } // constructor

    inline void resetEstimateClusterSize( )
    {
        vHorizontal.clear( );
        vVertical.clear( );
    } // method

    inline void addJumpToEstimateClusterSize( std::shared_ptr<SvJump> pJump )
    {
        if( !pJump->from_known( ) )
            // if the from position of a jump is unknown it indicates the maximal position of a one-sided call.
            vVertical.push_back( pJump->uiTo );
        else if( !pJump->to_known( ) )
            // if the to position of a jump is unknown it indicates the minimal position of a one-sided call.
            vHorizontal.push_back( pJump->uiFrom );
        else
        {
            // if both positions are known, we have a two-sided call, where the from pos indicates the horizontal
            // position and the to pos the vertical one
            vHorizontal.push_back( pJump->uiFrom );
            vVertical.push_back( pJump->uiTo );
        } // else
    } // method

    inline double getScore( ) const
    {
        return uiReferenceAmbiguity == 0 ? 0 : uiNumSuppReads / (double)uiReferenceAmbiguity;
    } // method

    bool supportedJumpsLoaded( ) const
    {
        return vSupportingJumps.size( ) == vSupportingJumpIds.size( );
    } // method

    bool insertedSequenceComputed( ) const
    {
        return pInsertedSequence != nullptr;
    } // method

    bool hasId( ) const
    {
        return iId != -1;
    } // method

#if 0
    Regex getDefaultRegex( )
    {
        assert( this->hasId( ) );
        return Regex( this->iId );
    } // method

    bool hasRegex( ) const
    {
        return !xRegex.sRegex.empty( );
    } // method

    Regex getRegex( )
    {
        if( !hasRegex( ) )
            return getDefaultRegex( );
        return xRegex;
    } // method
#endif

    void clear_jumps( )
    {
        vSupportingJumpIds.clear( );
        vSupportingJumps.clear( );
    } // method

    std::shared_ptr<SvJump> get_jump( size_t uiI )
    {
        return vSupportingJumps[ uiI ];
    } // method

    void add_jump( std::shared_ptr<SvJump> pJmp )
    {
        assert( ( !pJump->from_known( ) || !pJump->to_known( ) ) == this->bDummy );
        vSupportingJumpIds.push_back( pJmp->iId );
        vSupportingJumps.push_back( pJmp );
    } // method

    inline size_t estimateCoverage( )
    {
        return vSupportingJumpIds.size( );
    } // method


    inline void reEstimateClusterSize( )
    {
        std::sort( vHorizontal.begin( ), vHorizontal.end( ) );
        std::sort( vVertical.begin( ), vVertical.end( ) );

        if( bDummy )
        {
            nucSeqIndex uiMin = 1;
            nucSeqIndex uiMax = 0;
            size_t uiI = 0;
            size_t uiJ = vVertical.size( );
            while( uiMin > uiMax && uiI < vHorizontal.size( ) && uiJ > 0 )
            {
                uiMin = vHorizontal[ uiI ];
                uiMax = vVertical[ uiJ - 1 ];
                uiI++;
                uiJ--;
            } // while
            nucSeqIndex uiPos;
            if( uiI == vHorizontal.size( ) || uiJ == 0 )
            {
                if( uiI == vHorizontal.size( ) )
                    uiPos = vVertical[ vVertical.size( ) * 0.05 ];
                else
                    uiPos = vHorizontal[ vHorizontal.size( ) * 0.95 ];
            } // if
            else
                uiPos = ( uiMin + uiMax ) / 2;
            this->xXAxis.start( uiPos );
            this->xYAxis.start( uiPos );
        } // if
        else
        {
            this->xXAxis.start( vHorizontal[ vHorizontal.size( ) * ( bFromForward ? 0.95 : 0.05 ) ] );
            this->xYAxis.start( vVertical[ vVertical.size( ) * ( bToForward ? 0.05 : 0.95 ) ] );
        } // else


        this->xXAxis.size( 1 );
        this->xYAxis.size( 1 );
    } // method

    inline nucSeqIndex size( ) const
    {
        assert( this->supportedJumpsLoaded( ) );
        nucSeqIndex uiRefSize = this->xXAxis.start( ) > this->xYAxis.start( )
                                    ? this->xXAxis.start( ) - this->xYAxis.start( )
                                    : this->xYAxis.start( ) - this->xXAxis.start( );

        nucSeqIndex uiQuerySize = 0;
        for( auto pJump : vSupportingJumps )
            uiQuerySize += std::max( pJump->query_distance( ), pJump->ref_distance( ) );
        uiQuerySize /= vSupportingJumps.size( );

        return uiRefSize > uiQuerySize ? uiRefSize : uiQuerySize;
    } // method

    /**
     * joins two sv calls together,
     * in order to do so, they cannot have an ID in the database or the pInsertedSequence computed!
     */
    void join( SvCall& rOther )
    {
        assert( this->bDummy == rOther.bDummy );
        assert( this->bDummy || this->bFromForward == rOther.bFromForward );
        assert( this->bDummy || this->bToForward == rOther.bToForward );
        nucSeqIndex uiFromEnd = std::max( this->xXAxis.end( ), rOther.xXAxis.end( ) );
        nucSeqIndex uiToEnd = std::max( this->xYAxis.end( ), rOther.xYAxis.end( ) );
        this->xXAxis.start( std::min( this->xXAxis.start( ), rOther.xXAxis.start( ) ) );
        this->xYAxis.start( std::min( this->xYAxis.start( ), rOther.xYAxis.start( ) ) );
        this->xXAxis.size( uiFromEnd - this->xXAxis.start( ) );
        this->xYAxis.size( uiToEnd - this->xYAxis.start( ) );
        this->vSupportingJumpIds.insert( this->vSupportingJumpIds.end( ), rOther.vSupportingJumpIds.begin( ),
                                         rOther.vSupportingJumpIds.end( ) );
        this->uiNumSuppReads += rOther.uiNumSuppReads;
        this->uiReferenceAmbiguity = std::max( rOther.uiReferenceAmbiguity, this->uiReferenceAmbiguity );
        this->vSupportingJumps.insert( this->vSupportingJumps.end( ), rOther.vSupportingJumps.begin( ),
                                       rOther.vSupportingJumps.end( ) );
        this->uiOpenEdges += rOther.uiOpenEdges;
        assert( this->supportedJumpsLoaded( ) == rOther.supportedJumpsLoaded( ) );
        assert( !this->insertedSequenceComputed( ) );
        assert( !rOther.insertedSequenceComputed( ) );

        this->vHorizontal.insert( this->vHorizontal.end( ), rOther.vHorizontal.begin( ), rOther.vHorizontal.end( ) );
        this->vVertical.insert( this->vVertical.end( ), rOther.vVertical.begin( ), rOther.vVertical.end( ) );
    } // method

    size_t numJumps( ) const
    {
        return vSupportingJumpIds.size( );
    }
}; // class


class CompleteBipartiteSubgraphClusterVector : public libMS::Container
{
  public:
    std::vector<std::shared_ptr<SvCall>> vContent;
}; // class

}; // namespace libMSV

#ifdef WITH_PYTHON
void exportSVJump( libMS::SubmoduleOrganizer& xOrganizer );
#endif