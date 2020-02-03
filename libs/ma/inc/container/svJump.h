#pragma once

#include "container/nucSeq.h"
#include "container/seed.h"
#include "geom.h"
#include "util/parameter.h"
#include <cmath>
#include <limits>

namespace libMA
{

class SvJump : public Container
{
    static inline nucSeqIndex dist( nucSeqIndex uiA, nucSeqIndex uiB )
    {
        return uiA < uiB ? uiB - uiA : uiA - uiB;
    } // method

    // @todo this should not be here....
    const double s, s_neg;
    const double h;
    const double m;
    const nucSeqIndex uiSeedDirFuzziness;
    const nucSeqIndex uiSDFActivate;

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

    const nucSeqIndex uiFrom; // inclusive
    const nucSeqIndex uiTo; // inclusive
    const nucSeqIndex uiQueryFrom; // inclusive
    const nucSeqIndex uiQueryTo; // inclusive
    const bool bFromForward;
    const bool bToForward;
    const bool bFromSeedStart; // this should be call seed start or end of first seed
    nucSeqIndex uiNumSupportingNt;
    int64_t iId;
    int64_t iReadId;


#if DEBUG_LEVEL > 0
    size_t uiSeedAId = 0, uiSeedBId = 0;
#endif

    SvJump( std::shared_ptr<Presetting> pSelectedSetting,
            const nucSeqIndex uiFrom,
            const nucSeqIndex uiTo,
            const nucSeqIndex uiQueryFrom,
            const nucSeqIndex uiQueryTo,
            const bool bFromForward,
            const bool bToForward,
            const bool bFromSeedStart,
            const nucSeqIndex uiNumSupportingNt,
            int64_t iId = -1, /* -1 == no id obtained */
            int64_t iReadId = -1 /* -1 == no id obtained */ )
        : s( pSelectedSetting->xJumpS->get( ) ),
          s_neg( pSelectedSetting->xJumpSNeg->get( ) ),
          h( pSelectedSetting->xJumpH->get( ) ),
          m( pSelectedSetting->xJumpM->get( ) ),
          uiSeedDirFuzziness( pSelectedSetting->xSeedDirFuzziness->get( ) ),
          uiSDFActivate( uiSeedDirFuzziness * 2 ),
          uiFrom( uiFrom ),
          uiTo( uiTo ),
          uiQueryFrom( uiQueryFrom ),
          uiQueryTo( uiQueryTo ),
          bFromForward( bFromForward ),
          bToForward( bToForward ),
          bFromSeedStart( bFromSeedStart ),
          uiNumSupportingNt( uiNumSupportingNt ),
          iId( iId ),
          iReadId( iReadId )
    {
        assert( uiQueryFrom <= uiQueryTo );
        // necessary for mapping switch strand jumps rightwards
        assert( uiFrom * 2 + 1000 < static_cast<nucSeqIndex>( std::numeric_limits<int64_t>::max( ) ) );
    } // constructor

    SvJump( std::shared_ptr<Presetting> pSelectedSetting, const Seed& rA, const Seed& rB, const bool bFromSeedStart,
            int64_t iReadId )
        : SvJump(
              pSelectedSetting,
              /* uiFrom = */
              bFromSeedStart
                  ? rA.start_ref( )
                  : ( rA.bOnForwStrand ? rA.end_ref( ) - 1
                                       // @note rA's direction is mirrored on reference if rA is on rev comp strand
                                       : rA.start_ref( ) - rA.size( ) + 1 ),
              /* uiTo = */
              !bFromSeedStart
                  ? rB.start_ref( )
                  : ( rB.bOnForwStrand ? rB.end_ref( ) - 1
                                       // @note rB's direction is mirrored on reference if rB is on rev comp strand
                                       : rB.start_ref( ) - rB.size( ) + 1 ),
              /* uiQueryFrom = */
              std::min( bFromSeedStart ? rA.start( ) : rA.end( ) - 1, !bFromSeedStart ? rB.start( ) : rB.end( ) - 1 ),
              /* uiQueryTo = */
              std::max( bFromSeedStart ? rA.start( ) : rA.end( ) - 1, !bFromSeedStart ? rB.start( ) : rB.end( ) - 1 ),
              /* bFromForward = */ rA.bOnForwStrand,
              /* bToForward = */ rB.bOnForwStrand,
              /* bFromSeedStart = */ bFromSeedStart,
              /* uiNumSupportingNt = */ rA.size( ) + rB.size( ),
              /* iID */ -1,
              /* iReadId */ iReadId )
    {
#if DEBUG_LEVEL > 0
        uiSeedAId = rA.uiId;
        uiSeedBId = rB.uiId;
#endif
    } // constructor

    SvJump( std::shared_ptr<Presetting> pSelectedSetting, const Seed& rA, const nucSeqIndex qLen,
            const bool bFromSeedStart, int64_t iReadId, nucSeqIndex uiMaxJumpLen )
        : SvJump( pSelectedSetting,
                  /* uiFrom = */
                  bFromSeedStart
                      // if we jump to the start of the first seed we don't know where we are coming from
                      ? std::numeric_limits<uint32_t>::max( )
                      : ( rA.bOnForwStrand ? rA.end_ref( ) - 1
                                           // @note rA's direction is mirrored on reference if rA is on rev comp strand
                                           : rA.start_ref( ) - rA.size( ) + 1 ),
                  /* uiTo = */
                  !bFromSeedStart
                      // if we jump from the end of the last seed we don't know where we are going to
                      ? std::numeric_limits<uint32_t>::max( )
                      : rA.start_ref( ),
                  /* uiQueryFrom = */
                  bFromSeedStart ? ( rA.start( ) > uiMaxJumpLen ? rA.start( ) - uiMaxJumpLen : 0 ) : rA.end( ) - 1,
                  /* uiQueryTo = */
                  !bFromSeedStart ? ( rA.end( ) + uiMaxJumpLen < qLen ? rA.end( ) + uiMaxJumpLen : qLen - 1 )
                                  : rA.start( ),
                  /* bFromForward = */ rA.bOnForwStrand,
                  /* bToForward = */ rA.bOnForwStrand,
                  /* bFromSeedStart = */ bFromSeedStart,
                  /* uiNumSupportingNt = */ rA.size( ),
                  /* iID */ -1,
                  /* iReadId */ iReadId )
    {
#if DEBUG_LEVEL > 0
        uiSeedAId = rA.uiId;
        uiSeedBId = rA.uiId;
#endif
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

    bool from_fuzziness_is_rightwards( ) const
    {
        if( !from_known( ) )
            return false;
        if( !to_known( ) )
            return true;
        return !does_switch_strand( ) || bFromForward != bFromSeedStart;
    } // method

    nucSeqIndex fuzziness( ) const
    {
        double x = (double)std::max( dist( uiFrom, uiTo ), uiQueryTo - uiQueryFrom );
        assert( x >= 0 );
        double h_min = 0;
        return (nucSeqIndex)std::min(
            h, h_min + std::max( 0.0, x - ( uiTo >= uiFrom || uiQueryTo - uiQueryFrom >= uiFrom - uiTo ? s : s_neg ) ) *
                           m );
    } // method

    // down == left
    bool to_fuzziness_is_downwards( ) const
    {
        if( !from_known( ) )
            return true;
        if( !to_known( ) )
            return false;
        return !does_switch_strand( ) || bToForward != bFromSeedStart;
    } // method

    nucSeqIndex query_distance( ) const
    {
        return uiQueryTo - uiQueryFrom;
    } // method

    int64_t getSeedDirFuzziness( ) const
    {
        if( !from_known( ) || !to_known( ) )
            return query_distance( ) > uiSDFActivate ? (int64_t)uiSeedDirFuzziness : 0;
        return fuzziness( ) > uiSDFActivate ? (int64_t)uiSeedDirFuzziness : 0;
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

    int64_t from_start( ) const
    {
        return from_start_same_strand( ) +
               ( does_switch_strand( ) ? std::numeric_limits<int64_t>::max( ) / (int64_t)2 : 0 );
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

    int64_t to_start( ) const
    {
        if( !from_known( ) )
            return std::max( (int64_t)0, ( (int64_t)uiTo ) - (int64_t)query_distance( ) + getSeedDirFuzziness( ) );
        if( !to_known( ) )
            return std::max( (int64_t)0, ( (int64_t)uiFrom ) - getSeedDirFuzziness( ) );

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
            return std::numeric_limits<nucSeqIndex>::max( ) / (nucSeqIndex)2;
        return std::max( query_distance( ), ref_distance( ) );
    } // method

    nucSeqIndex numSupportingNt( ) const
    {
        return uiNumSupportingNt;
    } // method

    int64_t insert_ratio( ) const
    {
        if( !switch_strand_known( ) )
            return std::numeric_limits<int64_t>::max( ) / 2;
        return (int64_t)query_distance( ) - (int64_t)ref_distance( );
    } // method
}; // class

class SvCall : public Container, public geom::Rectangle<nucSeqIndex>
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

    bool bSwitchStrand;
    nucSeqIndex uiNumSuppReads;
    nucSeqIndex uiReferenceAmbiguity;
    std::vector<int64_t> vSupportingJumpIds;
    int64_t iId;
#if 0
    Regex xRegex;
#endif
    size_t uiOpenEdges = 0;
    std::vector<nucSeqIndex> vLeft;
    std::vector<nucSeqIndex> vRight;
    std::vector<nucSeqIndex> vUp;
    std::vector<nucSeqIndex> vDown;

    // these can be empty
    std::shared_ptr<NucSeq> pInsertedSequence;
    std::vector<std::shared_ptr<SvJump>> vSupportingJumps;

    SvCall( nucSeqIndex uiFromStart,
            nucSeqIndex uiToStart,
            nucSeqIndex uiFromSize,
            nucSeqIndex uiToSize,
            bool bSwitchStrand,
            nucSeqIndex uiNumSuppReads,
            std::vector<int64_t> vSupportingJumpIds = {},
            int64_t iId = -1 /* -1 == no id obtained */
#if 0
            , Regex xRegex = Regex( "", 0 )
#endif
            )
        : Rectangle( uiFromStart, uiToStart, uiFromSize, uiToSize ),
          bSwitchStrand( bSwitchStrand ),
          uiNumSuppReads( uiNumSuppReads ),
          uiReferenceAmbiguity( 1 ),
          vSupportingJumpIds( vSupportingJumpIds ),
          iId( iId )
#if 0
          ,
          xRegex( xRegex )
#endif
    {} // constructor

    SvCall( nucSeqIndex uiFromStart,
            nucSeqIndex uiToStart,
            nucSeqIndex uiFromSize,
            nucSeqIndex uiToSize,
            bool bSwitchStrand,
            nucSeqIndex uiNumSuppReads,
            uint32_t uiReferenceAmbiguity )
        : SvCall( uiFromStart, uiToStart, uiFromSize, uiToSize, bSwitchStrand, uiNumSuppReads )
    {
        this->uiReferenceAmbiguity = uiReferenceAmbiguity;
    } // constructor

    SvCall( std::shared_ptr<SvJump> pJump, bool bRememberJump = true )
        : SvCall( pJump->from_start_same_strand( ),
                  pJump->to_start( ),
                  pJump->from_size( ),
                  pJump->to_size( ),
                  pJump->does_switch_strand( ),
                  1,
                  std::vector<int64_t>{pJump->iId} )
    {
        if( bRememberJump )
            vSupportingJumps.push_back( pJump );
        addJumpToEstimateClusterSize( pJump );
        uiOpenEdges = 1;
    } // constructor

    inline void resetEstimateClusterSize( )
    {
        vLeft.clear( );
        vRight.clear( );
        vUp.clear( );
        vDown.clear( );
    } // method

    inline void addJumpToEstimateClusterSize( std::shared_ptr<SvJump> pJump )
    {
        nucSeqIndex uiF = pJump->uiFrom;
        nucSeqIndex uiT = pJump->uiTo;
        if( !pJump->from_known( ) )
            uiF = uiT;
        else if( !pJump->to_known( ) )
            uiT = uiF;
        nucSeqIndex uiDelToInsRatio = 1; // @todo check if 1 is a well enough approximation
        if( pJump->from_fuzziness_is_rightwards( ) )
        {
            this->vLeft.push_back( uiF );
            this->vRight.push_back( uiF + std::min( pJump->ref_distance( ) / uiDelToInsRatio,
                                                    pJump->query_distance( ) * uiDelToInsRatio ) );
        } // if
        else
        {
            this->vRight.push_back( uiF );
            nucSeqIndex uiX =
                std::min( pJump->ref_distance( ) / uiDelToInsRatio, pJump->query_distance( ) * uiDelToInsRatio );
            this->vLeft.push_back( uiF > uiX ? uiF - uiX : 0 );
        } // else
        if( pJump->to_fuzziness_is_downwards( ) )
        {
            this->vUp.push_back( uiT );
            nucSeqIndex uiX =
                std::min( pJump->ref_distance( ) / uiDelToInsRatio, pJump->query_distance( ) * uiDelToInsRatio );
            this->vDown.push_back( uiT > uiX ? uiT - uiX : 0 );
        } // if
        else
        {
            this->vDown.push_back( uiT );
            this->vUp.push_back( uiT + std::min( pJump->ref_distance( ) / uiDelToInsRatio,
                                                 pJump->query_distance( ) * uiDelToInsRatio ) );
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
        vSupportingJumpIds.push_back( pJmp->iId );
        vSupportingJumps.push_back( pJmp );
    } // method

    inline size_t estimateCoverage( )
    {
        return vSupportingJumpIds.size( );
    } // method

    static const size_t uiT = 5;
    static const nucSeqIndex uiT2 = 0;
    inline nucSeqIndex right( )
    {
        std::sort( vRight.begin( ), vRight.end( ) );
        return vRight[ vRight.size( ) / uiT ] + uiT2;
    } // method

    inline nucSeqIndex left( )
    {
        std::sort( vLeft.begin( ), vLeft.end( ), []( nucSeqIndex uiA, nucSeqIndex uiB ) { return uiA > uiB; } );
        nucSeqIndex uiX = vLeft[ vLeft.size( ) / uiT ];
        return uiX > uiT2 ? uiX - uiT2 : 0;
    } // method
    inline nucSeqIndex up( )
    {
        std::sort( vUp.begin( ), vUp.end( ) );
        return vUp[ vUp.size( ) / uiT ] + uiT2;
    } // method

    inline nucSeqIndex down( )
    {
        std::sort( vDown.begin( ), vDown.end( ), []( nucSeqIndex uiA, nucSeqIndex uiB ) { return uiA > uiB; } );
        nucSeqIndex uiX = vDown[ vDown.size( ) / uiT ];
        return uiX > uiT2 ? uiX - uiT2 : 0;
    } // method

    inline void reEstimateClusterSize( )
    {
        nucSeqIndex uiRight = this->right( );
        nucSeqIndex uiUp = this->up( );
        this->xXAxis.start( this->left( ) );
        this->xYAxis.start( this->down( ) );
        if( uiRight < this->xXAxis.start( ) )
        {
            this->xXAxis.start( ( this->xXAxis.start( ) + uiRight ) / 2 - 3 );
            this->xXAxis.size( 5 );
        } // if
        else
            this->xXAxis.size( uiRight - this->xXAxis.start( ) );
        if( uiUp < this->xYAxis.start( ) )
        {
            this->xYAxis.start( ( this->xYAxis.start( ) + uiUp ) / 2 - 3 );
            this->xYAxis.size( 5 );
        } // if
        else
            this->xYAxis.size( uiUp - this->xYAxis.start( ) );
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
        assert( this->bSwitchStrand == rOther.bSwitchStrand );
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

        this->vLeft.insert( this->vLeft.end( ), rOther.vLeft.begin( ), rOther.vLeft.end( ) );
        this->vRight.insert( this->vRight.end( ), rOther.vRight.begin( ), rOther.vRight.end( ) );
        this->vUp.insert( this->vUp.end( ), rOther.vUp.begin( ), rOther.vUp.end( ) );
        this->vDown.insert( this->vDown.end( ), rOther.vDown.begin( ), rOther.vDown.end( ) );
    } // method
}; // class

}; // namespace libMA

#ifdef WITH_PYTHON
void exportSVJump( py::module& rxPyModuleId );
#endif