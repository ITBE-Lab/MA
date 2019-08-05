#pragma once

#include "container/nucSeq.h"
#include "container/seed.h"
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

    const double s, s_neg;
    const double h;
    const double m;
    const nucSeqIndex uiSeedDirFuzziness;
    const nucSeqIndex uiSDFActivate;

  public:
    static bool validJump( const Seed& rA, const Seed& rB, const bool bFromSeedStart )
    {
        // do not create edges between seeds that are overlapping more than 3 nt on the query.
        if( rA.end( ) > rB.start( ) + 3 && rB.end( ) > rA.start( ) + 3 )
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
    {} // constructor

    SvJump( std::shared_ptr<Presetting> pSelectedSetting, const Seed& rA, const nucSeqIndex qLen,
            const bool bFromSeedStart, int64_t iReadId )
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
                  bFromSeedStart ? 0 : rA.end( ) - 1,
                  /* uiQueryTo = */
                  !bFromSeedStart ? qLen : rA.start( ),
                  /* bFromForward = */ rA.bOnForwStrand,
                  /* bToForward = */ rA.bOnForwStrand,
                  /* bFromSeedStart = */ bFromSeedStart,
                  /* uiNumSupportingNt = */ rA.size( ),
                  /* iID */ -1,
                  /* iReadId */ iReadId )
    {} // constructor

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
        double h_min = 1;
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
            return 1;

        return fuzziness( ) + getSeedDirFuzziness( );
    } // method

    int64_t to_end( ) const
    {
        return to_start( ) + to_size( );
    } // method

    int64_t sweep_end( ) const
    {
        return this->switch_strand_known( ) ? to_end( ) - 1 : to_start( ) + from_size( );
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

class SvCall : public Container
{
  public:
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

    nucSeqIndex uiFromStart;
    nucSeqIndex uiToStart;
    nucSeqIndex uiFromSize;
    nucSeqIndex uiToSize;
    bool bSwitchStrand;
    nucSeqIndex uiNumSuppNt;
    nucSeqIndex uiCoverage;
    std::vector<int64_t> vSupportingJumpIds;
    int64_t iId;
    Regex xRegex;
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
            nucSeqIndex uiNumSuppNt,
            std::vector<int64_t> vSupportingJumpIds = {},
            int64_t iId = -1, /* -1 == no id obtained */
            Regex xRegex = Regex( "", 0 ) )
        : uiFromStart( uiFromStart ),
          uiToStart( uiToStart ),
          uiFromSize( uiFromSize ),
          uiToSize( uiToSize ),
          bSwitchStrand( bSwitchStrand ),
          uiNumSuppNt( uiNumSuppNt ),
          uiCoverage( 1 ),
          vSupportingJumpIds( vSupportingJumpIds ),
          iId( iId ),
          xRegex( xRegex )
    {} // constructor

    SvCall( nucSeqIndex uiFromStart,
            nucSeqIndex uiToStart,
            nucSeqIndex uiFromSize,
            nucSeqIndex uiToSize,
            bool bSwitchStrand,
            nucSeqIndex uiNumSuppNt,
            uint32_t uiCoverage )
        : SvCall( uiFromStart, uiToStart, uiFromSize, uiToSize, bSwitchStrand, uiNumSuppNt )
    {
        this->uiCoverage = uiCoverage;
    } // constructor

    SvCall( std::shared_ptr<SvJump> pJump, bool bRememberJump = true )
        : SvCall( pJump->from_start_same_strand( ),
                  pJump->to_start( ),
                  pJump->from_size( ),
                  pJump->to_size( ),
                  pJump->does_switch_strand( ),
                  pJump->numSupportingNt( ),
                  std::vector<int64_t>{pJump->iId} )
    {
        if( bRememberJump )
            vSupportingJumps.push_back( pJump );
        addJumpToEstimateClusterSize( pJump );
        uiOpenEdges = 1;
    } // constructor

    inline void resetEstimateClusterSize()
    {
        vLeft.clear();
        vRight.clear();
        vUp.clear();
        vDown.clear();
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
    static const nucSeqIndex uiT2 = 5;
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
        this->uiFromStart = this->left( );
        this->uiToStart = this->down( );
        // make sure we never have a call smaller than 10x10:
        if( uiRight < 10 + this->uiFromStart )
        {
            this->uiFromStart = this->uiFromStart / 2 + uiRight / 2;
            this->uiFromSize = 10;
        } // if
        else
            this->uiFromSize = uiRight - uiFromStart;
        if( uiUp < 10 + this->uiToStart )
        {
            this->uiToStart = this->uiToStart / 2 + uiUp / 2;
            this->uiToSize = 10;
        } // if
        else
            this->uiToSize = uiUp - uiToStart;
    } // method

    inline nucSeqIndex size( ) const
    {
        assert( this->supportedJumpsLoaded( ) );
        nucSeqIndex uiRefSize = uiFromStart > uiToStart ? uiFromStart - uiToStart : uiToStart - uiFromStart;

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
        nucSeqIndex uiFromEnd =
            std::max( this->uiFromStart + this->uiFromSize, rOther.uiFromStart + rOther.uiFromSize );
        nucSeqIndex uiToEnd = std::max( this->uiToStart + this->uiToSize, rOther.uiToStart + rOther.uiToSize );
        this->uiFromStart = std::min( this->uiFromStart, rOther.uiFromStart );
        this->uiToStart = std::min( this->uiToStart, rOther.uiToStart );
        this->uiFromSize = uiFromEnd - this->uiFromStart;
        this->uiToSize = uiToEnd - this->uiToStart;
        this->vSupportingJumpIds.insert( this->vSupportingJumpIds.end( ), rOther.vSupportingJumpIds.begin( ),
                                         rOther.vSupportingJumpIds.end( ) );
        this->uiNumSuppNt += rOther.uiNumSuppNt;
        this->uiCoverage = std::max( rOther.uiCoverage, this->uiCoverage );
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