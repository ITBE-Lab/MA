#pragma once

#include "container/nucSeq.h"
#include "container/seed.h"
#include <cmath>

namespace libMA
{

class SvJump : public Container
{
    static inline nucSeqIndex dist( nucSeqIndex uiA, nucSeqIndex uiB )
    {
        return uiA < uiB ? uiB - uiA : uiA - uiB;
    } // method
    static const nucSeqIndex uiSeedDirFuzziness = 3;

  public:
    static bool validJump( const Seed& rA, const Seed& rB, const bool bFromSeedStart )
    {
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
    const bool bFromSeedStart;
    int64_t iId;
    int64_t iReadId;

    SvJump( const nucSeqIndex uiFrom,
            const nucSeqIndex uiTo,
            const nucSeqIndex uiQueryFrom,
            const nucSeqIndex uiQueryTo,
            const bool bFromForward,
            const bool bToForward,
            const bool bFromSeedStart,
            int64_t iId = -1, /* -1 == no id obtained */
            int64_t iReadId = -1 /* -1 == no id obtained */ )
        : uiFrom( uiFrom ),
          uiTo( uiTo ),
          uiQueryFrom( uiQueryFrom ),
          uiQueryTo( uiQueryTo ),
          bFromForward( bFromForward ),
          bToForward( bToForward ),
          bFromSeedStart( bFromSeedStart ),
          iId( iId ),
          iReadId( iReadId )
    {
        assert( uiQueryFrom <= uiQueryTo );
        // necessary for mapping switch strand jumps rightwards
        assert( uiFrom * 2 + 1000 < static_cast<nucSeqIndex>( std::numeric_limits<int64_t>::max( ) ) );
    } // constructor

    SvJump( const Seed& rA, const Seed& rB, const bool bFromSeedStart )
        : SvJump( /* uiFrom = */
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
                  std::min( bFromSeedStart ? rA.start( ) : rA.end( ) - 1,
                            !bFromSeedStart ? rB.start( ) : rB.end( ) - 1 ),
                  /* uiQueryTo = */
                  std::max( bFromSeedStart ? rA.start( ) : rA.end( ) - 1,
                            !bFromSeedStart ? rB.start( ) : rB.end( ) - 1 ),
                  /* bFromForward = */ rA.bOnForwStrand,
                  /* bToForward = */ rB.bOnForwStrand,
                  /* bFromSeedStart = */ bFromSeedStart )
    {} // constructor

    bool does_switch_strand( ) const
    {
        return bFromForward != bToForward;
    } // method

    bool from_fuzziness_is_rightwards( ) const
    {
        return !does_switch_strand( ) || bFromForward != bFromSeedStart;
    } // method

    nucSeqIndex fuzziness( ) const
    {
        return std::min( static_cast<nucSeqIndex>( 1 + std::pow( std::max( dist( uiFrom, uiTo ), //
                                                                           uiQueryTo - uiQueryFrom ), //
                                                                 1.5 ) //
                                                           / 1000 ),
                         (nucSeqIndex)1000 );
    } // method

    // down == left
    bool to_fuzziness_is_downwards( ) const
    {
        return !does_switch_strand( ) || bToForward != bFromSeedStart;
    } // method

    int64_t from_start_same_strand( ) const
    {
        if( from_fuzziness_is_rightwards( ) )
            return ( (int64_t)uiFrom ) - ( (int64_t)uiSeedDirFuzziness );
        return ( (int64_t)uiFrom ) - ( (int64_t)fuzziness( ) );
    } // method

    int64_t from_start( ) const
    {
        return from_start_same_strand( ) +
               ( does_switch_strand( ) ? std::numeric_limits<int64_t>::max( ) / (int64_t)2 : 0 );
    } // method

    nucSeqIndex from_size( ) const
    {
        return fuzziness( ) + uiSeedDirFuzziness;
    } // method

    int64_t from_end( ) const
    {
        return from_start( ) + (int64_t)from_size( ) - 1;
    } // method

    int64_t to_start( ) const
    {
        if( !to_fuzziness_is_downwards( ) )
            return (int64_t)uiTo - (int64_t)uiSeedDirFuzziness;
        return (int64_t)uiTo - (int64_t)fuzziness( );
    } // method

    nucSeqIndex to_size( ) const
    {
        return fuzziness( ) + uiSeedDirFuzziness;
    } // method

    int64_t to_end( ) const
    {
        return to_start( ) + to_size( ) - 1;
    } // method

    nucSeqIndex query_distance( ) const
    {
        return uiQueryTo - uiQueryFrom;
    } // method

    double score( ) const // @todo really necessary ?
    {
        return 0.08 * std::log( query_distance( ) + 1.5 );
    } // method
}; // class

class SvCall : public Container
{
  public:
    nucSeqIndex uiFromStart;
    nucSeqIndex uiToStart;
    nucSeqIndex uiFromSize;
    nucSeqIndex uiToSize;
    bool bSwitchStrand;
    double dScore;
    std::vector<int64_t> vSupportingJumpIds;
    int64_t iId;

    // these can be empty
    std::shared_ptr<NucSeq> pInsertedSequence;
    std::vector<SvJump> vSupportingJumps;

    SvCall( nucSeqIndex uiFromStart,
            nucSeqIndex uiToStart,
            nucSeqIndex uiFromSize,
            nucSeqIndex uiToSize,
            bool bSwitchStrand,
            double dScore,
            std::vector<int64_t> vSupportingJumpIds = {},
            int64_t iId = -1 /* -1 == no id obtained */ )
        : uiFromStart( uiFromStart ),
          uiToStart( uiToStart ),
          uiFromSize( uiFromSize ),
          uiToSize( uiToSize ),
          bSwitchStrand( bSwitchStrand ),
          dScore( dScore ),
          vSupportingJumpIds( vSupportingJumpIds ),
          iId( iId )
    {} // constructor

    SvCall( const SvJump& rJump, bool bRememberJump = false )
        : SvCall( rJump.from_start_same_strand( ),
                  rJump.to_start( ),
                  rJump.from_size( ),
                  rJump.to_size( ),
                  rJump.does_switch_strand( ),
                  rJump.score( ),
                  std::vector<int64_t>{rJump.iId} )
    {
        if(bRememberJump)
            vSupportingJumps.push_back(rJump);
    } // constructor

    SvCall( SvJump& rJump, bool bRememberJump = false )
        : SvCall( rJump.from_start_same_strand( ),
                  rJump.to_start( ),
                  rJump.from_size( ),
                  rJump.to_size( ),
                  rJump.does_switch_strand( ),
                  rJump.score( ),
                  std::vector<int64_t>{rJump.iId} )
    {
        if(bRememberJump)
            vSupportingJumps.push_back(rJump);
    } // constructor

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
        this->vSupportingJumpIds.insert(
            this->vSupportingJumpIds.end( ), rOther.vSupportingJumpIds.begin( ), rOther.vSupportingJumpIds.end( ) );
        this->dScore += rOther.dScore;
        assert( !this->supportedJumpsLoaded( ) );
        assert( !rOther.supportedJumpsLoaded( ) );
        assert( !this->insertedSequenceComputed( ) );
        assert( !rOther.insertedSequenceComputed( ) );
        assert( !this->hasId( ) );
        assert( !rOther.hasId( ) );
    } // method
}; // class

}; // namespace libMA

#ifdef WITH_PYTHON
void exportSVJump( py::module& rxPyModuleId );
#endif