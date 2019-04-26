#pragma once

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
        if( bFromSeedStart != !rB.bOnForwStrand ) // cases (0,2) (0,3) (3,0) (3,1)
            return true;
        if( !rA.bOnForwStrand != bFromSeedStart && rB.bOnForwStrand ) // cases (1,2) (2,1)
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

    SvJump( const nucSeqIndex uiFrom,
            const nucSeqIndex uiTo,
            const nucSeqIndex uiQueryFrom,
            const nucSeqIndex uiQueryTo,
            const bool bFromForward,
            const bool bToForward,
            const bool bFromSeedStart )
        : uiFrom( uiFrom ),
          uiTo( uiTo ),
          uiQueryFrom( uiQueryFrom ),
          uiQueryTo( uiQueryTo ),
          bFromForward( bFromForward ),
          bToForward( bToForward ),
          bFromSeedStart( bFromSeedStart )
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
        return !does_switch_strand() || bFromForward == bFromSeedStart;
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
        return !does_switch_strand() || bToForward == bFromSeedStart;
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

}; // namespace libMA

#ifdef WITH_PYTHON
void exportSVJump( py::module& rxPyModuleId );
#endif