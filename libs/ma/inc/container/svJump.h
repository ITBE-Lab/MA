#pragma once

#include "container/seed.h"

namespace libMA
{

class SvJump: public Container
{
  public:
    enum SeedOrientation
    {
        forwardToForward,
        reverseToReverse,
        forwardToReverse,
        reverseToForward
    }; // enum

  public:
    const nucSeqIndex uiFrom;
    const nucSeqIndex uiTo;
    const nucSeqIndex uiQueryDistance;
    const SeedOrientation xSeedOrientation;

    SvJump( const Seed& rA, const Seed& rB )
        : uiFrom( rA.end_ref( ) ),
          uiTo( rB.start_ref( ) ),
          uiQueryDistance( rB.start( ) - rA.end( ) ),
          xSeedOrientation(
              rA.bOnForwStrand && rB.bOnForwStrand
                  ? forwardToForward
                  : ( rA.bOnForwStrand && !rB.bOnForwStrand
                          ? forwardToReverse
                          : ( !rA.bOnForwStrand && rB.bOnForwStrand ? reverseToForward : reverseToReverse ) ) )
    {} // constructor

    SvJump( const nucSeqIndex uiFrom,
            const nucSeqIndex uiTo,
            const nucSeqIndex uiQueryDistance,
            const SeedOrientation xSeedOrientation )
        : uiFrom( uiFrom ), uiTo( uiTo ), uiQueryDistance( uiQueryDistance ), xSeedOrientation( xSeedOrientation )
    {} // constructor

    bool does_switch_strand( ) const
    {
        return xSeedOrientation == forwardToReverse || xSeedOrientation == reverseToForward;
    } // method

    nucSeqIndex from_start( ) const
    {
        return uiFrom - 3;
    } // method

    nucSeqIndex from_size( ) const
    {
        return 1;
    } // method

    nucSeqIndex from_end( ) const
    {
        return from_start( ) + from_size( ) - 1;
    } // method

    nucSeqIndex to_start( ) const
    {
        return 0;
    } // method

    nucSeqIndex to_size( ) const
    {
        return 1;
    } // method

    nucSeqIndex to_end( ) const
    {
        return to_start( ) + to_size( ) - 1;
    } // method

    float score( ) const
    {
        return 0;
    } // method
}; // class

}; // namespace libMA

#ifdef WITH_PYTHON
void exportSVJump( py::module& rxPyModuleId );
#endif