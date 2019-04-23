#include "container/svJump.h"

using namespace libMA;

#ifdef WITH_PYTHON

void exportSVJump( py::module& rxPyModuleId )
{
    // export the SvJump class
    py::class_<SvJump>( rxPyModuleId, "SvJump" )
        .def( py::init<const Seed&, const Seed&>( ) )
        .def( "does_switch_strand", &SvJump::does_switch_strand )
        .def( "from_start", &SvJump::from_start )
        .def( "from_size", &SvJump::from_size )
        .def( "from_end", &SvJump::from_end )
        .def( "to_start", &SvJump::to_start )
        .def( "to_size", &SvJump::to_size )
        .def( "to_end", &SvJump::to_end )
        .def( "score", &SvJump::score )
        .def_readonly( "query_distance", &SvJump::uiQueryDistance )
        .def_readonly( "seed_orientation", &SvJump::xSeedOrientation )
        .def_readonly( "from", &SvJump::uiFrom )
        .def_readonly( "to", &SvJump::uiTo );

    py::enum_<SvJump::SeedOrientation>(rxPyModuleId, "SeedOrientation")
        .value("forward_to_forward", SvJump::SeedOrientation::forwardToForward)
        .value("reverse_to_reverse", SvJump::SeedOrientation::reverseToReverse)
        .value("forward_to_reverse", SvJump::SeedOrientation::forwardToReverse)
        .value("reverse_to_forward", SvJump::SeedOrientation::reverseToForward);
} // function
#endif