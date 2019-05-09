#include "container/svJump.h"

using namespace libMA;

#ifdef WITH_PYTHON

void exportSVJump( py::module& rxPyModuleId )
{
    // export the SvJump class
    py::class_<SvJump>( rxPyModuleId, "SvJump" )
        .def( py::init<const Seed&, const Seed&, const bool>( ) )
        .def( "does_switch_strand", &SvJump::does_switch_strand )
        .def( "from_start", &SvJump::from_start )
        .def( "from_start_same_strand", &SvJump::from_start_same_strand )
        .def( "from_size", &SvJump::from_size )
        .def( "from_end", &SvJump::from_end )
        .def( "to_start", &SvJump::to_start )
        .def( "to_size", &SvJump::to_size )
        .def( "to_end", &SvJump::to_end )
        .def( "score", &SvJump::score )
        .def( "from_fuzziness_is_rightwards", &SvJump::from_fuzziness_is_rightwards )
        .def( "to_fuzziness_is_downwards", &SvJump::to_fuzziness_is_downwards )
        .def( "query_distance", &SvJump::query_distance )
        .def( "ref_distance", &SvJump::ref_distance )
        .def( "insert_ratio", &SvJump::insert_ratio )
        .def_readonly( "from_forward", &SvJump::bFromForward )
        .def_readonly( "to_froward", &SvJump::bToForward )
        .def_readonly( "from_seed_start", &SvJump::bFromSeedStart )
        .def_readonly( "from_pos", &SvJump::uiFrom )
        .def_readonly( "to_pos", &SvJump::uiTo )
        .def_readonly( "id", &SvJump::iId )
        .def_readonly( "read_id", &SvJump::iReadId );

    // export the SvCall class
    py::bind_vector<std::vector<int64_t>>(rxPyModuleId, "int64_tVector", "docstr");
    py::class_<SvCall>( rxPyModuleId, "SvCall" )
        .def( py::init<nucSeqIndex, nucSeqIndex, nucSeqIndex, nucSeqIndex, bool, double>( ) )
        .def( py::init<const SvJump&>( ) )
        .def( py::init<SvJump&>( ) )
        .def( "join", &SvCall::join )
        .def( "clear_jumps", &SvCall::clear_jumps )
        .def( "add_jump", &SvCall::add_jump )
        .def( "get_jump", &SvCall::get_jump )
        .def_readwrite( "score", &SvCall::dScore )
        .def_readwrite( "from_start", &SvCall::uiFromStart )
        .def_readwrite( "from_size", &SvCall::uiFromSize )
        .def_readwrite( "to_start", &SvCall::uiToStart )
        .def_readwrite( "to_size", &SvCall::uiToSize )
        .def_readwrite( "switch_strand", &SvCall::bSwitchStrand )
        .def_readwrite( "supporing_jump_ids", &SvCall::vSupportingJumpIds )
        .def_readonly( "id", &SvCall::iId );

} // function
#endif