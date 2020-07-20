#include "msv/container/svJump.h"
#include "ms/util/pybind11.h"

#ifdef WITH_PYTHON
using namespace libMSV;

void exportSVJump( libMS::SubmoduleOrganizer& xOrganizer )
{
    // export the SvJump class
    py::class_<SvJump, std::shared_ptr<SvJump>>( xOrganizer.util(), "SvJump" )
        .def( py::init<const Seed&, const Seed&, int64_t>( ) )
        .def( py::init<const Seed&, nucSeqIndex, const bool, int64_t, nucSeqIndex>( ) )
        .def( "does_switch_strand", &SvJump::does_switch_strand )
        .def( "switch_strand_known", &SvJump::switch_strand_known )
        .def( "from_start", &SvJump::from_start )
        .def( "from_start_same_strand", &SvJump::from_start_same_strand )
        .def( "from_size", &SvJump::from_size )
        .def( "from_end", &SvJump::from_end )
        .def( "to_start", &SvJump::to_start )
        .def( "to_size", &SvJump::to_size )
        .def( "to_end", &SvJump::to_end )
        .def( "num_supp_nt", &SvJump::numSupportingNt )
        .def( "fuzziness", &SvJump::fuzziness )
        .def( "from_fuzziness_is_rightwards", &SvJump::from_fuzziness_is_rightwards )
        .def( "to_fuzziness_is_downwards", &SvJump::to_fuzziness_is_downwards )
        .def( "query_distance", &SvJump::query_distance )
        .def( "ref_distance", &SvJump::ref_distance )
        .def( "insert_ratio", &SvJump::insert_ratio )
        .def( "from_known", &SvJump::from_known )
        .def( "to_known", &SvJump::to_known )
        .def( "ref_ambiguity", &SvJump::referenceAmbiguity )
        .def_readonly( "from_forward", &SvJump::bFromForward )
        .def_readonly( "to_forward", &SvJump::bToForward )
        .def_readonly( "was_mirrored", &SvJump::bWasMirrored )
        .def_readonly( "from_pos", &SvJump::uiFrom )
        .def_readonly( "to_pos", &SvJump::uiTo )
        .def_readonly( "id", &SvJump::iId )
        .def_readonly( "read_id", &SvJump::iReadId );

    // export the SvCall class
    py::bind_vector<std::vector<int64_t>>( xOrganizer._util(), "int64_tVector", "docstr" );
    py::class_<SvCall, std::shared_ptr<SvCall>>( xOrganizer.util(), "SvCall" )
        .def( py::init<nucSeqIndex, nucSeqIndex, nucSeqIndex, nucSeqIndex, bool, bool, uint32_t>( ) )
        .def( py::init<nucSeqIndex, nucSeqIndex, nucSeqIndex, nucSeqIndex, bool, bool, uint32_t, uint32_t>( ) )
        .def( py::init<std::shared_ptr<SvJump>>( ) )
        .def( "join", &SvCall::join )
        .def( "clear_jumps", &SvCall::clear_jumps )
        .def( "add_jump", &SvCall::add_jump )
        .def( "get_jump", &SvCall::get_jump )
        .def( "get_score", &SvCall::getScore )
        .def_readwrite( "num_supp_reads", &SvCall::uiNumSuppReads )
        .def_readwrite( "reference_ambiguity", &SvCall::uiReferenceAmbiguity )
        .def_readwrite( "x", &SvCall::xXAxis )
        .def_readwrite( "y", &SvCall::xYAxis )
        .def_readwrite( "from_forward", &SvCall::bFromForward )
        .def_readwrite( "to_forward", &SvCall::bToForward )
        .def_readwrite( "supporing_jump_ids", &SvCall::vSupportingJumpIds )
        .def_readwrite( "inserted_sequence", &SvCall::pInsertedSequence )
        .def_readonly( "id", &SvCall::iId );

} // function
#endif