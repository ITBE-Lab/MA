#include "container/alignment.h"
using namespace libLAuS;

void exportAlignment()
{
    boost::python::class_<
            Alignment, 
            boost::python::bases<Container>, 
            std::shared_ptr<Alignment>
        >(
                "Alignment",
                "contains the final output of the aligner\n"
            )
        .def(
                "at", 
                &Alignment::at,
                "arg1: self\n"
                "arg2: index at which to look\n"
                "returns: match type at the given position\n"
            )
        .def(boost::python::init<nucSeqIndex>())
        .def(
                "__getitem__", 
                &Alignment::at,
                "arg1: self\n"
                "arg2: index at which to look\n"
                "returns: match type at the given position\n"
            )
        .def(
                "append", 
                &Alignment::append_boost1,
                "arg1: self\n"
                "arg2: the matchtype to append\n"
                "arg3: how many times shall the matchtype be appended\n"
                "returns: nil\n"
            )
        .def(
                "append", 
                &Alignment::append_boost2,
                "arg1: self\n"
                "arg2: the matchtype to append\n"
                "returns: nil\n"
            )
        .def(
                "begin_on_ref", 
                &Alignment::beginOnRef,
                "arg1: self\n"
                "returns: starting position of the alignment on the reference\n"
            )
        .def(
                "end_on_ref", 
                &Alignment::endOnRef,
                "arg1: self\n"
                "returns: ending position of the alignment on the reference\n"
            )
        .def(
                "__len__", 
                &Alignment::length,
                "arg1: self\n"
                "returns: length of the alignmen\n"
            )
        .def(
                "length", 
                &Alignment::length,
                "arg1: self\n"
                "returns: length of the alignmen\n"
            )
        .def(
                "seed_coverage", 
                &Alignment::seedCoverage,
                "arg1: self\n"
            )
        .def(
                "get_score", 
                &Alignment::score,
                "arg1: self\n"
            )
        .def(
                "num_by_seeds", 
                &Alignment::numBySeeds,
                "arg1: self\n"
            )
    ;

    
    //export the matchType enum
    boost::python::enum_<MatchType>("MatchType")
        .value("match", MatchType::match)
        .value("seed", MatchType::seed)
        .value("missmatch", MatchType::missmatch)
        .value("insertion", MatchType::insertion)
        .value("deletion", MatchType::deletion);

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible< 
        std::shared_ptr<Alignment>,
        std::shared_ptr<Container> 
    >();

}