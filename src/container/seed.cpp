/**
 * @file seed.cpp
 * @author Markus Schmidt
 */
#include "container/seed.h"
#include "util/default_parameters.h"
using namespace libMA;

#ifdef _MSC_VER
// getting ambiguous abs otherwise
#include <cstdlib>
#endif

using namespace libMA::defaults;
extern int libMA::defaults::iGap;
extern int libMA::defaults::iExtend;
extern int libMA::defaults::iMatch;
extern int libMA::defaults::iMissMatch;

/*returns the sum off all scores within the list*/
nucSeqIndex Seeds::getScore( ) const
{
#if 0
    if(bConsistent)
    {
        nucSeqIndex iRet = 0;
        nucSeqIndex uiLastQPos = 0;
        nucSeqIndex uiLastRPos = 0;
        for(const Seed& rS : *this)
        {
            iRet += rS.getValue() * iMatch;
            if(uiLastQPos != 0 && uiLastRPos != 0)
            {
                nucSeqIndex uiQDist = rS.start() - uiLastQPos;
                nucSeqIndex uiRDist = rS.start_ref() - uiLastRPos;
                iRet -= iGap + iExtend * std::min(uiQDist, uiRDist);
                nucSeqIndex uiDist = uiQDist - uiRDist;
                if (uiRDist > uiQDist)
                    uiDist = uiRDist - uiQDist;
                iRet -= iMissMatch * uiDist;
                uiLastQPos = rS.start();
                uiLastRPos = rS.start_ref();
            }//if
        }//for
        return iRet;
    }//if
    else
    {
#endif
    nucSeqIndex iRet = 0;
    for( const Seed &rS : *this )
        iRet += rS.getValue( );
    return iRet;
    //    }//else
} // function

#ifdef WITH_PYTHON
void exportSeed( )
{
    exportInterval<nucSeqIndex>( );
    // export the Seed class
    boost::python::class_<Seed, boost::python::bases<Interval<nucSeqIndex>>>( "Seed" )
        .def_readwrite( "start_ref", &Seed::uiPosOnReference )
        .def( "__eq__", &Seed::operator==);


    // export the Seed class
    boost::python::class_<AlignmentStatistics>( "AlignmentStatistics", boost::python::init<>( ) )
        .def_readwrite( "index_of_strip", &AlignmentStatistics::index_of_strip )
        .def_readwrite( "num_seeds_in_strip", &AlignmentStatistics::num_seeds_in_strip )
        .def_readwrite( "anchor_size", &AlignmentStatistics::anchor_size )
        .def_readwrite( "anchor_ambiguity", &AlignmentStatistics::anchor_ambiguity )
        .def_readwrite( "name", &AlignmentStatistics::sName )
        .def_readwrite( "initial_q_beg", &AlignmentStatistics::uiInitialQueryBegin )
        .def_readwrite( "initial_r_beg", &AlignmentStatistics::uiInitialRefBegin )
        .def_readwrite( "initial_q_end", &AlignmentStatistics::uiInitialQueryEnd )
        .def_readwrite( "initial_r_end", &AlignmentStatistics::uiInitialRefEnd );

    // export the Seeds class
    boost::python::class_<Seeds, boost::python::bases<Container>,
                          boost::python::bases<std::list<Seed>>, std::shared_ptr<Seeds>>( "Seeds" )
        .def( boost::python::init<std::shared_ptr<Seeds>>( ) )
        .def(
            boost::python::vector_indexing_suite<Seeds,
                                                 /*
                                                  *    true = noproxy this means that the content of
                                                  * the vector is already exposed by boost python.
                                                  *    if this is kept as false, Container would be
                                                  * exposed a second time. the two Containers would
                                                  * be different and not inter castable.
                                                  */
                                                 true>( ) );

    // make vectors of container-pointers a thing
    IterableConverter( ).from_python<Seeds>( );

    // tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<std::shared_ptr<Seeds>, std::shared_ptr<Container>>( );
} // function
#endif