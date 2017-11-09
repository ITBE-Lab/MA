#include "sweepAllReturnBest.h"

void exportSweepAll()
{
    //export the LineSweepContainer class
	boost::python::class_<SweepAllReturnBest, boost::python::bases<CppModule>>(
        "SweepAllReturnBest",
        "Uses linesweeping to remove contradicting "
        "matches within several strips of consideration.\n",
        boost::python::init<CppModule>()
    );
}//function