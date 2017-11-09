#include "sweepAllReturnBest.h"

void exportSweepAll()
{
    //export the LineSweepContainer class
	boost::python::class_<
            SweepAllReturnBest, 
            boost::python::bases<CppModule>,
            std::shared_ptr<SweepAllReturnBest>
        >(
        "SweepAllReturnBest",
        "Uses linesweeping to remove contradicting "
        "matches within several strips of consideration.\n",
        boost::python::init<std::shared_ptr<CppModule>>()
    );
	boost::python::implicitly_convertible< 
		std::shared_ptr<SweepAllReturnBest>,
		std::shared_ptr<CppModule> 
	>();
}//function