#include "aligner.h"

//this creates our main function
BOOST_PYTHON_MODULE(LAuS)
{
        exportModule();
        exportContainer();
        exportFM_index();
        exportSequence();
        exportSegmentation();
        exportPack();
        exportIntervalTree();
        exportGraphicalMethod();
        exportExceptions();
        exportGetAnchors();
        exportAlignment();
        exportLinesweep();
        //exportNeedlemanWunsch();
}