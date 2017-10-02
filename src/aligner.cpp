#include "aligner.h"

/**
 * @brief The boost-python main method.
 *
 * A main function that exposes all Containers and Modules to python
 * by calling the respective export*() functions.
 * The main method is generated using boost-python.
 */
BOOST_PYTHON_MODULE(LAuS)
{
        exportContainer();
        exportModule();
        exportFM_index();
        exportSequence();
        exportSegmentation();
        exportPack();
        exportIntervalTree();
        exportExceptions();
        exportGetAnchors();
        exportAlignment();
        exportLinesweep();
        exportNeedlemanWunsch();
        exportSeed();
        exportBucketing();
}