#include "aligner.h"

/**
 * @brief The boost-python main method.
 *
 * A main function that exposes all Containers and Modules to python
 * by calling the respective export*() fonctions.
 * The main method is generated using boost-python.
 */
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
        exportNeedlemanWunsch();
}