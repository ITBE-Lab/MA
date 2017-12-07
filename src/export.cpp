#include "export.h"

/**
 * @brief The boost-python main method.
 *
 * A main function that exposes all Containers and Modules to python
 * by calling the respective export*() functions.
 * The main method is generated using boost-python.
 */
BOOST_PYTHON_MODULE(libLAuS)
{
        DEBUG_3(
                std::cout.setf(std::ios::unitbuf);
        )
        exportContainer();
        exportModule();
        exportFM_index();
        exportSequence();
        exportLongestNonEnclosedSegments();
        exportLongestLRSegments();
        exportPack();
        exportIntervalTree();
        exportExceptions();
        exportGetAnchors();
        exportAlignment();
        exportLinesweep();
        exportNeedlemanWunsch();
        exportSeed();
        exportStripOfConsideration();
        exportChaining();
        exportSMW();
        exportExtractAllSeeds();
        exportExecOnVector();
        exportReSeed();
}