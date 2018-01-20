#include "util/export.h"

/**
 * @brief The boost-python main method.
 *
 * A main function that exposes all Containers and Modules to python
 * by calling the respective export*() functions.
 * The main method is generated using boost-python.
 */
BOOST_PYTHON_MODULE(libMABS)
{
        DEBUG_3(
                std::cout.setf(std::ios::unitbuf);
        )
        exportContainer();
        exportModule();
        exportFM_index();
        exportSequence();
        exportBinarySeeding();
        exportPack();
        exportIntervalTree();
        exportExceptions();
        exportSeed();
        exportAlignment();
        exportLinesweep();
        exportNeedlemanWunsch();
        exportStripOfConsideration();
        exportChaining();
        exportSMW();
        exportExtractAllSeeds();
        exportExecOnVector();
        exportReSeed();
        exportFileReader();
        exportMappingQuality();
        exportPipe();
        exportPairedReads();
}//function

/**
 * @brief main function of executable
 * @details
 * main function if program is used as executable
 */
int main(int argc, char *argv[])
{
    std::cout << "HELLO WORLD" << std::endl;
    return 0;
}//main function