/** 
 * @file export.h
 * @brief Provides the boost-python main method.
 * @author Markus Schmidt
 * @details 
 * Calls the export functions of the various @ref Module "modules" and Containers.
 */

/**
 * @mainpage @ref libMA::Module "Modular" Aligner using @ref libMA::BinarySeeding "BinarySeeding" and @ref libMA::StripOfConsideration "Strips of Consideration"
 * @tableofcontents
 * @section intro_sec Introduction
 *
 * MA is a @ref libMA::Module "Modular" alignment tool build using C++11 and Boost Python.
 * The alignment process has been seperated into several @ref libMA::Module "modules".
 * The execution order of the @ref libMA::Module "modules" is set up using Python.
 * @ref libMA::Module "Modules" can be implemented in Python or C++. <br>
 * The @ref libMA::Pledge "Pledge" class allows setting up a @ref comp_graph_page "computational graph", 
 * that avoids unnecessary jumps between 
 * Python and C++. <br>
 *
 *
 * The general aligner structure is as follows:
 * - Seeding
 * - Seed Processing
 *   -# Filtering
 *   -# Coupeling
 * - Optimal Matching
 * 
 * 
 * A C++ Module is available for each of these tasks, respectively.
 * 
 * 
 * @ref Container "Containers" are used for the inputs and outputs of the @ref Module "modules".
 * The Expected in- and out- puts for each of the three main steps is as follows:
 * <table>
 * <caption>inputs and outputs for each main step in the alignment</caption>
 * <tr><th>Step <th>input <th>output
 * <tr><td>Seeding <td> libMA::FMIndex, libMA::NucSeq <td> libMA::SegmentVector
 * <tr><td>Seed Processing <td> libMA::SegmentVector <td> libMA::SegmentVector
 * <tr><td>Optimal Matching <td> libMA::SegmentVector,
 * libMA::NucSeq, libMA::Pack <td> libMA::Alignment
 * </table>
 *
 * @note The python classes can be easily identified by the prefix "MA." 
 * while C++ classes have "libMA::" as prefix.
 *
 * @section install_sec Installation
 * 
 * The easiest way to install the MA library is using pip:
 * @code{.sh}
 * pip install MA
 * @endcode
 * There is a git repo with the full source code here: .... <br>
 * 
 * @section quick_start_sec Quick start
 *
 * Here is some python code that sets up the three main @ref Module "modules" reqired for alignment:
 * @code{.py}
 * # A module that creates seeds.
 * seg = BinarySeeding()
 * # A module that removes inconsistent seeds.
 * sweep = SweepAllReturnBest()
 * # A module that creates local alignments in the gaps between the seeds.
 * nmw = NeedlemanWunsch()
 * # A module that prints the alignment to the console.
 * printer = AlignmentPrinter()
 * @endcode
 *
 * Here we set up the @ref Module "modules" we need for the alignment process.
 * The @ref Module "modules" themselves do not store data. They can be used multiple times.
 * Note that while BinarySeeding and NeedlemanWunsch are C++ modules,
 * SweepAllReturnBest and AlignmentPrinter are implemented in Python.
 *
 * @code{.py}
 * # Setup a container for the suffix array
 * fm_index = FMIndex()
 * # Load the array from a file.
 * fm_index.load("filename")
 * 
 * # Setup a container for the packed reference
 * ref = Pack()
 * # Load the pack from a file.
 * ref.load("filename")
 * 
 * # Create a query string.
 * query_string = "ACCTAA"
 * # Setup a container for the query sequence.
 * query = NucSeq(query_string)
 * @endcode
 *
 * Here we setup the @ref Container "containers" holding the data we need.
 *
 * @code{.py}
 * # Call the segmentation module.
 * segments = seg.execute(fm_index, query)
 * # Call the line sweep module.
 * seeds = sweep.execute(segments)
 * # Call the local alignment module.
 * alignment = nmw.execute(seeds, query, ref)
 * 
 * # Print the alignment
 * printer.execute(alignment)
 * @endcode
 * 
 * Here we perform the alignment process and print the results. <br>
 * Note that this setup does not perform well, since we skipped the seed filtering step.
 *
 * @note for a quick start on how to setup a computational graph,
 * this setup of @ref libMA::Module "modules" using a
 * computation graph can be seen @ref comp_graph_page "here".
 * 
 * @section todos TODOs
 * 
 * @todo make efficient NMW (have to wait for SMW)
 * @todo paired end alignments necessary
 * @todo different output formats (queue implementation...)
 * @todo code cleanup
 * @todo quality of SMW
 * @todo figure out where sequencer reads belong in the missmatch/indel-pic
 * @todo be able to deal with short reads
 * 
 * @section about_us_sec About Us
 * 
 * Here is the <a href="http://itbe.hanyang.ac.kr/">webpage of the IT-BE lab</a>.
 */

/**
 * @defgroup export
 * @brief functions that are used to export Container and @ref Module "module" classes to Python
 * @details
 * When the library is imported in python we need to tell python which classes functions etc.
 * we provide.
 */

#ifndef EXPORT_H
#define EXPORT_H

#include "module/needlemanWunsch.h"
#include "module/linesweep.h"
#include "module/binarySeeding.h"
#include "module/stripOfConsideration.h"
#include "module/chaining.h"
#include "module/smith_waterman.h"
#include "module/extractAllSeeds.h"
#include "module/execOnVector.h"
#include "module/reSeed.h"
#include "module/fileReader.h"
#include "module/fileWriter.h"
#include "module/mappingQuality.h"
#include "module/pairedReads.h"
#include "module/splitter.h"
#include "container/container.h"

std::vector<std::shared_ptr<libMA::Pledge>> EXPORTED setUpCompGraph(
    std::shared_ptr<libMA::Pledge> pPack,
    std::shared_ptr<libMA::Pledge> pFMDIndex,
    std::shared_ptr<libMA::Pledge> pQueries,
    std::shared_ptr<libMA::Module> pOut,
    unsigned int uiThreads,
    unsigned int uiMaxAmbiguity,
    unsigned int uiNumSOC,
    bool bPariedNormal,
    bool bPariedUniform,
    unsigned int uiPairedMean,
    double fPairedStd,
    double dPairedU,
    unsigned int reportNBest
);

#endif