/** 
 * @file export.h
 * @brief Provides the boost-python main method.
 * @author Markus Schmidt
 * @details 
 * Calls the export functions of the various @ref Module "modules" and Containers.
 * @copyright
Copyright 2018 Markus Schmidt, Arne Kutzner

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @mainpage @ref libMA::Module "Modular" Aligner using @ref libMA::BinarySeeding "BinarySeeding" and @ref libMA::StripOfConsideration "Strips of Consideration"
 * @tableofcontents
 * @section intro_sec Introduction
 *
 * @note This documentation is outdated!
 *
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
 * - Seed set assembling
 *   -# Filtering
 *   -# Strip of Consideration
 *   -# Harmonization
 * - Dynamic programming
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
 * <tr><td>Seed set assembling <td> libMA::SegmentVector <td> libMA::Seeds
 * <tr><td>Dynamic programming <td> libMA::Seeds, libMA::NucSeq, libMA::Pack <td> libMA::Alignment
 * </table>
 *
 * @note The python classes can be easily identified by the prefix "MA." 
 * while C++ classes have "libMA::" as prefix.
 *
 * @section install_sec Installation
 * 
 * The easiest way to install the MA library is using git hub:
 * @code{.sh}
 * git clone https://github.com/ItBeLab/ma
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
 * @note this setup of @ref libMA::Module "modules" using a
 * computation graph can be seen @ref comp_graph_page "here".
 * 
 * @section todos TODOs
 * 
 * @todo different output formats
 * @todo code cleanup
 * @todo remove boost::program_options
 * @todo remove boost::filesystem
 * @todo remove boost::iostreams
 * 
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
#include "module/extractAllSeeds.h"
#include "module/execOnVector.h"
#include "module/fileReader.h"
#include "module/fileWriter.h"
#include "module/mappingQuality.h"
#include "module/pairedReads.h"
#include "module/splitter.h"

std::vector<std::shared_ptr<libMA::Pledge>> EXPORTED setUpCompGraph(
    std::shared_ptr<libMA::Pledge> pPack,
    std::shared_ptr<libMA::Pledge> pFMDIndex,
    std::vector<std::shared_ptr<libMA::Pledge>> aQueries,
    std::shared_ptr<libMA::Module> pOut,
    unsigned int uiReportN,
    unsigned int uiThreads,
    bool bPariedNormal,
    bool bPariedUniform,
    unsigned int uiPairedMean,
    double fPairedStd,
    double dPairedU,
    bool bSeedSetPairs,
    float fGiveUp,
    unsigned int iMatch_,
    unsigned int iMisMatch_,
    unsigned int iGap_,
    unsigned int iExtend_
);

#endif