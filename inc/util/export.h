/** 
 * @file export.h
 * @brief Provides the boost-python main method.
 * @author Markus Schmidt
 * @details 
 * Calls the export functions of the various @ref Module "modules" and Containers.
 */

/**
 * @mainpage @ref libLAuS::Module "Modular" @ref libLAuS::LinearLineSweep "Linesweeping" Aligner using Segmentation
 * @tableofcontents
 * @section intro_sec Introduction
 *
 * LAuS is a @ref libLAuS::Module "Modular" alignment tool build using C++11 and Boost Python.
 * The alignment process has been seperated into several @ref libLAuS::Module "modules".
 * The execution order of the @ref libLAuS::Module "modules" is set up using Python.
 * @ref libLAuS::Module "Modules" can be implemented in Python or C++. <br>
 * The @ref libLAuS::Pledge "Pledge" class allows setting up a @ref comp_graph_sec "computational graph", 
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
 * <tr><td>Seeding <td> libLAuS::FMIndex, libLAuS::NucSeq <td> libLAuS::SegmentVector
 * <tr><td>Seed Processing <td> libLAuS::SegmentVector <td> libLAuS::SegmentVector
 * <tr><td>Optimal Matching <td> libLAuS::SegmentVector,
 * libLAuS::NucSeq, libLAuS::Pack <td> libLAuS::Alignment
 * </table>
 *
 * @note The python classes can be easily identified by the prefix "LAuS." 
 * while C++ classes have "libLAuS::" as prefix.
 *
 * @section install_sec Installation
 * 
 * The easiest way to install the LAuS library is using pip:
 * @code{.sh}
 * pip install LAuS
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
 *
 * @note for a quick start on how to setup a computational graph,
 * this setup of @ref libLAuS::Module "modules" using a
 * computation graph can be seen in the @ref comp_graph_sec "Pledge" class.
 * 
 * @section todos TODOs
 * 
 * - make banded efficient NMW (have to wait for SMW)
 * - STOP ADDING ITEMS TO THIS LIST
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

#include "module/getAnchors.h"
#include "module/needlemanWunsch.h"
#include "module/linesweep.h"
#include "module/binarySeeding.h"
#include "module/stripOfConsideration.h"
#include "module/chaining.h"
#include "module/smith_waterman.h"
#include "module/extractAllSeeds.h"
#include "module/execOnVector.h"
#include "module/reSeed.h"


#endif