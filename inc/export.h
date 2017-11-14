/** 
 * @file export.h
 * @brief Provides the boost-python main method.
 * @author Markus Schmidt
 * @details 
 * Calls the export functions of the various @ref CppModule "modules" and Containers.
 */

/**
 * @mainpage @ref CppModule "Modular" @ref LineSweep "Linesweeping" Aligner using Segmentation
 * @tableofcontents
 * @section intro_sec Introduction
 *
 * LAuS is a Modular alignment tool build using C++11 and Boost Python.
 * The alignment process has been seperated into several @ref CppModule "modules".
 * The execution order of the @ref CppModule "modules" is set up using Python.
 * @ref CppModule "Modules" can be implemented in Python or C++. <br>
 * The Pledge class allows setting up a @ref comp_graph_sec "computational graph", 
 * that avoids unnecessary jumps between 
 * Python and C++.
 *
 *
 * The general aligner structure is as follows:
 * - Create seeds
 * - Remove inconsistent seeds
 * - Create local alignments in the gaps between the seeds
 * 
 * 
 * A C++ Module is available for each of these tasks, respectively:
 * - LongestNonEnclosedSegments
 * - LineSweep
 * - NeedlemanWunsch
 * 
 * 
 * @ref Container "Containers" are used for the inputs and outputs of the @ref CppModule "modules".
 * The Expected in- and out- puts for each of the three main steps is as follows:
 * <table>
 * <caption>inputs and outputs for each main step in the alignment</caption>
 * <tr><th>Step <th>input <th>output
 * <tr><td>Create seeds <td> FM_Index, NucleotideSequence <td> SegmentTree
 * <tr><td>Remove inconsistent seeds <td> SegmentTree <td> SegmentTree
 * <tr><td>Create local alignments in the gaps between the seeds <td> SegmentTree,
 * NucleotideSequence, BWACompatiblePackedNucleotideSequencesCollection <td> Alignment
 * </table>
 *
 * @note The python classes can be easily identified by the prefix "LAuS." 
 * while C++ classes have no specific prefix.
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
 * Here is some python code that sets up the three main @ref CppModule "modules" reqired for alignment:
 * @code{.py}
 * # A module that creates seeds.
 * seg = LongestNonEnclosedSegments()
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
 * Note that while LongestNonEnclosedSegments and NeedlemanWunsch are C++ modules,
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
 * segments = seg.execute((fm_index, query)) #(*)
 * # Call the line sweep module.
 * seeds = sweep.execute((segments,)) #(*)
 * # Call the local alignment module.
 * alignment = nmw.execute((seeds, query, ref)) #(*)
 * 
 * # Print the alignment
 * printer.execute((alignment, )) #(*)
 * @endcode
 * 
 * Here we perform the alignment process and print the results. <br>
 * (*) All @ref Module "modules" use a tuple or list of @ref Container "containers" as input.
 * Therefore the @ref Module::execute "execute" function calls use the syntax shown above.
 *
 * @note for a quick start on how to setup a computational graph,
 * this setup of @ref CppModule "modules" using a
 * computation graph can be seen in the @ref comp_graph_sec "Pledge" class.
 * 
 * @section todos TODOs
 * 
 * - update this page... (the example has to be changed)
 * - rename classes and files aproproately
 * - move backwards extension into fm_index class
 * - STOP ADDING ITEMS TO THIS LIST
 * - reintroduce rotations to searchTree
 * - do runtime breakdown
 * 
 * @section about_us_sec About Us
 * 
 * Here is the <a href="http://itbe.hanyang.ac.kr/">webpage of the IT-BE lab</a>.
 */

/**
 * @defgroup export
 * @brief functions that are used to export Container and @ref CppModule "module" classes to Python
 * @details
 * When the library is imported in python we need to tell python which classes functions etc.
 * we provide.
 */

#ifndef EXPORT_H
#define EXPORT_H

#include "getAnchors.h"
#include "needlemanWunsch.h"
#include "linesweep.h"
#include "longestNonEnclosedSegments.h"
#include "bucketing.h"
#include "longestLRSegments.h"
#include "chaining.h"
#include "smith_waterman.h"
#include "extractAllSeeds.h"
#include "execOnVector.h"


#endif