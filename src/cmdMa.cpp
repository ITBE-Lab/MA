/** 
 * @file cmdMa.cpp
 * @brief Parses the program options
 * @details
 * Sets up the computational graph and executes it.
 * @author Markus Schmidt
 * @copyright
Copyright 2018 Markus Schmidt, Arne Kutzner

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "container/nucSeq.h"
#include "container/pack.h"
#include "container/fMIndex.h"
#include "module/fileReader.h"
#include "module/fileWriter.h"
#include "util/export.h"
#include "util/cxxopts.hpp"
#include "util/default_parameters.h"
#include "util/debug.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <thread>
#include <iostream>
#include <string.h>
/// @endcond

using namespace libMA;
using namespace cxxopts;

const std::string sHelp = 
"====================================== The Modular Aligner ======================================"
"\nGeneral options:"
"\n    -h, --help                     Display the complete help screen"
"\n        --genIndex                 Do FMD-index Generation. The -i and -x options specify the"
"\n                                   FASTA file used for index generation and the index prefix," 
"\n                                   respectively. If this option is not set, the aligner performs"
"\n                                   alignments. "
"\n"
"\nNecessary arguments for alignments:"
"\n    -x, --idx <prefix>             FMD-index used for alignments"
"\n    -i, --in <fname>               FASTA or FASTAQ input files."
"\n"
"\nAvailable alignment presettings:"
"\n     -m, --mode [str]              Set operation modus for MA."
"\n                                   fast:"
"\n                                       fast does...."
"\n                                       fast does...."
"\n                                       fast does...."
"\n                                   acc:"
"\n                                       acc does...."
"\n                                   pacBio:"
"\n                                       pacBio does...."
"\n                                   Default is 'fast'."
"\n"
"\nAlignments options:"
"\n    -o, --out <fname>              Filename used for SAM file output. Default output stream is"
"\n                                   standard output."
"\n    -t, --threads <num>            Use <num> threads. On startup MA checks the hardware and "
"\n                                   chooses this value accordingly."
"\n    -d, --noDP                     Switch that disables the final Dynamic Programming."
"\n    -n, --reportN <num>            Report up to <num> alignments; 0 means unlimited."
"\n                                   Default is 0."
"\n    -s, --seedSet [SMEMs/maxSpan]  Selects between the two seeding strategies super maximal"
"\n                                   extended matches 'SMEMs' and maximally spanning seeds"
"\n                                   'maxSpan'."
"\n                                   Default is 'maxSpan'."
"\n    -l, --minLen <num>             Seeds must have a minimum length of <num> nucleotides."
"\n                                   Default is 16."
"\n        --Match <num>              Sets the match score to <num>; <num> > 0."
"\n                                   Default is 3. "
"\n        --MisMatch <num>           Sets the mismatch penalty to <num>; <num> > 0."
"\n                                   Default is 4."
"\n        --Gap <num>                Sets the costs for opening a gap to <num>; <num> >= 0."
"\n                                   Default is 6."
"\n        --Extend <num>             Sets the costs for extending a gap to <num>; <num> > 0."
"\n                                   Default is 1"
"\n"
"\nPaired Reads options:"
"\n    -p, --paUni                    Enable paired alignments and model the distance as"
"\n                                   uniform distribution."
"\n                                   If set --in shall be used as follows:"
"\n                                   --in '<fname1>, <fname2>'."
"\n    -P, --paNorm                   Enable paired alignment and Model the distance as" 
"\n                                   normal distribution."
"\n                                   If set --in shall be used as follows:"
"\n                                   '--in <fname1>, <fname2>'."
"\n        --paIsolate <num>          Penalty for an unpaired read pair."
"\n                                   Default is 17."
"\n        --paMean <num>             Mean gap distance between read pairs."
"\n                                   Default is 400."
"\n        --paStd <num>              Standard deviation of gap distance between read pairs."
"\n                                   Default is 150."
"\n"
"\nAdvanced options:"
"\n        --giveUp <val>             Threshold with 0 <= <val> <= 1 used as give-up criteria."
"\n                                   SoC's with accumulative seed length smaller than "
"\n                                   'query_len * <val>' will be ignored."
"\n                                   Reducing this parameter will decrease runtime but allows"
"\n                                   the aligner to discover more dissimilar matches."
"\n                                   Increasing this parameter will increase runtime but might"
"\n                                   cause the aligner to miss the correct reference location."
"\n                                   Default is 0.002."
"\n        --maxTries <num>           At most the best <num> SoC's will be inspected."
"\n                                   Generally the best alignment is found in the best scored SoC."
"\n                                   However, if the best alignment is from a very repetitive"
"\n                                   region, we might have to inspect several SoC's to find the"
"\n                                   optimal one."
"\n                                   Default is 50."
"\n        --minRefSize <num>         If the reference is smaller than <num> nt we disable post SoC"
"\n                                   heuristics."
"\n                                   Default is 10,000,000."
"\n        --minSecToPrimRatio <num>  Only output secondary alignments if their score is larger"
"\n                                   or equal than <num> * score of primary alignment."
"\n                                   Default is 0.75."
"\n"
"\nVersion 0.1.0 (alpha)"
"\nBy Markus Schmidt & Arne Kutzner"
"\nFor more information visit: https://github.com/ITBE-Lab/ma"
;

/**
 * main function
 */
int main(int argc, char* argv[])
{

    Options options("MA", "\t\t===== The Modular Aligner =====");
    options.add_options()
        ("h,help", "Display the complete help screen")
        ("t,threads", "Number of threads", 
           value<unsigned int>()->default_value(std::to_string(std::thread::hardware_concurrency()))
           , "arg     "
        )
    ;

    if (argc <= 1)
    {
        std::cout << "Invalid usage; you need to supply -x and -i at least!" << std::endl;
        std::cout << sHelp << std::endl;
        return 0;
    }//if

    try
    {
        defaults::configureFast();
        for(int i = 0; i < argc-1; i++)
            if(
                    strcmp(argv[i], "-m") == 0 ||
                    strcmp(argv[i], "--mode") == 0
                )
            {
                if(strcmp(argv[i+1], "acc") == 0)
                    defaults::configureAccurate();
                if(strcmp(argv[i+1], "fast") == 0)
                    defaults::configureFast();
                if(strcmp(argv[i+1], "pacBio") == 0)
                    defaults::configurePacBio();
            }// if

        options.add_options("Alignment options (requires -a)")
            ("i,in", "Input file(s) as (multi-)fasta(-q)", 
                value<std::vector<std::string>>(), "args"
            )
            ("o,out", "Output file as SAM",value<std::string>()->default_value("stdout"))
            ("m,mode", "Pre-setting [fast/acc]", 
                value<std::string>()->default_value(defaults::sParameterSet)
            )
            ("s,seedSet", "Seeding strategy [SMEMs/maxSpanning]",
                value<std::string>()->default_value(defaults::sSeedSet)
            )
            ("n,reportN", "Report up to N alignments; 0: unlimited",
                value<unsigned int>()->default_value(defaults::uiReportN)
            )
            ("l,minLen", "Minimum seed length",
                value<unsigned int>()->default_value(defaults::uiMinLen)
            )
            ("giveUp", "Minimum SoC score (relative to query length)",
                value<double>()->default_value(defaults::fGiveUp)
            )
            ("minSecToPrimRatio", "Minimum sec score ratio",
                value<double>()->default_value(defaults::fMinSecScoreRatio)
            )
            ("maxTries", "Max num SoC",
                value<unsigned int>()->default_value(defaults::uiMaxTries)
            )
            ("minRefSize", "ref size switch",
                value<unsigned long long>()->default_value(defaults::uiGenomeSizeDisable)
            )
            ("d,noDP", "Disable DP")
            ("Match", "DP match score.", value<unsigned int>()->default_value(defaults::uiMatch))
            ("MisMatch", "DP mismatch penalty.", 
                value<unsigned int>()->default_value(defaults::uiMissMatch)
            )
            ("Gap", "DP gap open penalty.", 
                value<unsigned int>()->default_value(defaults::uiOpen)
            )
            ("Extend", "DP gap extend penalty.", 
                value<unsigned int>()->default_value(defaults::uiExtend)
            )
            ("Gap2", "DP gap open penalty.", 
                value<unsigned int>()->default_value(defaults::uiOpen2)
            )
            ("Extend2", "DP gap extend penalty.", 
                value<unsigned int>()->default_value(defaults::uiExtend2)
            )
            ("x,idx", "Do FMD-index generation", value<std::string>())
            ("genIndex", "Do FMD-index generation")
            ("SoCWidth", "SoC width", value<unsigned int>()->default_value("0"))
            ("disableHeuristics", "disable all heuristics", value<bool>()->default_value(defaults::bDisableHeuristics))
        ;

        options.add_options("Paired Reads options (requires either -U or -N)")
            ("p,paUni", "Enable paired alignment; Distance as uniform distribution")
            ("P,paNorm", "Enable paired alignment; Distance as normal distribution")
            ("paIsolate", "Penalty for unpaired alignments", 
                value<double>()->default_value(defaults::uiUnpaired), "arg    "
            )
            ("paMean", "Gap distance mean", 
                value<unsigned int>()->default_value(defaults::uiMean)
            )
            ("paStd", "Gap distance standard deviation", 
                value<double>()->default_value(defaults::uiStd)
            )
        ;

        auto result = options.parse(argc, argv);
        
        if (result.count("help"))
        {
            std::cout << sHelp << std::endl;
            //@todo cmake version number: https://stackoverflow.com/questions/27395120/correct-way-to-encode-embed-version-number-in-program-code
            DEBUG(
                std::cout << "DEBUG LEVEL: " << DEBUG_LEVEL << std::endl;
            )
            return 0;
        }//if

        auto uiT =                result["threads"].      as<unsigned int>();
        auto bPariedNormal =      result.count("paNorm")  > 0;
        auto bPariedUniform =     result.count("paUni") > 0;
        auto uiPairedMean =       result["paMean"].         as<unsigned int>();
        auto fPairedStd =         result["paStd"].          as<double>();
        auto dPairedU =           result["paIsolate"].     as<double>();
        auto fMinSecScoreRatio =  result["minSecToPrimRatio"].as<double>();
        auto uiReportN =          result["reportN"].      as<unsigned int>();
        auto sParameterSet =      result["mode"]. as<std::string>();
        auto sSeedSet =           result["seedSet"].      as<std::string>();
        auto uiMinLen =           result["minLen"].       as<unsigned int>();
        auto uiSoCWidth =         result["SoCWidth"].     as<unsigned int>();
        auto bDisableHeuristics = result["disableHeuristics"].     as<bool>();
        if(bPariedNormal && bPariedUniform)
        {
            std::cerr << "--normal and --uniform are exclusive." << std::endl;
            return 1;
        }// else if
        std::string sGenome;
        if( result.count("idx") > 0 )
            sGenome =           result["idx"].       as<std::string>();
        else
        {
            std::cerr << "error: --idx is compulsory" << std::endl;
            return 1;
        }// else if
        std::vector<std::string> aIn;
        if( result.count("in") > 0 )
            aIn =          result["in"].      as<std::vector<std::string>>();
        else
        {
            std::cerr << "error: --in is compulsory" << std::endl;
            return 1;
        }// else if
        if(aIn.size() != 1 && !(bPariedNormal || bPariedUniform) && result.count("genIndex") == 0)
        {
            std::cerr << "error: --in takes one argument in unpaired mode" << std::endl;
            return 1;
        }// if
        else if(aIn.size() != 2 && (bPariedNormal || bPariedUniform) && result.count("genIndex") == 0)
        {
            std::cerr << "error: --in takes two arguments in paired mode" << std::endl;
            return 1;
        }// else if
        auto sOut =        result["out"].     as<std::string>();
        auto bFindMode =        result.count("noDP") > 0;
        auto fGiveUp =          result["giveUp"].       as<double>();
        if( fGiveUp < 0 || fGiveUp > 1 )
        {
            std::cerr << "error: --giveUp <val>; with 0 <= <val> <= 1" << std::endl;
            return 1;
        }// else if
        auto iMatch =           result["Match"].        as<unsigned int>();
        if( iMatch == 0 )
        {
            std::cerr << "error: --Match must be larger than 0" << std::endl;
            return 1;
        }// else if
        auto iExtend =          result["Extend"].       as<unsigned int>();
        if( iExtend == 0 )
        {
            std::cerr << "error: --Extend must be larger than 0" << std::endl;
            return 1;
        }// else if
        auto iExtend2 =         result["Extend2"].      as<unsigned int>();
        auto iMissMatch =       result["MisMatch"].     as<unsigned int>();
        if( iMissMatch == 0 )
        {
            std::cerr << "error: --MisMatch must be larger than 0" << std::endl;
            return 1;
        }// else if
        auto iGap =             result["Gap"].          as<unsigned int>();
        auto iGap2 =             result["Gap2"].          as<unsigned int>();
        auto maxTries =         result["maxTries"].     as<unsigned int>();
        auto uiGenomeSizeDisable = result["minRefSize"].as<unsigned long long>();

        if(result.count("genIndex"))
        {
            std::shared_ptr<Pack> pPack(new Pack());
            //create the pack
            for(std::string sFileName : aIn)
                pPack->vAppendFASTA(sFileName.c_str());
            //store the pack
            pPack->vStoreCollection(sGenome);
            //create the fmd index
            FMIndex xFMDIndex(pPack);
            //store the fmd index
            xFMDIndex.vStoreFMIndex(sGenome.c_str());
        }//if
        else
        {
            // padding parameter is disabled at the moment
            //if(uiPadding <= 1)//input is for local
            //    std::cerr 
            //        << "WARNING: Relative padding should be larger or equal to one"
            //        << std::endl;
            if(sSeedSet != "SMEMs" && sSeedSet != "maxSpan")
                std::cerr << "WARNING: selected invalid seed set; using maxSpan" << std::endl;
            /*
            *
            * Alignment starts here
            *
            */
            //setup the alignment input
            std::shared_ptr<Pledge> pPack(new Pledge(std::shared_ptr<Container>(new Pack())));
            std::shared_ptr<Pack> pPack_(new Pack());
            pPack_->vLoadCollection(sGenome);
            pPack->set(pPack_);
            std::shared_ptr<Pledge> pFMDIndex(new Pledge(std::shared_ptr<Container>(new FMIndex())));
            std::shared_ptr<FMIndex> pFMDIndex_(new FMIndex());
            pFMDIndex_->vLoadFMIndex(sGenome);
            pFMDIndex->set(pFMDIndex_);
            std::vector<std::shared_ptr<Pledge>> aQueries;
            std::shared_ptr<Pledge> pNil(new Pledge(std::shared_ptr<Container>(new Nil())));
            pNil->set(std::shared_ptr<Container>(new Nil()));
            std::shared_ptr<FileReader> pReader;
            for(std::string sFileName : aIn)
            {
                pReader = std::shared_ptr<FileReader>(new FileReader(sFileName));
                aQueries.push_back(Module::promiseMe(
                    pReader, 
                    std::vector<std::shared_ptr<Pledge>>{pNil}
                ));
            }//for
            std::shared_ptr<Module> pOut;
            if(bFindMode)
                pOut.reset( new SeedSetFileWriter(sOut) );
            else
                pOut.reset( new FileWriter(sOut, pPack_) );
            //setup the graph
            std::vector<std::shared_ptr<Pledge>> aGraphSinks = setUpCompGraph(
                pPack,
                pFMDIndex,
                aQueries,
                pOut,
                uiReportN,
                uiT,//num threads
                bPariedNormal,
                bPariedUniform,
                uiPairedMean,
                fPairedStd,
                dPairedU,
                sSeedSet != "SMEMs",
                fGiveUp,
                iMatch,
                iMissMatch,
                iGap,
                iExtend,
                iGap2,
                iExtend2,
                bFindMode,
                std::stoi(defaults::uiMaxGapArea),
                std::stoi(defaults::uiPadding),
                maxTries,
                std::stoi(defaults::uiMinSeedSizeDrop),
                std::stoi(defaults::uiMinAmbiguity),
                std::stoi(defaults::uiMaxAmbiguity),
                uiMinLen,
                std::stoi(defaults::uiMaxEqualScoreLookahead),
                std::stof(defaults::fRelMinSeedSizeAmount),
                std::stof(defaults::fScoreDiffTolerance),
                std::stof(defaults::fMinimalQueryCoverage),
                std::stof(defaults::fScoreTolerace),
                std::stoi(defaults::uiSwitchQLen),
                defaults::bOptimisticGapEstimation == "true",
                std::stof(defaults::fSoCScoreMinimum),
                defaults::bSkipLongBWTIntervals == "true",
                std::stoi(defaults::uiCurrHarmScoreMin),
                uiGenomeSizeDisable,
                uiSoCWidth,
                bDisableHeuristics,
                fMinSecScoreRatio
            );
            // this is a hidden option
            if(result.count("info") > 0)
            {
                std::cout << "threads: " << uiT << std::endl;
                std::cout << "computational Graph:" << std::endl;
                for(auto pledge : aGraphSinks)
                    std::cout << pledge->getGraphDesc() << std::endl;
                return 0;
            }//if
            //run the alignment
            size_t uiLastProg = 0;
            std::mutex xPrintMutex;
            Pledge::simultaneousGet(
                aGraphSinks,
                [&]
                ()
                {
                    std::lock_guard<std::mutex> xGuard(xPrintMutex);
                    size_t uiCurrProg = (1000 * pReader->getCurrPosInFile()) / pReader->getFileSize();
                    if(uiCurrProg > uiLastProg)
                    {
                        std::cerr 
                            << static_cast<double>(uiCurrProg) / 10 
                            << "% aligned.     " 
                            << '\r' 
                            << std::flush;
                        uiLastProg = uiCurrProg;
                    }// if
                }// lambda
                );
        }//if
    }//try
    catch (const OptionException &ex)
    {
        std::cerr << ex.what() << std::endl;
        std::cout << sHelp << std::endl;
    }//catch
    catch (std::runtime_error &ex)
    {
        std::cerr << ex.what() << std::endl;
    }//catch
    //catch (...)
    //{
    //    std::cerr << "unknown exception encountered" << std::endl;
    //}//catch
    return 0;
}//main function