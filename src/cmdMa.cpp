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
"""
=========================================== The Modular Aligner ===========================================
General options:
    -h, --help                     Display the complete help screen
        --genIndex                 Do FMD-index Generation. The -i and -x options specify the FASTA
                                   file used for index generation and the index prefix, respectively.
                                   If this option is not set, the aligner performs alignments. 

Necessary arguments for alignments:
    -x, --idx <prefix>             FMD-index used for alignments
    -i, --in <fname>               FASTA or FASTAQ input files.

Alignments options:
    -o, --out <fname>              Filename used for SAM file output. Default output stream is
                                   standard output.
    -t, --threads <num>            Use <num> threads. On startup MA checks the hardware and chooses 
                                   this value accordingly.
    -m, --mode [fast/acc]           Set operation modus for MA. 
                                   Default is 'fast'.
    -d, --noDP                     Switch that disables the final Dynamic Programming.
    -n, --reportN <num>            Report up to <num> alignments; 0 means unlimited.
                                   Default is 1.
    -s, --seedSet [SMEMs/maxSpan]  Selects between the two seeding strategies super maximal extend matches
                                   'SMEMs' and maximally spanning seeds 'maxSpan'. 
                                   Default is 'maxSpan'.
    -l, --minLen <num>             Seeds must have a minimum length of <num> nucleotides.
                                   Default is 16.
        --Match <num>              Sets the match score to <num>; <num> > 0.
                                   Default is 3. 
        --MissMatch <num>          Sets the mismatch penalty to <num>; <num> > 0.
                                   Default is 4.
        --Gap <num>                Sets the costs for opening a gap to <num>; <num> >= 0.
                                   Default is 6.
        --Extend <num>             Sets the costs for extending a gap to <num>; <num> > 0.
                                   Default is 1

Paired Reads options:
    -p, --paUni                    Enable paired alignments and model the distance as uniform distribution.
                                   If enabled --in shall be used as follows: --in '<fname1>, <fname2>'.
    -P, --paNorm                   Enable paired alignment and Model the distance as normal distribution.
                                   If enabled --in shall be used as follows: '--in <fname1>, <fname2>'.
        --paIsolate <num>          Penalty for an unpaired read pair.
                                   Default is 17.
        --paMean <num>             Mean gap distance between read pairs.
                                   Default is 400.
        --paStd <num>              Standard deviation of gap distance between read pairs.
                                   Default is 150.

Advanced options:
        --giveUp <val>             Threshold with 0 <= <val> <= 1 used as give-up criteria.
                                   SoC's with accumulative seed length smaller than 
                                   'query_len * <val>' will be ignored.
                                   Reducing this parameter will decrease runtime, but allow
                                   the aligner to discover more dissimilar matches.
                                   Increasing this parameter will increase runtime, but might cause
                                   the aligner to miss the correct reference location.
                                   Default is 0.002.
        --maxTries <num>           At most the best <num> SoC's will be inspected.
                                   Generally the best alignment is found in the best scored SoC.
                                   However, if the best alignment is from a very repetitive region,
                                   we might have to inspect several SoC's to find the optimal one.
                                   Default is 50.
        --minRefSize <num>         If the reference is smaller than <num> nt we disable all heuristics.
                                   Default is 10000000.

Version 0.1.0 (alpha)
By Markus Schmidt & Arne Kutzner
For more information visit: https://github.com/ITBE-Lab/ma
""";

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
        ("x,idx", "Do FMD-index generation")
    ;

    if (argc <= 1)
    {
        std::cout << options.help({"", "General"}) << std::endl;
        return 0;
    }//if

    try
    {
        defaults::configureFast();
        for(int i = 0; i < argc-1; i++)
            if(
                    strcmp(argv[i], "-p") == 0 ||
                    strcmp(argv[i], "--parameterSet") == 0
                )
            {
                if(strcmp(argv[i+1], "accurate") == 0)
                    defaults::configureAccurate();
                if(strcmp(argv[i+1], "fast") == 0)
                    defaults::configureFast();
            }// if

        options.add_options("Alignment options (requires -a)")
            ("i,in", "Input file(s) as (multi-)fasta(-q)", 
                value<std::vector<std::string>>(), "args"
            )
            ("o,out", "Output file as SAM",value<std::string>()->default_value("stdout"))
            ("m,mode", "Pre-setting [fast/accurate]", 
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
            ("maxTries", "Max num SoC",
                value<unsigned int>()->default_value(defaults::uiMaxTries)
            )
            ("minRefSize", "ref size switch",
                value<unsigned long long>()->default_value(defaults::uiGenomeSizeDisable)
            )
            ("d,noDP", "Disable DP", value<bool>()->default_value(defaults::bFindMode))
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
        ;

        options.add_options("Paired Reads options (requires either -U or -N)")
            ("p,paUni", "Enable paired alignment; Distance as uniform distribution")
            ("P,paNormal", "Enable paired alignment; Distance as normal distribution")
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


        auto uiT =              result["threads"].      as<unsigned int>();
        auto bPariedNormal =    result.count("normal")  > 0;
        auto bPariedUniform =   result.count("uniform") > 0;
        auto uiPairedMean =     result["mean"].         as<unsigned int>();
        auto fPairedStd =       result["std"].          as<double>();
        auto dPairedU =         result["unpaired"].     as<double>();
        auto uiReportN =        result["reportN"].      as<unsigned int>();
        auto sParameterSet =    result["parameterSet"]. as<std::string>();
        auto sSeedSet =         result["seedSet"].      as<std::string>();
        auto uiMinLen =         result["minLen"].       as<unsigned int>();
        if(bPariedNormal && bPariedUniform)
        {
            std::cerr << "--normal and --uniform are exclusive." << std::endl;
            return 1;
        }// else if
        std::string sGenome;
        if( result.count("genome") > 0 )
            sGenome =           result["genome"].       as<std::string>();
        else if(result.count("align") > 0)
        {
            std::cerr << "error: --genome is compulsory if --align is set" << std::endl;
            return 1;
        }// else if
        std::vector<std::string> aIndexIn;
        if( result.count("indexIn") > 0 )
            aIndexIn =          result["indexIn"].      as<std::vector<std::string>>();
        else if(result.count("fmdIndex")  > 0)
        {
            std::cerr << "error: --indexIn is compulsory if --fmdIndex is set" << std::endl;
            return 1;
        }// else if
        std::string sIndexOut;
        if( result.count("indexOut") > 0 )
            sIndexOut =         result["indexOut"].     as<std::string>();
        else if(result.count("fmdIndex")  > 0)
        {
            std::cerr << "error: --indexOut is compulsory if --fmdIndex is set" << std::endl;
            return 1;
        }// else if
        auto sAlignOut =        result["alignOut"].     as<std::string>();
        std::vector<std::string> aAlignIn;
        if( result.count("alignIn") > 0 )
            aAlignIn =          result["alignIn"].      as<std::vector<std::string>>();
        else if(result.count("align")  > 0)
        {
            std::cerr << "error: --alignIn is compulsory if --align is set" << std::endl;
            return 1;
        }// else if
        auto bFindMode =        result.count("basicMode") > 0;
        auto fGiveUp =          result["giveUp"].       as<double>();
        auto iMatch =           result["Match"].        as<unsigned int>();
        auto iExtend =          result["Extend"].       as<unsigned int>();
        auto iMissMatch =       result["MisMatch"].     as<unsigned int>();
        auto iGap =             result["Gap"].          as<unsigned int>();
        auto maxTries =         result["maxTries"].     as<unsigned int>();
        auto uiGenomeSizeDisable = result["minRefSize"].as<unsigned long long>();

        bool bDoneSth = false;

        if (result.count("help"))
        {
            std::cout << sHelp << std::endl;
            //@todo cmake version number: https://stackoverflow.com/questions/27395120/correct-way-to-encode-embed-version-number-in-program-code
            DEBUG(
                std::cout << "DEBUG LEVEL: " << DEBUG_LEVEL << std::endl;
            )
            bDoneSth = true;
        }//if
        if(result.count("fmdIndex"))
        {
            std::shared_ptr<Pack> pPack(new Pack());
            //create the pack
            for(std::string sFileName : aIndexIn)
                pPack->vAppendFASTA(sFileName.c_str());
            //store the pack
            pPack->vStoreCollection(sIndexOut);
            //create the fmd index
            FMIndex xFMDIndex(pPack);
            //store the fmd index
            xFMDIndex.vStoreFMIndex(sIndexOut.c_str());
            bDoneSth = true;
        }//if
        if(result.count("align"))
        {
            // padding parameter is disabled at the moment
            //if(uiPadding <= 1)//input is for local
            //    std::cerr 
            //        << "WARNING: Relative padding should be larger or equal to one"
            //        << std::endl;
            if(sSeedSet != "SMEMs" && sSeedSet != "maxSpanning")
                std::cerr << "WARNING: selected invalid seed set; using maxSpanning" << std::endl;
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
            for(std::string sFileName : aAlignIn)
            {
                std::shared_ptr<FileReader> pReader(new FileReader(sFileName));
                aQueries.push_back(Module::promiseMe(
                    pReader, 
                    std::vector<std::shared_ptr<Pledge>>{pNil}
                ));
            }//for
            std::shared_ptr<Module> pOut;
            if(bFindMode)
                pOut.reset( new SeedSetFileWriter(sAlignOut) );
            else
                pOut.reset( new FileWriter(sAlignOut) );
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
                uiGenomeSizeDisable
            );
            if(result.count("info") > 0)
            {
                std::cout << "threads: " << uiT << std::endl;
                std::cout << "computational Graph:" << std::endl;
                for(auto pledge : aGraphSinks)
                    std::cout << pledge->getGraphDesc() << std::endl;
                return 0;
            }//if
            //run the alignment
            Pledge::simultaneousGet(aGraphSinks, true);
            bDoneSth = true;
        }//if
        if(!bDoneSth)
            std::cout << "No task was specified. Use one of the following: -h, -a, -f" << std::endl;
    }//try
    catch (const OptionException &ex)
    {
        std::cout << options.help({"", "General Options"}) << std::endl;
        std::cerr << ex.what() << std::endl;
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