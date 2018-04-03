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

#include <iostream>
#include <boost/program_options.hpp>
#include "container/nucSeq.h"
#include "container/pack.h"
#include "container/fMIndex.h"
#include "module/fileReader.h"
#include "module/fileWriter.h"
#include "util/export.h"
#include <thread>

using namespace libMA;
using namespace boost::program_options;

/**
 * main function
 */
int main(int argc, char*argv[])
{
    unsigned int uiT;
    unsigned int uiMaxAmbiguity;
    unsigned int uiNumSOC;
    unsigned int uiNumNW;
    unsigned int uiReportNBest;
    nucSeqIndex uiMaxGapArea;
    bool bPariedNormal;
    bool bPariedUniform;
    unsigned int uiPairedMean;
    unsigned int uiMinAmbiguity;
    double fPairedStd;
    double dPairedU;
    std::string sIndexOut;
    std::string sParameterSet;
    std::string sSeedSet;
    std::string sGenome;
    std::vector<std::string> aIndexIn;
    std::string sAlignOut;
    std::vector<std::string> aAlignIn;
    float fGiveUp;
    float fMappingQualMin;
    float fRelativePadding;
    unsigned int iMatch;
    unsigned int iExtend;
    unsigned int iMissMatch;
    unsigned int iGap;
    float fMinGlob;
    bool bGlobal;

    options_description gen_desc{"General Options"};
    gen_desc.add_options()
        ("help,h", "Complete help screen")
        ("align,a", "Sequence alignment")
        ("threads,t", value<unsigned int>(&uiT)->default_value(std::thread::hardware_concurrency()), "Used concurency")
        ("fmdIndex,f", "FMD-index generation")
    ;

    if (argc <= 1)
    {
        std::cout << "\t\t===== MA THE MODULAR ALIGNER =====" << std::endl;
        std::cout << gen_desc << std::endl;
    }//if

    try
    {
        //precheck which parameterset is selected 
        options_description p_set_desc{"Parameter Set Options (hidden)"};
        p_set_desc.add_options()
            ("parameterset,p", value<std::string>(&sParameterSet)->default_value("fast"), "Predefined parameters [fast/accurate]")
        ;
        variables_map p_set_vm;
        parsed_options parsed = command_line_parser(argc, argv).
            options(p_set_desc).allow_unregistered().run();
        store(parsed, p_set_vm);
        notify(p_set_vm);

        //adjust the default values according to the selected parameterset
        bool bAccurate = p_set_vm["parameterset"].as<std::string>() == std::string("accurate");

        options_description align_desc{"Alignment Options"};
        align_desc.add_options()
            ("alignIn,i", value<std::vector<std::string>>(&aAlignIn)->composing(), "Input file paths [*.fasta/*.fastaq/*]")
            ("alignOut,o", value<std::string>(&sAlignOut)->default_value("stdout"), "Output file path [*.sam/*.bam/*]")
            ("reportN,n", value<unsigned int>(&uiReportNBest)->default_value(1), "Report the N best Alignments")
            ("genome,g", value<std::string>(&sGenome), "FMD-index input file prefix")
            ("parameterset,p", value<std::string>(&sParameterSet)->default_value("fast"), "Predefined parameters [fast/accurate]")
            ("maxAmbiguity,A", value<unsigned int>(&uiMaxAmbiguity)->default_value(bAccurate ? 1000 : 10), "Maximal ambiguity")
            ("minAmbiguity,B", value<unsigned int>(&uiMinAmbiguity)->default_value(bAccurate ? 3 : 0), "Minimal ambiguity")
            ("seedSet,s", value<std::string>(&sSeedSet)->default_value(bAccurate ? "SMEMs" : "maxSpanning"), "Used seed set [SMEMs/maxSpanning]")
            ("soc,S", value<unsigned int>(&uiNumSOC)->default_value(bAccurate ? 10 : 3), "Strip of consideration amount")
            ("numNw,w", value<unsigned int>(&uiNumNW)->default_value(bAccurate ? 10 : 3), "apply NW to x SoC")
            ("nwLimit,l", value<nucSeqIndex>(&uiMaxGapArea)->default_value(1000000), "Maximal DP matrix size")
            ("giveUp,v", value<float>(&fGiveUp)->default_value(bAccurate ? 0.03 : 0.01), "Give up if the best SoC score is lower")
            ("toGlobal,T", value<float>(&fMappingQualMin)->default_value(bAccurate ? 0 : 0.05), "Transform lower mapping qual local alignments to global")
            ("padding,P", value<float>(&fRelativePadding)->default_value(1.1), "Relative padding for global alignments.")
            ("Match", value<unsigned int>(&iMatch)->default_value(2), "NW match score.")
            ("MissMatch", value<unsigned int>(&iMissMatch)->default_value(4), "NW missmatch penalty.")
            ("Gap", value<unsigned int>(&iGap)->default_value(6), "NW gap open penalty.")
            ("Extend", value<unsigned int>(&iExtend)->default_value(1), "NW gap extend penalty.")
            ("socMinScore", value<float>(&fMinGlob)->default_value(-2.0), "Minimum score for SoC width.")
            ("global,G", value<bool>(&bGlobal)->default_value(bAccurate ? true : false), "Perform global alignment")
        ;

        options_description p_desc{"Paired Reads Options"};
        p_desc.add_options()
            ("uniform,U", "Enable paired alignment; Gaps modeled as uniform distribution")
            ("normal,N", "Enable paired alignment; Gaps modeled as normal distribution")
            ("unpaired,u", value<double>(&dPairedU)->default_value(17), "Penalty for unpaired alignments")
            ("mean,m", value<unsigned int>(&uiPairedMean)->default_value(0), "Gap distance mean")
            ("std,d", value<double>(&fPairedStd)->default_value(1.0), "Gap distance standard deviation")
        ;

        options_description index_desc{"FMD-Index Generation Options"};
        index_desc.add_options()
            ("indexIn,I", value<std::vector<std::string>>(&aIndexIn)->composing(), "FASTA input file paths")
            ("indexOut,O", value<std::string>(&sIndexOut), "FMD-index output file prefix")
        ;

        options_description hidden_desc{"Hidden Options"};
        hidden_desc.add_options()
            ("info", "Print comp. graph info & abort")
        ;
        
        options_description all_desc{"All Options"};
        all_desc.add(gen_desc).add(align_desc).add(p_desc).add(index_desc);

        options_description comp_desc{"Complete Options"};
        comp_desc.add(gen_desc).add(align_desc).add(p_desc).add(index_desc).add(hidden_desc);

        variables_map vm;
        store(parse_command_line(argc, argv, comp_desc), vm);
        notify(vm);


        bPariedNormal = vm.count("normal") != 0;
        bPariedUniform = vm.count("uniform") != 0;
        bool bDoneSth = false;

        if (vm.count("help"))
        {
            std::cout << "\t\t===== MA THE MODULAR ALIGNER =====" << std::endl;
            std::cout << all_desc << std::endl;
            std::cout << "For more information visit: http://itbe.hanyang.ac.kr" << std::endl;
            bDoneSth = true;
        }//if
        if(vm.count("fmdIndex"))
        {
            if(vm.count("indexIn") == 0 || vm.count("indexOut") == 0 )
                std::cerr << "--indexIn and --indexOut are compulsory if --fmdIndex is set" << std::endl;
            else
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
            }//else
            bDoneSth = true;
        }//if
        if(vm.count("align"))
        {
            if(vm.count("alignIn") == 0 || vm.count("genome") == 0 )
                std::cerr << "--alignIn and --genome are compulsory if --align is set" << std::endl;
            else
            {
                if(fRelativePadding <= 1)//input is for local
                    std::cerr 
                        << "WARNING: Relative padding should be larger or equal to one"
                        << std::endl;
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
                std::shared_ptr<Module> pOut(new FileWriter(sAlignOut));
                //setup the graph
                std::vector<std::shared_ptr<Pledge>> aGraphSinks = setUpCompGraph(
                    pPack,
                    pFMDIndex,
                    aQueries,
                    pOut,
                    uiT,//num threads
                    uiMaxAmbiguity,
                    uiNumSOC,
                    bPariedNormal,
                    bPariedUniform,
                    uiPairedMean,
                    fPairedStd,
                    dPairedU,
                    sSeedSet != "SMEMs",
                    uiReportNBest,
                    !bGlobal,//input is for local
                    uiMaxGapArea,
                    iMatch,
                    iMissMatch,
                    iGap,
                    iExtend,
                    uiMinAmbiguity,
                    fMinGlob, // socMinimal score
                    uiNumNW,
                    fGiveUp,
                    fMappingQualMin,
                    fRelativePadding
                );
                if(vm.count("info") > 0)
                {
                    std::cout << "computational Graph:" << std::endl;
                    for(auto pledge : aGraphSinks)
                        std::cout << pledge->getGraphDesc() << std::endl;
                    return 0;
                }//if
                //run the alignment
                Pledge::simultaneousGet(aGraphSinks, true);
                bDoneSth = true;
            }//else
        }//if
        if(!bDoneSth)
            std::cout << "No task was specified. Use one of the following: -h, -a, -f" << std::endl;
    }//try
    catch (const error &ex)
    {
        std::cout << gen_desc << std::endl;
        std::cerr << ex.what() << std::endl;
    }//catch
    catch (std::runtime_error &ex)
    {
        std::cerr << ex.what() << std::endl;
    }//catch
    catch (...)
    {
        std::cerr << "unknown exception encountered" << std::endl;
    }//catch
    return 0;
}//main function