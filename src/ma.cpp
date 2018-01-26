/** 
 * @file ma.cpp
 * @brief Parses the program options
 * @details
 * Sets up the computational graph and executes it.
 * @author Markus Schmidt
 */

#include <iostream>
#include <boost/program_options.hpp>
#include "container/nucSeq.h"
#include "container/pack.h"
#include "container/fMIndex.h"
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
    bool bPariedNormal;
    bool bPariedUniform;
    unsigned int uiPairedMean;
    double fPairedStd;
    double dPairedU;
    std::string sIndexOut;
    std::string sParameterSet;
    std::string sSeedSet;
    std::string sGenome;
    std::vector<std::string> vIndexIn;
    std::vector<std::string> aAlignOut;
    std::vector<std::string> vAlignIn;

    options_description gen_desc{"General Options"};
    gen_desc.add_options()
        ("help,h", "Complete help screen")
        ("align,a", "Sequence alignment")
        ("threads,t", value<unsigned int>(&uiT)->default_value(std::thread::hardware_concurrency()), "Used concurency")
        ("fmdIndex,f", "FMD-index generation")
    ;

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
            ("alignIn,i", value<std::vector<std::string>>(&vAlignIn)->composing(), "Input file paths [*.fasta/*.fastaq/*]")
            ("alignOut,o", value<std::vector<std::string>>(&aAlignOut)->composing(), "Output file paths [*.sam/*.bam/*]")
            ("genome,g", value<std::string>(&sGenome), "FMD-index input file prefix")
            ("parameterset,p", value<std::string>(&sParameterSet)->default_value("fast"), "Predefined parameters [fast/accurate]")
            ("maxAmbiguity,A", value<unsigned int>(&uiMaxAmbiguity)->default_value(bAccurate ? 100 : 5), "Maximal ambiguity")
            ("seedSet,s", value<std::string>(&sSeedSet)->default_value(bAccurate ? "complete" : "pairs"), "Used seed set [complete/pairs]")
            ("soc,S", value<unsigned int>(&uiNumSOC)->default_value(bAccurate ? 10 : 5), "Strip of consideration amount")
            ("global,G", "Perform global alignment")
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
            ("indexIn,I", value<std::vector<std::string>>(&vIndexIn)->composing(), "FASTA input file paths")
            ("indexOut,O", value<std::string>(&sIndexOut), "FMD-index output file prefix")
        ;
        
        options_description all_desc{"All Options"};
        all_desc.add(gen_desc).add(align_desc).add(p_desc).add(index_desc);

        variables_map vm;
        store(parse_command_line(argc, argv, all_desc), vm);
        notify(vm);

        bPariedNormal = vm.count("normal") != 0;
        bPariedUniform = vm.count("uniform") != 0;

        if (argc <= 1)
            std::cout << gen_desc << std::endl;
        if(vm.count("fmdIndex"))
        {
            if(vm.count("indexIn") == 0 || vm.count("indexOut") == 0 )
                std::cerr << "--indexIn and --indexOut are compulsory if --fmdIndex is set" << std::endl;
            else
            {
                std::shared_ptr<Pack> pPack(new Pack());
                //create the pack
                for(std::string sFileName : vIndexIn)
                    pPack->vAppendFastaFile(sFileName.c_str());
                //store the pack
                pPack->vStoreCollection(sIndexOut);
                //create the fmd index
                FMIndex xFMDIndex(pPack);
                //store the fmd index
                xFMDIndex.vStoreFMIndex(sIndexOut.c_str());
            }//else
        }//if
        if(vm.count("align"))
        {
            if(vm.count("alignIn") == 0 || vm.count("alignOut") == 0 || vm.count("genome") == 0 )
                std::cerr << "--alignIn, --alignOut and --genome are compulsory if --align is set" << std::endl;
            else
            {
                //if paired alignment
                if(bPariedNormal || bPariedUniform)
                {

                }//if
                else
                {

                }//else
            }//else
        }//if
        if (vm.count("help"))
        {
            std::cout << "\t\t===== MA THE MODULAR ALIGNER =====" << std::endl;
            std::cout << all_desc << std::endl;
            std::cout << "For more information visit: http://itbe.hanyang.ac.kr" << std::endl;
        }//if
    }//try
    catch (const error &ex)
    {
        std::cout << gen_desc << std::endl;
        std::cerr << ex.what() << std::endl;
    }//catch
    return 0;
}//main function