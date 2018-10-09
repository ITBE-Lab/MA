/**
 * @file cmdMa.cpp
 * @brief Parses the program options
 * @details
 * Sets up the computational graph and executes it.
 * @author Markus Schmidt
 * @copyright
Copyright 2018 Markus Schmidt, Arne Kutzner

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "container/fMIndex.h"
#include "container/nucSeq.h"
#include "container/pack.h"
#include "module/dbWriter.h"
#include "module/fileReader.h"
#include "module/fileWriter.h"
#include "util/cxxopts.hpp"
#include "util/debug.h"
#include "util/default_parameters.h"
#include "util/export.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <iostream>
#include <string.h>
#include <thread>
/// @endcond

using namespace libMA;
using namespace cxxopts;

const std::string sHelp =
    "====================================== The Modular Aligner ======================================"
    "\nGeneral options:"
    "\n    -h, --help                     Display the complete help screen"
    "\n        --genIndex                 Do FMD-index Generation. The -i and -x options specify "
    "the"
    "\n                                   FASTA file used for index generation and the index "
    "prefix,"
    "\n                                   respectively. If this option is not set, the aligner "
    "performs"
    "\n                                   alignments. "
    "\n"
    "\nNecessary arguments for alignments:"
    "\n    -x, --idx <prefix>             FMD-index used for alignments"
    "\n    -i, --in <fname>               FASTA or FASTAQ input files."
    "\n"
    "\nAvailable presettings:"
    "\n     -m, --mode [str]              Operation mode for MA. (default is 'fast')"
    "\n                                   fast:"
    "\n                                       Best compromise between performance and accuracy."
    "\n                                       Recommended for PacBio reads."
    "\n                                   acc:"
    "\n                                       Better accuracy than in fast mode but worse runtimes."
    "\n                                       Particularly effective for short reads."
    "\n"
    "\nAlignments options:"
    "\n    -o, --out <fname>              Filename used for SAM file output. Default output stream is"
    "\n                                   standard output."
    "\n    -t, --threads <num>            Use <num> threads. On startup MA checks the hardware and "
    "\n                                   chooses this value accordingly."
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
    "\n                                   --in <fname1> --in <fname2>."
    "\n    -P, --paNorm                   Enable paired alignment and Model the distance as"
    "\n                                   normal distribution."
    "\n                                   If set --in shall be used as follows:"
    "\n                                   --in <fname1> --in <fname2>."
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
    "\n        --maxTries <num>           Maximally <num> many SoC's are evaluated for a single"
    "\n                                   alignment. Generally the best alignment is found in the "
    "best"
    "\n                                   scored SoC. However, if the best alignment is from a very"
    "\n                                   repetitive region, we might have to inspect several "
    "SoC's to"
    "\n                                   find the optimal one."
    "\n                                   Default is 30."
    "\n        --minTries <num>           At least <num> many SoC's are evaluated for a single"
    "\n                                   alignment."
    "\n                                   Default is 2."
    "\n        --minRefSize <num>         If the reference is smaller than <num> nt we disable "
    "post SoC"
    "\n                                   heuristics."
    "\n                                   Default is 10,000,000."
    "\n        --minSecToPrimRatio <num>  Limit output of secondary alignments to alignments with:"
    "\n                                   score-secondary-alignment >= <num> * "
    "score-primary-alignment."
    "\n                                   Default is 0.25."
    "\n        --maxOverlapSupp <num>     A secondary alignment becomes supplementary, if it "
    "overlaps"
    "\n                                   less than <num> percent with the primary alignment."
    "\n                                   Default is 0.1."
    "\n        --disableHeuristics        All performance optimizing heuristics are turned off."
    "\n"
    "\nVersion 0.1.2 (Sep 2018)"
    "\nBy Markus Schmidt & Arne Kutzner"
    "\nFor more information visit: https://github.com/ITBE-Lab/ma";

/**
 * main function
 * @todo order of parameters seems to matter.. fix that
 */
int main( int argc, char* argv[] )
{
    Options options( "MA", "\t\t===== The Modular Aligner =====" );
    options.add_options( )( "h,help", "Display the complete help screen" )(
        "t,threads", "Number of threads",
        value<unsigned int>( )->default_value( std::to_string( std::thread::hardware_concurrency( ) ) ), "arg     " );

    if( argc <= 1 )
    {
        std::cout << "Invalid usage; you need to supply -x and -i at least!" << std::endl;
        std::cout << sHelp << std::endl;
        return 0;
    } // if

    try
    {
        defaults::configureFast( );
        for( int i = 0; i < argc - 1; i++ )
            if( strcmp( argv[ i ], "-m" ) == 0 || strcmp( argv[ i ], "--mode" ) == 0 )
            {
                if( strcmp( argv[ i + 1 ], "acc" ) == 0 )
                    defaults::configureAccurate( );
                if( strcmp( argv[ i + 1 ], "fast" ) == 0 )
                    defaults::configureFast( );
            } // if

        options.add_options( "Alignment options (requires -a)" )( "i,in", "Input file(s) as (multi-)fasta(-q)",
                                                                  value<std::vector<std::string>>( ), "args" )(
            "o,out", "Output file as SAM", value<std::string>( )->default_value( "stdout" ) )(
            "m,mode", "Pre-setting [fast/acc]", value<std::string>( )->default_value( defaults::sParameterSet ) )(
            "s,seedSet", "Seeding strategy [SMEMs/maxSpanning]",
            value<std::string>( )->default_value( defaults::sSeedSet ) )(
            "n,reportN", "Report up to N alignments; 0: unlimited",
            value<unsigned int>( )->default_value( std::to_string( defaults::uiReportN ) ) )(
            "l,minLen", "Minimum seed length",
            value<unsigned int>( )->default_value( std::to_string( defaults::uiMinLen ) ) )(
            "giveUp", "Minimum SoC score (relative to query length)",
            value<double>( )->default_value( std::to_string( defaults::fGiveUp ) ) )(
            "minSecToPrimRatio", "Minimum sec score ratio",
            value<double>( )->default_value( std::to_string( defaults::fMinSecScoreRatio ) ) )(
            "maxTries", "Max num SoC",
            value<unsigned int>( )->default_value( std::to_string( defaults::uiMaxTries ) ) )(
            "minTries", "Min num SoC",
            value<unsigned int>( )->default_value( std::to_string( defaults::uiMinTries ) ) )(
            "minRefSize", "ref size switch",
            value<unsigned long long>( )->default_value( std::to_string( defaults::uiGenomeSizeDisable ) ) )(
            "d,noDP", "Disable DP" )( "Match", "DP match score.",
                                      value<int>( )->default_value( std::to_string( defaults::iMatch ) ) )(
            "MisMatch", "DP mismatch penalty.",
            value<int>( )->default_value( std::to_string( defaults::iMissMatch ) ) )(
            "Gap", "DP gap open penalty.", value<int>( )->default_value( std::to_string( defaults::iGap ) ) )(
            "Extend", "DP gap extend penalty.", value<int>( )->default_value( std::to_string( defaults::iExtend ) ) )(
            "Gap2", "DP gap open penalty.", value<int>( )->default_value( std::to_string( defaults::iGap2 ) ) )(
            "Extend2", "DP gap extend penalty.", value<int>( )->default_value( std::to_string( defaults::iExtend2 ) ) )(
            "x,idx", "Do FMD-index generation", value<std::string>( ) )( "genIndex", "Do FMD-index generation" )(
            "SoCWidth", "SoC width", value<unsigned int>( )->default_value( std::to_string( defaults::uiSoCWidth ) ) )(
            "disableHeuristics", "disable all heuristics",
            value<bool>( )->default_value( defaults::bDisableHeuristics ? "true" : "false" ) )(
            "maxDeltaDist", "", value<double>( )->default_value( std::to_string( defaults::dMaxDeltaDist ) ) )(
            "minDeltaDist", "", value<uint64_t>( )->default_value( std::to_string( defaults::uiMinDeltaDist ) ) )(
            "maxOverlapSupp", "",
            value<double>( )->default_value( std::to_string( defaults::dMaxOverlapSupplementary ) ) );

        options.add_options( "Paired Reads options (requires either -U or -N)" )(
            "p,paUni", "Enable paired alignment; Distance as uniform distribution" )(
            "P,paNorm", "Enable paired alignment; Distance as normal distribution" )(
            "paIsolate", "Penalty for unpaired alignments",
            value<unsigned int>( )->default_value( std::to_string( defaults::uiUnpaired ) ),
            "arg    " )( "paMean", "Gap distance mean",
                         value<unsigned int>( )->default_value( std::to_string( defaults::uiMean ) ) )(
            "paStd", "Gap distance standard deviation",
            value<double>( )->default_value( std::to_string( defaults::fStd ) ) )(
            "db_conninfo", "db_conninfo", value<std::string>( )->default_value( "" ) )(
            "run_id", "run_id", value<int32_t>( )->default_value( "-1" ) );

        auto result = options.parse( argc, argv );

        if( result.count( "help" ) )
        {
            std::cout << sHelp << std::endl;
            //@todo cmake version number:
            // https://stackoverflow.com/questions/27395120/correct-way-to-encode-embed-version-number-in-program-code
            DEBUG( std::cout << "DEBUG LEVEL: " << DEBUG_LEVEL << std::endl; )
            return 0;
        } // if

        auto uiT = result[ "threads" ].as<unsigned int>( );
        defaults::bNormalDist = result.count( "paNorm" ) > 0;
        defaults::bUniformDist = result.count( "paUni" ) > 0;
        defaults::uiMean = result[ "paMean" ].as<unsigned int>( );
        defaults::fStd = result[ "paStd" ].as<double>( );
        defaults::uiUnpaired = result[ "paIsolate" ].as<unsigned int>( );
        defaults::fMinSecScoreRatio = result[ "minSecToPrimRatio" ].as<float>( );
        defaults::uiReportN = result[ "reportN" ].as<unsigned int>( );
        defaults::sParameterSet = result[ "mode" ].as<std::string>( );
        defaults::sSeedSet = result[ "seedSet" ].as<std::string>( );
        defaults::uiMinLen = result[ "minLen" ].as<unsigned int>( );
        defaults::uiSoCWidth = result[ "SoCWidth" ].as<unsigned int>( );
        defaults::bDisableHeuristics = result[ "disableHeuristics" ].as<bool>( );
        defaults::uiMinDeltaDist = result[ "minDeltaDist" ].as<uint64_t>( );
        defaults::dMaxDeltaDist = result[ "maxDeltaDist" ].as<double>( );
        defaults::dMaxOverlapSupplementary = result[ "maxOverlapSupp" ].as<double>( );
        if( defaults::bNormalDist && defaults::bUniformDist )
        {
            std::cerr << "--normal and --uniform are exclusive." << std::endl;
            return 1;
        } // else if
        std::string sGenome;
        if( result.count( "idx" ) > 0 )
            sGenome = result[ "idx" ].as<std::string>( );
        else
        {
            std::cerr << "error: --idx is compulsory" << std::endl;
            return 1;
        } // else if
        std::vector<std::string> aIn;
        if( result.count( "in" ) > 0 )
            aIn = result[ "in" ].as<std::vector<std::string>>( );
        else
        {
            std::cerr << "error: --in is compulsory" << std::endl;
            return 1;
        } // else if
        if( aIn.size( ) != 1 && !( defaults::bNormalDist || defaults::bUniformDist ) &&
            result.count( "genIndex" ) == 0 )
        {
            std::cerr << "error: --in takes one argument in unpaired mode" << std::endl;
            return 1;
        } // if
        else if( aIn.size( ) != 2 && ( defaults::bNormalDist || defaults::bUniformDist ) &&
                 result.count( "genIndex" ) == 0 )
        {
            std::cerr << "error: --in takes two arguments in paired mode" << std::endl;
            return 1;
        } // else if
        auto sOut = result[ "out" ].as<std::string>( );
        defaults::fGiveUp = result[ "giveUp" ].as<float>( );
        if( defaults::fGiveUp < 0 || defaults::fGiveUp > 1 )
        {
            std::cerr << "error: --giveUp <val>; with 0 <= <val> <= 1" << std::endl;
            return 1;
        } // else if
        defaults::iMatch = result[ "Match" ].as<int>( );
        if( defaults::iMatch == 0 )
        {
            std::cerr << "error: --Match must be larger than 0" << std::endl;
            return 1;
        } // else if
        defaults::iExtend = result[ "Extend" ].as<int>( );
        if( defaults::iExtend == 0 )
        {
            std::cerr << "error: --Extend must be larger than 0" << std::endl;
            return 1;
        } // else if
        defaults::iExtend2 = result[ "Extend2" ].as<int>( );
        defaults::iMissMatch = result[ "MisMatch" ].as<int>( );
        if( defaults::iMissMatch == 0 )
        {
            std::cerr << "error: --MisMatch must be larger than 0" << std::endl;
            return 1;
        } // else if
        defaults::iGap = result[ "Gap" ].as<int>( );
        defaults::iGap2 = result[ "Gap2" ].as<int>( );
        defaults::uiMaxTries = result[ "maxTries" ].as<unsigned int>( );
        defaults::uiMinTries = result[ "minTries" ].as<unsigned int>( );
        defaults::uiGenomeSizeDisable = result[ "minRefSize" ].as<unsigned long long>( );
        std::string sBbOutput = result[ "db_conninfo" ].as<std::string>( );
#ifdef WITH_POSTGRES
        //int32_t iRunId = result[ "run_id" ].as<int32_t>( );
#endif

        if( result.count( "genIndex" ) )
        {
            std::shared_ptr<Pack> pPack( new Pack( ) );
            // create the pack
            for( std::string sFileName : aIn )
                pPack->vAppendFASTA( sFileName.c_str( ) );
            // store the pack
            pPack->vStoreCollection( sGenome );
            // create the fmd index
            FMIndex xFMDIndex( pPack );
            // store the fmd index
            xFMDIndex.vStoreFMIndex( sGenome.c_str( ) );
        } // if
        else
        {
            // padding parameter is disabled at the moment
            // if(uiPadding <= 1)//input is for local
            //    std::cerr
            //        << "WARNING: Relative padding should be larger or equal to one"
            //        << std::endl;
            if( defaults::sSeedSet != "SMEMs" && defaults::sSeedSet != "maxSpan" )
                std::cerr << "WARNING: selected invalid seed set; using maxSpan" << std::endl;
/*
 *
 * Alignment starts here
 *
 */
// setup the alignment input
#if 1
            auto pPack = makePledge<Pack>( sGenome );
            auto pFMDIndex = makePledge<FMIndex>( sGenome );

            auto pFileReader = std::make_shared<FileReader>( aIn[ 0 ] );
            std::shared_ptr<WriterModule> pFileWriter = std::static_pointer_cast<WriterModule>(
				std::make_shared<FileWriter>( sOut, pPack->get( ) ));

            auto pQueries = promiseMe( pFileReader );
            // if( aIn.size( ) == 1 )
            // {
            //     pReader = std::shared_ptr<FileReader>( new FileReader( aIn[ 0 ] ) );
            //     pQueries = Module::promiseMe( pReader, std::vector<std::shared_ptr<Pledge>>{pNil} );
            // } // if
            // else if( aIn.size( ) == 2 )
            // {
            //     pReader = std::shared_ptr<PairedFileReader>( new PairedFileReader( aIn[ 0 ], aIn[ 1 ] ) );
            //     pQueries = Module::promiseMe( pReader, std::vector<std::shared_ptr<Pledge>>{pNil} );
            // }
            // else
            // {
            //     throw AnnotatedException( "Cannot have more than two inputs!" );
            // } // else
            // std::vector<std::shared_ptr<Module>> vOut;
#ifdef WITH_POSTGRES
            // if( sBbOutput.size( ) > 0 )
            // {
            //     if( iRunId == -1 )
            //     {
            //         DbRunConnection xConn( sBbOutput );
            //         auto xRes = xConn.exec( "INSERT INTO run (aligner_name, header_id) VALUES (\'MA\', 0) RETURNING "
            //                                 "id" );
            //         iRunId = std::stoi( xRes.get( 0, 0 ) );
            //     } // setupConn scope
            //     for( size_t uiI = 0; uiI < uiT; uiI++ )
            //         vOut.emplace_back( new DbWriter( sBbOutput, iRunId ) );
            // } // if
            // else
#endif
            // {
            //     vOut.emplace_back( new FileWriter( sOut, pPack_ ) );
            // } // else
            // bool bPaired = defaults::bNormalDist || defaults::bUniformDist;
            // std::vector<std::shared_ptr<Pledge>> aGraphSinks;
            // if( bPaired )
            //    aGraphSinks = setUpCompGraphPaired( pPack, pFMDIndex, pQueries, vOut,
            //                                        uiT // num threads
            //    );
            // else
            // setup the graph
            auto aGraphSinks = setUpCompGraph( pPack, pFMDIndex, pQueries, pFileWriter, uiT );
            // run the alignment
            size_t uiLastProg = 0;
            std::mutex xPrintMutex;
            BasePledge::simultaneousGet( aGraphSinks,
                                         [&]( ) {
                                             std::lock_guard<std::mutex> xGuard( xPrintMutex );
                                             size_t uiCurrProg = ( 1000 * pFileReader->getCurrPosInFile( ) ) /
                                                                 pFileReader->getFileSize( );
                                             if( uiCurrProg > uiLastProg )
                                             {
                                                 std::cerr << " " << static_cast<double>( uiCurrProg ) / 10
                                                           << "% aligned.     " << '\r' << std::flush;
                                                 uiLastProg = uiCurrProg;
                                             } // if
                                         } // lambda
            );
#endif
            std::cerr << "100% aligned.     " << std::endl;
        } // if
    } // try
    catch( const OptionException& ex )
    {
        std::cerr << ex.what( ) << std::endl;
        std::cout << sHelp << std::endl;
    } // catch
    catch( std::runtime_error& ex )
    {
        std::cerr << ex.what( ) << std::endl;
    } // catch
    catch( std::exception& ex )
    {
        std::cerr << ex.what( ) << std::endl;
    } // catch
    // catch (...)
    //{
    //    std::cerr << "unknown exception encountered" << std::endl;
    //}//catch
    return 0;
} // main function