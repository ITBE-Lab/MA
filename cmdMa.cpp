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
#include "version.h"

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
    "\n        --genIndex                 Do FMD-index Generation. The -i and -x options specify the"
    "\n                                   FASTA file used for index generation and the index prefix,"
    "\n                                   respectively. If this option is not set, the aligner performs"
    "\n                                   alignments. "
    "\n"
    "\nNecessary arguments for alignments:"
    "\n    -x, --idx <prefix>             FMD-index used for alignments"
    "\n    -i, --in <fname>               FASTA or FASTAQ input files. Multiple files can be given as follows:"
    "\n                                   --in <fname1> --in <fname2>."
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
    "\n    -p, --Paired <fname>           Enable paired alignment."
    "\n                                   <val> shall be the filename of the mate reads."
    "\n        --paIsolate <num>          Penalty for an unpaired read pair."
    "\n                                   Default is 0.3."
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
    "\n                                   alignment. Generally the best alignment is found in the best"
    "\n                                   scored SoC. However, if the best alignment is from a very"
    "\n                                   repetitive region, we might have to inspect several SoC's to"
    "\n                                   find the optimal one."
    "\n                                   Default is 30."
    "\n        --minTries <num>           At least <num> many SoC's are evaluated for a single"
    "\n                                   alignment."
    "\n                                   Default is 2."
    "\n        --minRefSize <num>         If the reference is smaller than <num> nt we disable post SoC"
    "\n                                   heuristics."
    "\n                                   Default is 10,000,000."
    "\n        --noSecondary              Do not output secondary alignments."
    "\n                                   Off by default."
    "\n        --noSupplementary          Do not output supplementary alignments."
    "\n                                   Off by default."
    "\n        --maxOverlapSupp <num>     A secondary alignment becomes supplementary, if it overlaps"
    "\n                                   less than <num> percent with the primary alignment."
    "\n                                   Default is 0.1."
    "\n        --disableHeuristics        All performance optimizing heuristics are turned off."
    "\n"
    "\nVersion " MA_VERSION "\nBy Markus Schmidt & Arne Kutzner"
    "\nFor more information visit: https://github.com/ITBE-Lab/ma"
#if DEBUG_LEVEL > 0
    "\nDEBUG MODE"
#endif
    ;

void generateHelpMessage( ParameterSetManager& rManager )
{
    const std::string sHeader =
        "========================================= The Modular Aligner =========================================";
    const std::string sIndentOptions = "    ";
    std::string sIndentDesc;
    for( size_t uiI = 0; uiI < sHeader.size( ) / 2; uiI++ )
        sIndentDesc += " ";
    std::cout << sHeader << std::endl;

    // presettings
    std::cout << "Available presettings:" << std::endl;
    std::cout << sIndentOptions << "-m, --mode <";
    std::string sOptions;
    for( auto& xPair : rManager.xParametersSets )
    {
        std::string sOut = xPair.first;
        std::replace( sOut.begin( ), sOut.end( ), ' ', '_' );
        sOptions += sOut + "/";
    }
    sOptions.pop_back( );
    std::cout << sOptions << "> [" << rManager.xParametersSets.begin( )->first << "]" << std::endl << sIndentDesc;

    std::string sDescription =
        "Operation mode for MA. The selected mode will change the parameters. However, parameters you "
        "set manually will overwrite the presetting.";
    std::istringstream xStream( sDescription );
    size_t uiCharCount = 0;
    const size_t uiMaxCharCnt = sHeader.size( ) - sIndentDesc.size( );
    for( std::string sWord; xStream >> sWord; )
    {
        if( uiCharCount + sWord.size( ) >= uiMaxCharCnt )
        {
            uiCharCount = 0;
            std::cout << std::endl << sIndentDesc;
        } // if
        std::cout << sWord << " ";
        uiCharCount += sWord.size( ) + 1;
    } // for
    std::cout << std::endl << std::endl;


    // general options
    std::cout << "General options: (these options are not affected by presettings)" << std::endl;
    for( auto xTup : rManager.xGlobalParameterSet.xpAllParameters )
    {
        auto pParameter = xTup.second;
        std::cout << sIndentOptions;
        if( pParameter->cShort != pParameter->NO_SHORT_DEFINED )
            std::cout << "-" << pParameter->cShort << ", ";
        std::string sName = pParameter->sName;
        std::replace( sName.begin( ), sName.end( ), ' ', '_' );
        std::cout << "--" << sName << " <" << pParameter->type_name( ) << "> [";
        std::cout << pParameter->asText( ) << "]";
        std::cout << std::endl << sIndentDesc;

        std::istringstream xStream( pParameter->sDescription );
        size_t uiCharCount = 0;
        const size_t uiMaxCharCnt = sHeader.size( ) - sIndentDesc.size( );
        for( std::string sWord; xStream >> sWord; )
        {
            if( uiCharCount + sWord.size( ) >= uiMaxCharCnt )
            {
                uiCharCount = 0;
                std::cout << std::endl << sIndentDesc;
            } // if
            std::cout << sWord << " ";
            uiCharCount += sWord.size( ) + 1;
        } // for
        std::cout << std::endl << std::endl;
    } // for

    // other options
    for( auto xPair : rManager.getSelected( )->xpParametersByCategory )
    {
        std::cout << xPair.first.second << " options:" << std::endl;
        for( auto pParameter : xPair.second )
        {
            std::cout << sIndentOptions;
            if( pParameter->cShort != pParameter->NO_SHORT_DEFINED )
                std::cout << "-" << pParameter->cShort << ", ";
            std::string sName = pParameter->sName;
            std::replace( sName.begin( ), sName.end( ), ' ', '_' );
            std::cout << "--" << sName << " <" << pParameter->type_name( ) << "> [";
            std::cout << pParameter->asText( ) << "]";
            std::cout << std::endl << sIndentDesc;

            std::istringstream xStream( pParameter->sDescription );
            size_t uiCharCount = 0;
            const size_t uiMaxCharCnt = sHeader.size( ) - sIndentDesc.size( );
            for( std::string sWord; xStream >> sWord; )
            {
                if( uiCharCount + sWord.size( ) >= uiMaxCharCnt )
                {
                    uiCharCount = 0;
                    std::cout << std::endl << sIndentDesc;
                } // if
                std::cout << sWord << " ";
                uiCharCount += sWord.size( ) + 1;
            } // for
            std::cout << std::endl << std::endl;
        } // for
    } // for

    std::cout << "Version " << MA_VERSION << "\nBy Markus Schmidt & Arne Kutzner" << std::endl;
    std::cout << "For more information visit: https://github.com/ITBE-Lab/ma" << std::endl;
#if DEBUG_LEVEL > 0
    std::cout << "DEBUG MODE" << std::endl;
#endif
} // function

/**
 * main function
 * @todo order of parameters seems to matter.. fix that
 */
int main( int argc, char* argv[] )
{
    if( MA_VERSION != sLibMaVersion )
    {
        std::cerr << "Fatal error: cmbMA verion \"" << MA_VERSION << "\" does not match libMA version \""
                  << sLibMaVersion << "\". Something went wrong during building/linking." << std::endl;
        return 1;
    } // if
    ParameterSetManager xParameterManager;
    xParameterManager.xGlobalParameterSet.bSAMOutputInReadsFolder->set( false );

    // set the mode...
    for( int iI = 2; iI < argc; iI += 2 )
    {
        std::string sOptionName = argv[ iI - 1 ];
        std::string sOptionValue = argv[ iI ];
        std::replace( sOptionValue.begin( ), sOptionValue.end( ), '_', ' ' );
        if( sOptionName == "-m" || sOptionName == "--mode" )
            xParameterManager.setSelected( sOptionValue );
    } // for

    if( argc <= 1 )
    {
        generateHelpMessage( xParameterManager );
        return 0;
    } // if

    try
    {
        for( int iI = 1; iI < argc; iI++ )
        {
            std::string sOptionName = argv[ iI ];

            if( sOptionName == "-m" || sOptionName == "--mode" ) // we did this already
            {
                iI++; // also ignore the following argument
                continue;
            }// if

            std::replace( sOptionName.begin( ), sOptionName.end( ), '_', ' ' );

            if( iI + 1 < argc && argv[ iI + 1 ][ 0 ] != '-' )
            {
                std::string sOptionValue = argv[ iI ];
                iI++; // have key value pair so next element is certainly no key

                if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] != '-' && sOptionName.size( ) == 2 )
                    xParameterManager.byShort( sOptionName[ 1 ] )->setByText( sOptionValue );

                else if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] == '-' && sOptionName.size( ) > 2 )
                    xParameterManager.byName( sOptionName.substr( 2, sOptionName.size( ) - 2 ) )
                        ->setByText( sOptionValue );

                else
                    throw AnnotatedException(
                        std::string( "unknown option type: " )
                            .append( sOptionName )
                            .append( ". Did you forget to add the '-' or '--' at the beginning?" ) );
            } // if
            else // boolean flag option
            {
                if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] != '-' && sOptionName.size( ) == 2 )
                {
                    auto pX = std::dynamic_pointer_cast<AlignerParameter<bool>>(
                        xParameterManager.byShort( sOptionName[ 1 ] ) );

                    if( pX == nullptr )
                        throw AnnotatedException( "Parameters need to be provided as key value pairs" );
                    pX->set( true );
                } // if
                else if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] == '-' && sOptionName.size( ) > 2 )
                {
                    auto pX = std::dynamic_pointer_cast<AlignerParameter<bool>>(
                        xParameterManager.byName( sOptionName.substr( 2, sOptionName.size( ) - 2 ) ) );
                    if( pX == nullptr )
                        throw AnnotatedException( "Parameters need to be provided as key value pairs" );
                    pX->set( true );
                } // else if
                else
                    throw AnnotatedException(
                        std::string( "unknown option type: " )
                            .append( sOptionName )
                            .append( ". Did you forget to add the '-' or '--' at the beginning?" ) );
            }
        } // for
        if( xParameterManager.xGlobalParameterSet.pbPrintHelpMessage->get( ) )
            generateHelpMessage( xParameterManager );
    } // try
    catch( const AnnotatedException& ex )
    {
        std::cerr << "Error:\n" << ex.what( ) << std::endl;
    } // catch
    catch( std::runtime_error& ex )
    {
        std::cerr << "Error:\n" << ex.what( ) << std::endl;
    } // catch
    catch( std::exception& ex )
    {
        std::cerr << "Error:\n" << ex.what( ) << std::endl;
    } // catch
    catch( ... )
    {
        std::cerr << "Error:\n"
                  << "unknown exception encountered" << std::endl;
    } // catch
#if 0
    try
    {


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
            std::vector<std::shared_ptr<BasePledge>> aGraphSinks;
            std::shared_ptr<Reader> pReader;

            // setup the graph
            if( result.count( "p" ) == 0 )
            {
                std::shared_ptr<TP_WRITER> pFileWriter;
#ifdef WITH_POSTGRES
                if( sBbOutput.size( ) > 0 )
                {
                    if( iRunId == -1 )
                    {
                        DbRunConnection xConn( sBbOutput );
                        auto xRes =
                            xConn.exec( "INSERT INTO run (aligner_name, header_id) VALUES (\'MA\', 0) RETURNING "
                                        "id" );
                        iRunId = std::stoi( xRes.get( 0, 0 ) );
                    } // setupConn scope
                    pFileWriter = std::make_shared<DbWriter>( xParameterManager, sBbOutput, iRunId );
                } // if
                else
#endif
                {
                    pFileWriter = std::make_shared<FileWriter>( xParameterManager, sOut, pPack->get( ) );
                } // else or scope
                auto pFileReader = std::make_shared<FileListReader>( aIn );
                pReader = pFileReader;

                auto pQueries = promiseMe( pFileReader );
                aGraphSinks = setUpCompGraph( xParameterManager, pPack, pFMDIndex, pQueries, pFileWriter, uiT );
            } // if
            else
            {
                std::shared_ptr<TP_PAIRED_WRITER> pFileWriter;
#ifdef WITH_POSTGRES
                if( sBbOutput.size( ) > 0 )
                {
                    if( iRunId == -1 )
                    {
                        DbRunConnection xConn( sBbOutput );
                        auto xRes =
                            xConn.exec( "INSERT INTO run (aligner_name, header_id) VALUES (\'MA\', 0) RETURNING "
                                        "id" );
                        iRunId = std::stoi( xRes.get( 0, 0 ) );
                    } // setupConn scope
                    pFileWriter = std::make_shared<PairedDbWriter>( xParameterManager, sBbOutput, iRunId );
                } // if
                else
#endif
                {
                    pFileWriter = std::make_shared<PairedFileWriter>( xParameterManager, sOut, pPack->get( ) );
                } // else or scope
                std::vector<std::string> aIn2 = result[ "p" ].as<std::vector<std::string>>( );
                auto pFileReader = std::make_shared<PairedFileReader>( xParameterManager, aIn, aIn2 );
                pReader = pFileReader;

                auto pQueries = promiseMe( pFileReader );
                aGraphSinks = setUpCompGraphPaired( xParameterManager, pPack, pFMDIndex, pQueries, pFileWriter, uiT );
            } // else

            // run the alignment
            size_t uiLastProg = 0;
            size_t uiLastFile = 0;
            std::mutex xPrintMutex;

            BasePledge::simultaneousGet( aGraphSinks,
                                         [&]( ) {
                                             std::lock_guard<std::mutex> xGuard( xPrintMutex );
                                             size_t uiCurrProg =
                                                 ( 1000 * pReader->getCurrPosInFile( ) ) / pReader->getFileSize( );
                                             if( pReader->getCurrFileIndex( ) > uiLastFile || uiCurrProg > uiLastProg )
                                             {
                                                 std::cerr << "File " << pReader->getCurrFileIndex( ) << "/"
                                                           << pReader->getNumFiles( ) << ": "
                                                           << static_cast<double>( uiCurrProg ) / 10
                                                           << "% aligned.          " << '\r' << std::flush;
                                                 uiLastProg = uiCurrProg;
                                             } // if
                                             return true;
                                         } // lambda
            );
#endif
            std::cerr << "done.                " << std::endl;
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
    catch( ... )
    {
        std::cerr << "unknown exception encountered" << std::endl;
    } // catch
#endif
    return 0;
} // main function