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
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#endif

#include "container/fMIndex.h"
#include "container/nucSeq.h"
#include "container/pack.h"
#include "module/fileReader.h"
#include "module/fileWriter.h"
#include "debug.h"
#include "util/execution-context.h"
#include "util/export.h"
#include "version.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <iostream>
#include <string.h>
#include <thread>
/// @endcond

using namespace libMA;

const std::string sHeader =
    "========================================= The Modular Aligner =========================================";
const std::string sIndentOptions = "    ";
void printOption( std::string sName,
                  const char cShort,
                  const std::string& sTypeName,
                  const std::string& sDefaultVal,
                  const std::string& sDescription,
                  const std::string& sIndentDesc )
{
    std::string sOptionHead = sIndentOptions;
    if( cShort != AlignerParameterBase::NO_SHORT_DEFINED )
    {
        sOptionHead.append( "-" );
        sOptionHead.push_back( cShort );
        sOptionHead.append( ", " );
    } // if
    std::replace( sName.begin( ), sName.end( ), ' ', '_' );
    sOptionHead.append( "--" );
    sOptionHead.append( sName );
    sOptionHead.append( " <" );
    sOptionHead.append( sTypeName );
    sOptionHead.append( "> [" );
    sOptionHead.append( sDefaultVal );
    sOptionHead.append( "]" );
    if( sOptionHead.size( ) < sIndentDesc.size( ) - 4 )
    {
        std::cout << sOptionHead;
        for( size_t i = sOptionHead.size( ); i < sIndentDesc.size( ); i++ )
            std::cout << " ";
    } // if
    else
        std::cout << sOptionHead << std::endl << sIndentDesc;

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
} // function

void generateHelpMessage( ParameterSetManager& rManager, bool bFull = true )
{
    std::string sIndentDesc;
    for( size_t uiI = 0; uiI < sHeader.size( ) / 2; uiI++ )
        sIndentDesc += " ";
    std::cout << sHeader << std::endl;

    // presettings

    std::cout << "Available presettings:" << std::endl;
    std::string sOptions = "'";
    for( auto& xPair : rManager.xParametersSets )
    {
        std::string sOut = xPair.second->sName;
        std::replace( sOut.begin( ), sOut.end( ), ' ', '_' );
        sOptions += sOut + "', '";
    } // for
    sOptions.pop_back( );
    sOptions.pop_back( );
    sOptions.pop_back( );
    printOption(
        "Presetting",
        'p',
        "name",
        rManager.xParametersSets.begin( )->first,
        "Optimize aligner parameters for a selected sequencing technique. Available presettings are: " + sOptions + ".",
        sIndentDesc );
    // general options
    std::cout << "General options: (these options are not affected by presettings)" << std::endl;
    printOption( "Index",
                 'x',
                 "file_name",
                 "",
                 "Filename of FMD-index. (A FMD-index can be generated via the --Create_Index option.) This option "
                 "must be set.",
                 sIndentDesc );
    printOption( "In",
                 'i',
                 "file_name",
                 "",
                 "Filenames of Fasta/Fastq files containing reads. gz-compressed files are automatically decompressed. "
                 "Multiple files can be specified by a comma separated list. One file name must be provided at least.",
                 sIndentDesc );
    printOption( "Mate_In",
                 'm',
                 "file_name",
                 "",
                 "Filenames of the mates in the case of paired reads. If this option is set, the aligner switches to "
                 "paired mode automatically. The number of reads given as mates must match the accumulated "
                 "number of reads provided via the 'in'-option.",
                 sIndentDesc );
    printOption( "Create_Index",
                 'X',
                 "fasta_file_name,output_folder,index_name",
                 "",
                 "Generate a FMD-index for a Fasta file. 'fasta_file_name' has to be the file-path of the Fasta file "
                 "holding the genome used for index creation. 'output_folder' is the folder-path of the location used "
                 "for index storage. 'index_name' is the name used for identifying the new FMD-Index. In the context "
                 "of alignments, the genome-name is used for FMD-index selection.",
                 sIndentDesc );
    for( auto xPair : rManager.pGlobalParameterSet->xpParametersByCategory )
    {
        for( auto pParameter : xPair.second )
            printOption( pParameter->sName,
                         pParameter->cShort,
                         pParameter->type_name( ),
                         pParameter->asText( ),
                         pParameter->sDescription,
                         sIndentDesc );
    } // for

    if( bFull )
    {
        // other options
        for( auto xPair : rManager.getSelected( )->xpParametersByCategory )
        {
            std::cout << xPair.first.second << " options:" << std::endl;
            for( auto pParameter : xPair.second )
                printOption( pParameter->sName,
                             pParameter->cShort,
                             pParameter->type_name( ),
                             pParameter->asText( ),
                             pParameter->sDescription,
                             sIndentDesc );
        } // for
    } // if

    std::cout << "Version " << MA_VERSION << "\nBy Markus Schmidt & Arne Kutzner" << std::endl;
    std::cout << "Compiled with following switches:";
    if( bLibMaWithPython )
        std::cout << " WITH_PYTHON";
#ifdef WITH_POSTGRES
    std::cout << " WITH_POSTGRES";
#endif
#ifdef WITH_ZLIB
    std::cout << " WITH_ZLIB";
#endif
#if DEBUG_LEVEL > 0
    std::cout << " DEBUG_MODE";
#endif
    std::cout << "\nFor more information visit: https://github.com/ITBE-Lab/ma" << std::endl;
} // function


std::vector<fs::path> fsSplit( const std::string& sSubject, const std::string sRegex )
{
    std::vector<fs::path> vVector;
    for( std::string sPath : split( sSubject, sRegex ) )
        vVector.push_back( fs::path( sPath ) );
    return vVector;
} // function

/**
 * main function
 */
int main( int argc, char* argv[] )
{
    if( MA_VERSION != sLibMaVersion )
    {
        std::cerr << "Fatal error: cmbMA verion \"" << MA_VERSION << "\" does not match libMA version \""
                  << sLibMaVersion << "\". Something went wrong during building/linking." << std::endl;
        return 1;
    } // if
    ExecutionContext xExecutionContext;
    // change the way output works to a simple -o for the command line aligner.
    // Also disable Use Max Hardware concurrency parameter and set -t to max hardware_concurrency by default.
    xExecutionContext.xParameterSetManager.pGlobalParameterSet->xSAMOutputTypeChoice->uiSelection = 2;
    xExecutionContext.xParameterSetManager.pGlobalParameterSet->pbUseMaxHardareConcurrency->set( false );
    xExecutionContext.xParameterSetManager.pGlobalParameterSet->piNumberOfThreads->set(
        std::thread::hardware_concurrency( ) );
    // remove not with respect to pbUseMaxHardareConcurrency in description...
    xExecutionContext.xParameterSetManager.pGlobalParameterSet->piNumberOfThreads->sDescription =
        "Number of threads used in the context of alignments.";
    xExecutionContext.xParameterSetManager.pGlobalParameterSet->unregisterParameter(
        xExecutionContext.xParameterSetManager.pGlobalParameterSet->xSAMOutputTypeChoice.pContent );
    xExecutionContext.xParameterSetManager.pGlobalParameterSet->unregisterParameter(
        xExecutionContext.xParameterSetManager.pGlobalParameterSet->xSAMOutputPath.pContent );
    xExecutionContext.xParameterSetManager.pGlobalParameterSet->unregisterParameter(
        xExecutionContext.xParameterSetManager.pGlobalParameterSet->pbUseMaxHardareConcurrency.pContent );

    // set the mode...
    for( int iI = 2; iI < argc; iI += 2 )
    {
        std::string sOptionName = argv[ iI - 1 ];
        std::string sOptionValue = argv[ iI ];
        if( sOptionName == "-p" || ParameterSetBase::uniqueParameterName(sOptionName) == "--presetting" )
            xExecutionContext.xParameterSetManager.setSelected( sOptionValue );
    } // for

    if( argc <= 1 )
    {
        generateHelpMessage( xExecutionContext.xParameterSetManager, false );
        return 0;
    } // if

    try
    {
        for( int iI = 1; iI < argc; iI++ )
        {
            std::string sOptionName = argv[ iI ];

            // we did this already
            if( sOptionName == "-p" || ParameterSetBase::uniqueParameterName(sOptionName) == "--presetting" )
            {
                iI++; // also ignore the following argument
                continue;
            } // if

            if( sOptionName == "-x" || ParameterSetBase::uniqueParameterName(sOptionName) == "--index" )
            {
                std::string sOptionValue = argv[ iI + 1 ];
                const std::string s = xExecutionContext.xGenomeManager.loadGenome( sOptionValue );
                if( !s.empty( ) )
                    throw std::runtime_error( s );
                iI++; // also ignore the following argument
                continue;
            } // if

            if( sOptionName == "-i" || ParameterSetBase::uniqueParameterName(sOptionName) == "--in" )
            {
                std::string sOptionValue = argv[ iI + 1 ];
                xExecutionContext.xReadsManager.vsPrimaryQueryFullFileName = fsSplit( sOptionValue, "," );
                iI++; // also ignore the following argument
                continue;
            } // if

            if( sOptionName == "-m" || ParameterSetBase::uniqueParameterName(sOptionName) == "--matein" )
            {
                std::string sOptionValue = argv[ iI + 1 ];
                xExecutionContext.xReadsManager.vsMateQueryFullFileName = fsSplit( sOptionValue, "," );
                xExecutionContext.xParameterSetManager.getSelected( )->xUsePairedReads->set( true );
                iI++; // also ignore the following argument
                continue;
            } // if

            if( sOptionName == "-X" || ParameterSetBase::uniqueParameterName(sOptionName) == "--createindex" )
            {
                std::string sOptionValue = argv[ iI + 1 ];
                auto vsStrings = split( sOptionValue, "," );
                if( vsStrings.size( ) != 3 )
                    throw std::runtime_error( "--Index needs exactly three parameters" );
                xExecutionContext.xGenomeManager.makeIndexAndPackForGenome(
                    fs::path( vsStrings[ 1 ] ), //
                    fs::path( vsStrings[ 0 ] ), //
                    vsStrings[ 2 ], //
                    []( const std::string s ) { std::cout << s << std::endl; } // lambda
                );
                return 0;
            } // if

            if( iI + 1 < argc && ( argv[ iI + 1 ][ 0 ] != '-' || is_number( std::string( argv[ iI + 1 ] ) ) ) )
            {
                std::string sOptionValue = argv[ iI + 1 ];
                iI++; // have key value pair so next element is certainly no key

                if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] != '-' && sOptionName.size( ) == 2 )
                    xExecutionContext.xParameterSetManager.byShort( sOptionName[ 1 ] )->setByText( sOptionValue );

                else if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] == '-' && sOptionName.size( ) > 2 )
                    xExecutionContext.xParameterSetManager.byName( sOptionName.substr( 2, sOptionName.size( ) - 2 ) )
                        ->setByText( sOptionValue );

                else
                    throw std::runtime_error(
                        std::string( "unknown option type: " )
                            .append( sOptionName )
                            .append( ". Did you forget to add the '-' or '--' at the beginning?" ) );
            } // if
            else // boolean flag option
            {
                if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] != '-' && sOptionName.size( ) == 2 )
                {
                    auto pX = std::dynamic_pointer_cast<AlignerParameter<bool>>(
                        xExecutionContext.xParameterSetManager.byShort( sOptionName[ 1 ] ) );

                    if( pX == nullptr )
                        throw std::runtime_error( "Parameters need to be provided as key value pairs" );
                    pX->set( true );
                } // if
                else if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] == '-' && sOptionName.size( ) > 2 )
                {
                    auto pX = std::dynamic_pointer_cast<AlignerParameter<bool>>(
                        xExecutionContext.xParameterSetManager.byName(
                            sOptionName.substr( 2, sOptionName.size( ) - 2 ) ) );
                    if( pX == nullptr )
                        throw std::runtime_error( "Parameters need to be provided as key value pairs" );
                    pX->set( true );
                } // else if
                else
                    throw std::runtime_error(
                        std::string( "unknown option type: " )
                            .append( sOptionName )
                            .append( ". Did you forget to add the '-' or '--' at the beginning?" ) );
            } // else
        } // for
        if( xExecutionContext.xParameterSetManager.pGlobalParameterSet->pbPrintHelpMessage->get( ) )
        {
            generateHelpMessage( xExecutionContext.xParameterSetManager );
            return 0;
        } // if

        std::pair<int, double> xPreviousProgress = std::make_pair( -1, 0 );
        std::cout << "starting alignment." << std::endl;
        xExecutionContext.doAlign( [&] //
                                   ( double dProgress, int iCurrFile, int iFilesTotal ) //
                                   {
                                       dProgress = (int)( dProgress * 10 );
                                       dProgress /= 10;
                                       std::pair<int, double> xProgress = std::make_pair( iCurrFile, dProgress );
                                       if( xProgress > xPreviousProgress )
                                       {
                                           std::cerr << "\rFile " << xProgress.first + 1 << " of " << iFilesTotal
                                                     << ": " << xProgress.second << "% aligned.             "
                                                     << std::flush;
                                           xPreviousProgress = xProgress;
                                       } // if
                                       return true; // always continue the alignment
                                   } // lambda
        );
        std::cerr << "\rdone.                         " << std::endl;
    } // try
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
    return 0;
} // main function