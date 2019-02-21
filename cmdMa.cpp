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
#include "util/debug.h"
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
    std::cout << sIndentOptions;
    if( cShort != AlignerParameterBase::NO_SHORT_DEFINED )
        std::cout << "-" << cShort << ", ";
    std::replace( sName.begin( ), sName.end( ), ' ', '_' );
    std::cout << "--" << sName << " <" << sTypeName << "> [";
    std::cout << sDefaultVal << "]";
    std::cout << std::endl << sIndentDesc;

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

void generateHelpMessage( ParameterSetManager& rManager )
{
    std::string sIndentDesc;
    for( size_t uiI = 0; uiI < sHeader.size( ) / 2; uiI++ )
        sIndentDesc += " ";
    std::cout << sHeader << std::endl;

    // presettings

    std::cout << "Available presettings:" << std::endl;
    std::string sOptions;
    for( auto& xPair : rManager.xParametersSets )
    {
        std::string sOut = xPair.first;
        std::replace( sOut.begin( ), sOut.end( ), ' ', '_' );
        sOptions += sOut + "/";
    }
    sOptions.pop_back( );
    printOption( "Presetting",
                 'p',
                 sOptions,
                 rManager.xParametersSets.begin( )->first,
                 "Operation mode for MA. The selected mode will change the parameters. However, parameters you "
                 "set manually will overwrite the presetting.",
                 sIndentDesc );

    // general options
    std::cout << "General options: (these options are not affected by presettings)" << std::endl;
    printOption( "genome", 'x', "path", "", "Path to the genome file. This option must always be set.", sIndentDesc );
    printOption(
        "in",
        'i',
        "path",
        "",
        "Path to the input file. For multiple files: provide comma seperated list. This option must always be set.",
        sIndentDesc );
    printOption( "mate_in",
                 'm',
                 "path",
                 "",
                 "Path to the input file of mates. Automatically activates paired alignment mode. For multiple files: "
                 "provide comma seperated list.",
                 sIndentDesc );
    for( auto xTup : rManager.xGlobalParameterSet.xpAllParameters )
    {
        auto pParameter = xTup.second;
        printOption( pParameter->sName,
                     pParameter->cShort,
                     pParameter->type_name( ),
                     pParameter->asText( ),
                     pParameter->sDescription,
                     sIndentDesc );
    } // for

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

    std::cout << "Version " << MA_VERSION << "\nBy Markus Schmidt & Arne Kutzner" << std::endl;
    std::cout << "For more information visit: https://github.com/ITBE-Lab/ma" << std::endl;
#if DEBUG_LEVEL > 0
    std::cout << "DEBUG MODE" << std::endl;
#endif
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
    ExecutionContext xExecutionContext;
    xExecutionContext.xParameterSetManager.xGlobalParameterSet.bSAMOutputInReadsFolder->set( false );

    // set the mode...
    for( int iI = 2; iI < argc; iI += 2 )
    {
        std::string sOptionName = argv[ iI - 1 ];
        std::string sOptionValue = argv[ iI ];
        std::replace( sOptionValue.begin( ), sOptionValue.end( ), '_', ' ' );
        if( sOptionName == "-m" || sOptionName == "--mode" )
            xExecutionContext.xParameterSetManager.setSelected( sOptionValue );
    } // for

    if( argc <= 1 )
    {
        generateHelpMessage( xExecutionContext.xParameterSetManager );
        return 0;
    } // if

    try
    {
        for( int iI = 1; iI < argc; iI++ )
        {
            std::string sOptionName = argv[ iI ];

            if( sOptionName == "-p" || sOptionName == "--Presetting" ) // we did this already
            {
                iI++; // also ignore the following argument
                continue;
            } // if

            if( sOptionName == "-x" || sOptionName == "--genome" )
            {
                std::string sOptionValue = argv[ iI + 1 ];
                xExecutionContext.xGenomeManager.loadGenome( sOptionValue + ".json" );
                iI++; // also ignore the following argument
                continue;
            } // if

            if( sOptionName == "-i" || sOptionName == "--in" )
            {
                std::string sOptionValue = argv[ iI + 1 ];
                xExecutionContext.xReadsManager.vsPrimaryQueryFullFileName = fsSplit( sOptionValue, "," );
                iI++; // also ignore the following argument
                continue;
            } // if

            if( sOptionName == "-m" || sOptionName == "--mate_in" )
            {
                std::string sOptionValue = argv[ iI + 1 ];
                xExecutionContext.xReadsManager.vsMateQueryFullFileName = fsSplit( sOptionValue, "," );
                xExecutionContext.xParameterSetManager.getSelected( )->xUsePairedReads->set( true );
                iI++; // also ignore the following argument
                continue;
            } // if

            std::replace( sOptionName.begin( ), sOptionName.end( ), '_', ' ' );

            if( iI + 1 < argc && argv[ iI + 1 ][ 0 ] != '-' )
            {
                std::string sOptionValue = argv[ iI + 1 ];
                iI++; // have key value pair so next element is certainly no key

                if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] != '-' && sOptionName.size( ) == 2 )
                    xExecutionContext.xParameterSetManager.byShort( sOptionName[ 1 ] )->setByText( sOptionValue );

                else if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] == '-' && sOptionName.size( ) > 2 )
                    xExecutionContext.xParameterSetManager.byName( sOptionName.substr( 2, sOptionName.size( ) - 2 ) )
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
                        xExecutionContext.xParameterSetManager.byShort( sOptionName[ 1 ] ) );

                    if( pX == nullptr )
                        throw AnnotatedException( "Parameters need to be provided as key value pairs" );
                    pX->set( true );
                } // if
                else if( sOptionName[ 0 ] == '-' && sOptionName[ 1 ] == '-' && sOptionName.size( ) > 2 )
                {
                    auto pX = std::dynamic_pointer_cast<AlignerParameter<bool>>(
                        xExecutionContext.xParameterSetManager.byName(
                            sOptionName.substr( 2, sOptionName.size( ) - 2 ) ) );
                    if( pX == nullptr )
                        throw AnnotatedException( "Parameters need to be provided as key value pairs" );
                    pX->set( true );
                } // else if
                else
                    throw AnnotatedException(
                        std::string( "unknown option type: " )
                            .append( sOptionName )
                            .append( ". Did you forget to add the '-' or '--' at the beginning?" ) );
            } // else
        } // for
        if( xExecutionContext.xParameterSetManager.xGlobalParameterSet.pbPrintHelpMessage->get( ) )
            generateHelpMessage( xExecutionContext.xParameterSetManager );

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
    return 0;
} // main function