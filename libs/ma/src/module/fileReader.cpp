/**
 * @file fileReader.cpp
 * @author Markus Schmidt
 */
#include "module/fileReader.h"
#include "util/pybind11.h"
#include <cctype>

using namespace libMA;

bool validNuc( char c )
{
    for( char c2 : {'A', 'C', 'G', 'T', 'N', 'U', 'R', 'Y', 'K', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V'} )
        if( c2 == toupper( c ) )
            return true;
    return false;
} // method

size_t len( std::string& sLine )
{
    size_t uiLineSize = sLine.length( );
    while( uiLineSize > 0 && !validNuc( sLine[ uiLineSize - 1 ] ) )
        uiLineSize--;
    return uiLineSize;
} // function

size_t lenq( std::string& sLine )
{
    size_t uiLineSize = sLine.length( );
    while( uiLineSize > 0 && ( sLine[ uiLineSize - 1 ] == '\n' || sLine[ uiLineSize - 1 ] == '\r' ) )
        uiLineSize--;
    return uiLineSize;
} // function

std::shared_ptr<NucSeq> FileReader::execute( std::shared_ptr<FileStream> pStream )
{
    std::unique_lock<std::mutex> xLock( pStream->xMutex );
    pStream->peek( );
    if( pStream->eof( ) ) // eof case
        return nullptr;
    auto pRet = std::make_shared<NucSeq>( );
    DEBUG( pRet->uiFromLine = pStream->uiNumLinesRead; )
    // FASTA format
    if( pStream->peek( ) == '>' )
    {
        std::string sLine = "";
        pStream->safeGetLine( sLine );
        if( sLine.size( ) == 0 )
            throw std::runtime_error( "Invalid line in fasta" );
        // make sure that the name contains no spaces
        // in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr( 1, sLine.find( ' ' ) );

        while( !pStream->eof( ) && pStream->peek( ) != '>' && pStream->peek( ) != ' ' )
        {
            sLine = ""; // in the case that we hit an empty line getline does nothing...
            pStream->safeGetLine( sLine );
            if( sLine.size( ) == 0 )
                continue;
            DEBUG( for( auto character
                        : sLine ) {
                bool bOkay = false;
                if( character == 'N' || character == 'n' )
                {
                    if( pStream->uiNumLinesWithNs == 0 )
                        std::cerr << "WARNING: " << sLine << " contains Ns! line: " << pStream->uiNumLinesRead
                                  << " (this is only printed once)" << std::endl;
                    pStream->uiNumLinesWithNs++;
                    continue;
                }
                for( char c : std::vector<char>{'A', 'C', 'T', 'G', 'a', 'c', 't', 'g'} )
                    if( c == character )
                        bOkay = true;
                if( !bOkay )
                {
                    std::cerr << "Invalid symbol in fasta: " << sLine << std::endl;
                    throw std::runtime_error( "Invalid symbol in fasta" );
                } // if
            } // for
                   ) // DEBUG
            pRet->vAppend( (const uint8_t*)sLine.c_str( ), len( sLine ) );
        } // while
        pRet->vTranslateToNumericFormUsingTable( pRet->xNucleotideTranslationTable, 0 );

        // run self tests for the nucSeq
        DEBUG_2( std::cout << pRet->fastaq( ) << std::endl; ) // DEBUG_2
        DEBUG( pRet->check( ); )
        if( pRet->length( ) == 0 )
            throw std::runtime_error( "found empty read: " + pRet->sName );
        return pRet;
    } // if
#if WITH_QUALITY == 1
    // FASTAQ format
    if( pStream->peek( ) == '@' )
    {
        pRet->addQuality( );
        std::string sLine;
        pStream->safeGetLine( sLine );
        if( sLine.size( ) == 0 )
            throw std::runtime_error( "Invalid line in fasta" );
        // make sure that the name contains no spaces
        // in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr( 1, sLine.find( ' ' ) );
        while( !pStream->eof( ) && pStream->peek( ) != '+' && pStream->peek( ) != ' ' )
        {
            sLine = "";
            pStream->safeGetLine( sLine );
            if( sLine.size( ) == 0 )
                continue;
            size_t uiLineSize = len( sLine );
            pRet->vAppend( (const uint8_t*)sLine.c_str( ), uiLineSize );
        } // while
        pRet->vTranslateToNumericFormUsingTable( pRet->xNucleotideTranslationTable, 0 );
        pStream->safeGetLine( sLine );
        // quality
        if( sLine[ 0 ] == '+' )
        {
            size_t uiPos = 0;
            while( !pStream->eof( ) && ( pStream->peek( ) != '@' || uiPos == 0 ) )
            {
                pStream->safeGetLine( sLine );
                if( sLine.size( ) == 0 )
                    continue;
                size_t uiLineSize = lenq( sLine );
                for( size_t i = 0; i < uiLineSize; i++ )
                    pRet->quality( i + uiPos ) = (uint8_t)sLine[ i ];
                uiPos += uiLineSize;
            } // while
        }
        if( pRet->length( ) == 0 )
            throw std::runtime_error( "found empty read: " + pRet->sName );
        return pRet;
    } // if
#else
    // FASTAQ format
    if( pStream->peek( ) == '@' )
    {
        std::string sLine = "";
        pStream->safeGetLine( sLine );
        if( sLine.size( ) == 0 )
            throw std::runtime_error( "Invalid line in fastq" );
        // make sure that the name contains no spaces
        // in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr( 1, sLine.find( ' ' ) );
        size_t uiNumChars = 0;
        while( !pStream->eof( ) && pStream->peek( ) != '+' && pStream->peek( ) != ' ' )
        {
            sLine = "";
            pStream->safeGetLine( sLine );
            if( sLine.size( ) == 0 )
                continue;
            DEBUG( for( auto character
                        : sLine ) {
                bool bOkay = false;
                if( character == 'N' || character == 'n' )
                {
                    if( pStream->uiNumLinesWithNs == 0 )
                        std::cerr << "WARNING: " << sLine << " contains Ns! line: " << pStream->uiNumLinesRead
                                  << " (this is only printed once)" << std::endl;
                    pStream->uiNumLinesWithNs++;
                    continue;
                }
                for( char c : std::vector<char>{'A', 'C', 'T', 'G', 'a', 'c', 't', 'g'} )
                    if( c == character )
                        bOkay = true;
                if( !bOkay )
                {
                    std::cerr << "Invalid symbol in fasta: " << sLine << std::endl;
                    throw std::runtime_error( "Invalid symbol in fastq" );
                } // if
            } // for
                   ) // DEBUG
            size_t uiLineSize = len( sLine );
            uiNumChars += uiLineSize;
            pRet->vAppend( (const uint8_t*)sLine.c_str( ), uiLineSize );
        } // while
        pRet->vTranslateToNumericFormUsingTable( pRet->xNucleotideTranslationTable, 0 );
        // quality
        pStream->safeGetLine( sLine );
        if( sLine.size( ) != 1 || sLine[ 0 ] != '+' )
            throw std::runtime_error( "Invalid line in fastq" );
        while( !pStream->eof( ) && uiNumChars > 0 )
        {
            pStream->safeGetLine( sLine );
            uiNumChars -= sLine.size( );
        } // while
        if( pRet->length( ) == 0 )
            throw std::runtime_error( "found empty read: " + pRet->sName );
        return pRet;
    } // if
#endif
    // if we reach this point we do not know how to interpret the next line in the file.
    throw std::runtime_error(
        "Error while reading file.\nIs your input really in FASTA/Q format?\nError occurred in file: " +
        pStream->fileName( ) );
} // function

#ifdef WITH_PYTHON

void exportFileReader( py::module& rxPyModuleId )
{
    py::class_<fs::path>( rxPyModuleId, "path" ).def( py::init<std::string>( ) );
    py::bind_vector<std::vector<fs::path>>( rxPyModuleId, "filePathVector", "docstr" );
    py::class_<FileStream, Container, std::shared_ptr<FileStream>>( rxPyModuleId, "FileStream" );
    py::class_<FileStreamFromPath, FileStream, std::shared_ptr<FileStreamFromPath>>( rxPyModuleId,
                                                                                     "FileStreamFromPath" )
        .def( py::init<fs::path>( ) )
        .def( py::init<std::string>( ) );
    py::class_<StringStream, FileStream, std::shared_ptr<StringStream>>( rxPyModuleId, "StringStream" )
        .def( py::init<std::string>( ) );

    py::bind_vector_ext<PairedReadsContainer, Container, std::shared_ptr<PairedReadsContainer>>(
        rxPyModuleId, "ContainerVectorNucSeq", "docstr" );

    exportCyclicQueue<FileStream>( rxPyModuleId, "File" );
    exportCyclicQueue<PairedFileStream>( rxPyModuleId, "PairedFile" );

    // export the FileReader class
    exportModule<FileReader>( rxPyModuleId, "FileReader" );
    // export the PairedFileReader class
    exportModule<PairedFileReader>( rxPyModuleId, "PairedFileReader" );

    rxPyModuleId.def( "combine_file_streams", &combineFileStreams );
} // function
#endif