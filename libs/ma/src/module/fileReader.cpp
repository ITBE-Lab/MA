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

std::shared_ptr<NucSeq> FileReader::execute( )
{
    // std::lock_guard<std::mutex> xGuard(*pSynchronizeReading);
    std::shared_ptr<NucSeq> pRet( new NucSeq( ) );
    DEBUG( pRet->uiFromLine = pFile->uiNumLinesRead; )
    // FASTA format
    if( !pFile->eof( ) && pFile->peek( ) == '>' )
    {
        std::string sLine = "";
        pFile->safeGetLine( sLine );
        if( sLine.size( ) == 0 )
            throw AnnotatedException( "Invalid line in fasta" );
        ;
        // make sure that the name contains no spaces
        // in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr( 1, sLine.find( ' ' ) );

        while( !pFile->eof( ) && pFile->peek( ) != '>' && pFile->peek( ) != ' ' )
        {
            sLine = ""; // in the case that we hit an empty line getline does nothing...
            pFile->safeGetLine( sLine );
            if( sLine.size( ) == 0 )
                continue;
            DEBUG( for( auto character
                        : sLine ) {
                bool bOkay = false;
                if( character == 'N' || character == 'n' )
                {
                    if( uiNumLinesWithNs == 0 )
                        std::cerr << "WARNING: " << sLine << " contains Ns! line: " << pFile->uiNumLinesRead
                                  << " (this is only printed once)" << std::endl;
                    uiNumLinesWithNs++;
                    continue;
                }
                for( char c : std::vector<char>{'A', 'C', 'T', 'G', 'a', 'c', 't', 'g'} )
                    if( c == character )
                        bOkay = true;
                if( !bOkay )
                {
                    std::cerr << "Invalid symbol in fasta: " << sLine << std::endl;
                    throw AnnotatedException( "Invalid symbol in fasta" );
                } // if
            } // for
                   ) // DEBUG
            pRet->vAppend( (const uint8_t*)sLine.c_str( ), len( sLine ) );
        } // while
        pRet->vTranslateToNumericFormUsingTable( pRet->xNucleotideTranslationTable, 0 );

        // run self tests for the nucSeq
        DEBUG_2( std::cout << pRet->fastaq( ) << std::endl; ) // DEBUG_2
        DEBUG( pRet->check( ); )
        pFile->peek( ); // peek is necessary since eof() depends on last stream operation
        if( pFile->eof( ) )
            this->setFinished( );
        if( pRet->length( ) == 0 )
            throw std::runtime_error( "found empty read: " + pRet->sName );
        return pRet;
    } // if
#if WITH_QUALITY == 1
    // FASTAQ format
    if( !pFile->eof( ) && pFile->peek( ) == '@' )
    {
        pRet->addQuality( );
        std::string sLine;
        pFile->safeGetLine( sLine );
        if( sLine.size( ) == 0 )
            throw AnnotatedException( "Invalid line in fasta" );
        // make sure that the name contains no spaces
        // in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr( 1, sLine.find( ' ' ) );
        while( !pFile->eof( ) && pFile->peek( ) != '+' && pFile->peek( ) != ' ' )
        {
            sLine = "";
            pFile->safeGetLine( sLine );
            if( sLine.size( ) == 0 )
                continue;
            size_t uiLineSize = len( sLine );
            pRet->vAppend( (const uint8_t*)sLine.c_str( ), uiLineSize );
        } // while
        pRet->vTranslateToNumericFormUsingTable( pRet->xNucleotideTranslationTable, 0 );
        pFile->safeGetLine( sLine );
        // quality
        if( sLine[ 0 ] == '+' )
        {
            size_t uiPos = 0;
            while( !pFile->eof( ) && ( pFile->peek( ) != '@' || uiPos == 0 ) )
            {
                pFile->safeGetLine( sLine );
                if( sLine.size( ) == 0 )
                    continue;
                size_t uiLineSize = lenq( sLine );
                for( size_t i = 0; i < uiLineSize; i++ )
                    pRet->quality( i + uiPos ) = (uint8_t)sLine[ i ];
                uiPos += uiLineSize;
            } // while
        }
        pFile->peek( ); // peek is necessary since eof() depends on last stream operation
        if( pFile->eof( ) )
            this->setFinished( );
        if( pRet->length( ) == 0 )
            throw std::runtime_error( "found empty read: " + pRet->sName );
        return pRet;
    } // if
#else
    // FASTAQ format
    if( !pFile->eof( ) && pFile->peek( ) == '@' )
    {
        std::string sLine = "";
        pFile->safeGetLine( sLine );
        if( sLine.size( ) == 0 )
            throw AnnotatedException( "Invalid line in fastq" );
        // make sure that the name contains no spaces
        // in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr( 1, sLine.find( ' ' ) );
        size_t uiNumChars = 0;
        while( !pFile->eof( ) && pFile->peek( ) != '+' && pFile->peek( ) != ' ' )
        {
            sLine = "";
            pFile->safeGetLine( sLine );
            if( sLine.size( ) == 0 )
                continue;
            DEBUG( for( auto character
                        : sLine ) {
                bool bOkay = false;
                if( character == 'N' || character == 'n' )
                {
                    if( uiNumLinesWithNs == 0 )
                        std::cerr << "WARNING: " << sLine << " contains Ns! line: " << pFile->uiNumLinesRead
                                  << " (this is only printed once)" << std::endl;
                    uiNumLinesWithNs++;
                    continue;
                }
                for( char c : std::vector<char>{'A', 'C', 'T', 'G', 'a', 'c', 't', 'g'} )
                    if( c == character )
                        bOkay = true;
                if( !bOkay )
                {
                    std::cerr << "Invalid symbol in fasta: " << sLine << std::endl;
                    throw AnnotatedException( "Invalid symbol in fastq" );
                } // if
            } // for
                   ) // DEBUG
            size_t uiLineSize = len( sLine );
            uiNumChars += uiLineSize;
            pRet->vAppend( (const uint8_t*)sLine.c_str( ), uiLineSize );
        } // while
        pRet->vTranslateToNumericFormUsingTable( pRet->xNucleotideTranslationTable, 0 );
        // quality
        pFile->safeGetLine( sLine );
        if( sLine.size( ) != 1 || sLine[ 0 ] != '+' )
            throw AnnotatedException( "Invalid line in fastq" );
        while( !pFile->eof( ) && uiNumChars > 0 )
        {
            pFile->safeGetLine( sLine );
            uiNumChars -= sLine.size( );
        } // while
        // if(pFile->good() && !pFile->eof() && pFile->peek() != '@')
        //    throw AnnotatedException("Invalid line in fastq");
        pFile->peek( ); // peek is necessary since eof() depends on last stream operation
        if( pFile->eof( ) )
            this->setFinished( );
        if( pRet->length( ) == 0 )
            throw std::runtime_error( "found empty read: " + pRet->sName );
        return pRet;
    } // if
#endif
    std::string sLine = "EoF";
    if( !pFile->eof( ) )
        pFile->safeGetLine( sLine );
    else
        this->setFinished( );
    // if we reach this point we have read all content of the file
    throw AnnotatedException(
        "Error while reading file.\nIs your input really in FASTA/Q format?\nError occurred in file: " +
        xFileName.string( ) );
} // function

std::shared_ptr<TP_PAIRED_READS> PairedFileReader::execute( )
{
    auto pRet = std::make_shared<TP_PAIRED_READS>( );
    pRet->push_back( pF1->execute( ) );
    pRet->push_back( pF2->execute( ) );
    if( pF1->getCurrFileIndex( ) != pF2->getCurrFileIndex( ) )
        throw std::runtime_error(
            "Cannot perfrom paired alignment on files with different amounts of reads. FileReader Status: " +
            this->status( ) );
    // forward the finished flags...
    if( pF1->isFinished( ) || pF2->isFinished( ) )
    {
        /*
         * Print a warning if the fasta files have a different number of queries.
         */
        if( !pF1->isFinished( ) || !pF2->isFinished( ) )
            throw std::runtime_error(
                "You cannot perform paired alignment with a different amount of primary queries and mate queries." );
        this->setFinished( );
    }
    if( ( *pRet )[ 0 ] == nullptr )
        return nullptr;
    if( ( *pRet )[ 1 ] == nullptr )
        return nullptr;

    if( bRevCompMate )
    {
        pRet->back( )->vReverse( );
        pRet->back( )->vSwitchAllBasePairsToComplement( );
    } // if
    return pRet;
} // function

void PairedFileReader::checkPaired( )
{
    while( !pF1->isFinished( ) && !pF2->isFinished( ) )
    {
        auto pQ1 = pF1->execute( );
        auto pQ2 = pF2->execute( );
        // if( pQ1->sName != pQ2->sName )
        //     throw std::runtime_error( "paired queries with different names: " + pQ1->sName + " != " + pQ2->sName +
        //                               " FileReader Status:" + this->status( ) );
        if( pF1->getCurrFileIndex( ) != pF2->getCurrFileIndex( ) )
            throw std::runtime_error(
                "Cannot perfrom paired alignment on files with different amounts of reads. FileReader Status: " +
                this->status( ) );
    } // while
    if( !pF1->isFinished( ) || !pF2->isFinished( ) )
        throw std::runtime_error(
            "Cannot perfrom paired alignment on files with different amounts of reads. FileReader Status: " +
            this->status( ) );
    pF1->reset( );
    pF2->reset( );
} // method

#ifdef WITH_PYTHON

void exportFileReader( py::module& rxPyModuleId )
{
    py::class_<fs::path>( rxPyModuleId, "path" ).def( py::init<std::string>( ) );
    //.def( "__str__", &fs::path::string );

    py::bind_vector<std::vector<fs::path>>( rxPyModuleId, "filePathVector", "docstr" );
    // export the FileReader class
    exportModule<FileReader, fs::path>( rxPyModuleId, "FileReader" );
    exportModule<FileListReader, std::vector<fs::path>>(
        rxPyModuleId, "FileListReader", []( auto&& x ) { x.def( "read_all", &FileListReader::read_all ); } );

    py::bind_vector_ext<TP_PAIRED_READS, Container, std::shared_ptr<TP_PAIRED_READS>>(
        rxPyModuleId, "ContainerVectorNucSeq", "docstr" );

    // export the PairedFileReader class
    exportModule<PairedFileReader, std::vector<fs::path>, std::vector<fs::path>>( rxPyModuleId, "PairedFileReader" );
} // function
#endif