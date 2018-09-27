/**
 * @file fileReader.cpp
 * @author Markus Schmidt
 */
#include "module/fileReader.h"

using namespace libMA;

ContainerVector FileReader::getInputType( ) const
{
    return ContainerVector{std::shared_ptr<Container>( new Nil( ) )};
} // function

std::shared_ptr<Container> FileReader::getOutputType( ) const
{
    return std::shared_ptr<Container>( new NucSeq( ) );
} // function

size_t len( std::string &sLine )
{
    size_t uiLineSize = sLine.length( );
    while( uiLineSize > 0 && sLine[ uiLineSize - 1 ] != 'A' && sLine[ uiLineSize - 1 ] != 'C' &&
           sLine[ uiLineSize - 1 ] != 'T' && sLine[ uiLineSize - 1 ] != 'G' &&
           sLine[ uiLineSize - 1 ] != 'N' && sLine[ uiLineSize - 1 ] != 'a' &&
           sLine[ uiLineSize - 1 ] != 'c' && sLine[ uiLineSize - 1 ] != 't' &&
           sLine[ uiLineSize - 1 ] != 'g' && sLine[ uiLineSize - 1 ] != 'n' )
        uiLineSize--;
    return uiLineSize;
}

#if USE_BUFFERED_ASYNC_READER == 1
std::shared_ptr<Container> FileReader::execute( std::shared_ptr<ContainerVector> vpInput )
{
    /*
     * Has next and next require synchronized access.
     * This is done by the module synchronization.
     */
    if( pFile->hasNext( ) )
        return pFile->next( );

    // if we reach this point we have read all content of the file
    return Nil::pEoFContainer;
} // function
#else

std::shared_ptr<Container> FileReader::execute( std::shared_ptr<ContainerVector> vpInput )
{
    // std::lock_guard<std::mutex> xGuard(*pSynchronizeReading);
    std::shared_ptr<NucSeq> pRet( new NucSeq( ) );
    DEBUG( pRet->uiFromLine = uiNumLinesRead; )
    // FASTA format
    if( pFile->good( ) && !pFile->eof( ) && pFile->peek( ) == '>' )
    {
        std::string sLine = "";
        safeGetline( sLine );
        if( sLine.size( ) == 0 )
            throw AlignerException( "Invalid line in fasta" );
        ;
        // make sure that the name contains no spaces
        // in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr( 1, sLine.find( ' ' ) );

        while( pFile->good( ) && !pFile->eof( ) && pFile->peek( ) != '>' && pFile->peek( ) != ' ' )
        {
            sLine = ""; // in the case that we hit an empty line getline does nothing...
            safeGetline( sLine );
            if( sLine.size( ) == 0 )
                continue;
            DEBUG( for( auto character
                        : sLine ) {
                bool bOkay = false;
                if( character == 'N' || character == 'n' )
                {
                    if( uiNumLinesWithNs == 0 )
                        std::cerr << "WARNING: " << sLine
                                  << " contains Ns! line: " << uiNumLinesRead
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
                    throw AlignerException( "Invalid symbol in fasta" );
                } // if
            } // for
                   ) // DEBUG
            size_t uiLineSize = len( sLine );
#if WITH_QUALITY == 1
            // uiLineSize uint8_t's with value 127
            std::vector<uint8_t> xQuality( uiLineSize, 126 );
#endif
            pRet->vAppend( (const uint8_t *)sLine.c_str( ),
#if WITH_QUALITY == 1
                           xQuality.data( ),
#endif
                           uiLineSize );
        } // while
        pRet->vTranslateToNumericFormUsingTable( pRet->xNucleotideTranslationTable, 0 );

        // run self tests for the nucSeq
        DEBUG_2( std::cout << pRet->fastaq( ) << std::endl; ) // DEBUG_2
        DEBUG( pRet->check( ); )
        return pRet;
    } // if
#if WITH_QUALITY == 1
    // FASTAQ format
    if( pFile->good( ) && !pFile->eof( ) && pFile->peek( ) == '@' )
    {
        std::string sLine;
        safeGetline( sLine );
        if( sLine.size( ) == 0 )
            throw AlignerException( "Invalid line in fasta" );
        ;
        // make sure that the name contains no spaces
        // in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr( 1, sLine.find( ' ' ) );
        while( pFile->good( ) && !pFile->eof( ) && pFile->peek( ) != '+' && pFile->peek( ) != ' ' )
        {
            sLine = "";
            safeGetline( sLine );
            if( sLine.size( ) == 0 )
                continue;
            size_t uiLineSize = len( sLine );
            std::vector<uint8_t> xQuality( uiLineSize, 126 ); // uiLineSize uint8_t's with value 127
            pRet->vAppend( (const uint8_t *)sLine.c_str( ), xQuality.data( ), uiLineSize );
        } // while
        pRet->vTranslateToNumericFormUsingTable( pRet->xNucleotideTranslationTable, 0 );
        // quality
        unsigned int uiPos = 0;
        while( pFile->good( ) && !pFile->eof( ) && pFile->peek( ) != '@' )
        {
            safeGetline( sLine );
            size_t uiLineSize = len( sLine );
            for( size_t i = 0; i < uiLineSize; i++ )
                pRet->quality( i + uiPos ) = (uint8_t)sLine[ i ];
            uiPos += uiLineSize;
        } // while
        return pRet;
    } // if
#else
    // FASTAQ format
    if( pFile->good( ) && !pFile->eof( ) && pFile->peek( ) == '@' )
    {
        std::string sLine = "";
        safeGetline( sLine );
        if( sLine.size( ) == 0 )
            throw AlignerException( "Invalid line in fastq" );
        // make sure that the name contains no spaces
        // in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr( 1, sLine.find( ' ' ) );
        size_t uiNumChars = 0;
        while( pFile->good( ) && !pFile->eof( ) && pFile->peek( ) != '+' && pFile->peek( ) != ' ' )
        {
            sLine = "";
            safeGetline( sLine );
            if( sLine.size( ) == 0 )
                continue;
            DEBUG( for( auto character
                        : sLine ) {
                bool bOkay = false;
                if( character == 'N' || character == 'n' )
                {
                    if( uiNumLinesWithNs == 0 )
                        std::cerr << "WARNING: " << sLine
                                  << " contains Ns! line: " << uiNumLinesRead
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
                    throw AlignerException( "Invalid symbol in fastq" );
                } // if
            } // for
                   ) // DEBUG
            size_t uiLineSize = len( sLine );
            uiNumChars += uiLineSize;
            pRet->vAppend( (const uint8_t *)sLine.c_str( ), uiLineSize );
        } // while
        pRet->vTranslateToNumericFormUsingTable( pRet->xNucleotideTranslationTable, 0 );
        // quality
        safeGetline( sLine );
        if( sLine.size( ) != 1 || sLine[ 0 ] != '+' )
            throw AlignerException( "Invalid line in fastq" );
        while( pFile->good( ) && !pFile->eof( ) && uiNumChars > 0 )
        {
            safeGetline( sLine );
            uiNumChars -= sLine.size( );
        } // while
        // if(pFile->good() && !pFile->eof() && pFile->peek() != '@')
        //    throw AlignerException("Invalid line in fastq");
        return pRet;
    } // if
#endif
    // if we reach this point we have read all content of the file
    return Nil::pEoFContainer;
} // function
#endif

ContainerVector PairedFileReader::getInputType( ) const
{
    return ContainerVector{std::shared_ptr<Container>( new Nil( ) )};
} // function

std::shared_ptr<Container> PairedFileReader::getOutputType( ) const
{
    return std::make_shared<ContainerVector>( std::make_shared<NucSeq>( ) );
} // function


std::shared_ptr<Container> PairedFileReader::execute( std::shared_ptr<ContainerVector> vpInput )
{
    auto pRet = std::make_shared<ContainerVector>( );
    pRet->push_back( xF1.execute( vpInput ) );
    pRet->push_back( xF2.execute( vpInput ) );
    if( pRet->front( ) == Nil::pEoFContainer )
        return Nil::pEoFContainer;
    if( pRet->back( ) == Nil::pEoFContainer )
        return Nil::pEoFContainer;
    return pRet;
} // function

#ifdef WITH_PYTHON
void exportFileReader( )
{
    // export the FileReader class
    boost::python::class_<FileReader, boost::python::bases<Module>, std::shared_ptr<FileReader>>(
        "FileReader", boost::python::init<std::string>( ) )
#if USE_BUFFERED_ASYNC_READER == 1
        DEBUG(.def( "testBufReader", &FileReader::testBufReader ).staticmethod( "testBufReader" ) )
#endif
            ;
    boost::python::implicitly_convertible<std::shared_ptr<FileReader>, std::shared_ptr<Module>>( );


    boost::python::
        class_<PairedFileReader, boost::python::bases<Module>, std::shared_ptr<PairedFileReader>>(
            "PairedFileReader", boost::python::init<std::string, std::string>( ) );
    boost::python::implicitly_convertible<std::shared_ptr<PairedFileReader>,
                                          std::shared_ptr<Module>>( );

} // function
#endif