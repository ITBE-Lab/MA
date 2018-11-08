/**
 * @file fileReader.cpp
 * @author Markus Schmidt
 */
#include "module/fileReader.h"
#include "util/pybind11.h"

using namespace libMA;

size_t len( std::string& sLine )
{
    size_t uiLineSize = sLine.length( );
    while( uiLineSize > 0 && sLine[ uiLineSize - 1 ] != 'A' && sLine[ uiLineSize - 1 ] != 'C' &&
           sLine[ uiLineSize - 1 ] != 'T' && sLine[ uiLineSize - 1 ] != 'G' && sLine[ uiLineSize - 1 ] != 'N' &&
           sLine[ uiLineSize - 1 ] != 'a' && sLine[ uiLineSize - 1 ] != 'c' && sLine[ uiLineSize - 1 ] != 't' &&
           sLine[ uiLineSize - 1 ] != 'g' && sLine[ uiLineSize - 1 ] != 'n' )
        uiLineSize--;
    return uiLineSize;
} // function

std::shared_ptr<NucSeq> FileReader::execute( )
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
            throw AnnotatedException( "Invalid line in fasta" );
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
                        std::cerr << "WARNING: " << sLine << " contains Ns! line: " << uiNumLinesRead
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
            size_t uiLineSize = len( sLine );
#if WITH_QUALITY == 1
            // uiLineSize uint8_t's with value 127
            std::vector<uint8_t> xQuality( uiLineSize, 126 );
#endif
            pRet->vAppend( (const uint8_t*)sLine.c_str( ),
#if WITH_QUALITY == 1
                           xQuality.data( ),
#endif
                           uiLineSize );
        } // while
        pRet->vTranslateToNumericFormUsingTable( pRet->xNucleotideTranslationTable, 0 );

        // run self tests for the nucSeq
        DEBUG_2( std::cout << pRet->fastaq( ) << std::endl; ) // DEBUG_2
        DEBUG( pRet->check( ); )
        pFile->peek( ); // peek is necessary since eof() depends on last stream operation
        if( pFile->eof( ) )
            this->setFinished( );
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
            pRet->vAppend( (const uint8_t*)sLine.c_str( ), xQuality.data( ), uiLineSize );
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
        pFile->peek( ); // peek is necessary since eof() depends on last stream operation
        if( pFile->eof( ) )
            this->setFinished( );
        return pRet;
    } // if
#else
    // FASTAQ format
    if( pFile->good( ) && !pFile->eof( ) && pFile->peek( ) == '@' )
    {
        std::string sLine = "";
        safeGetline( sLine );
        if( sLine.size( ) == 0 )
            throw AnnotatedException( "Invalid line in fastq" );
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
                        std::cerr << "WARNING: " << sLine << " contains Ns! line: " << uiNumLinesRead
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
        safeGetline( sLine );
        if( sLine.size( ) != 1 || sLine[ 0 ] != '+' )
            throw AnnotatedException( "Invalid line in fastq" );
        while( pFile->good( ) && !pFile->eof( ) && uiNumChars > 0 )
        {
            safeGetline( sLine );
            uiNumChars -= sLine.size( );
        } // while
        // if(pFile->good() && !pFile->eof() && pFile->peek() != '@')
        //    throw AnnotatedException("Invalid line in fastq");
        pFile->peek( ); // peek is necessary since eof() depends on last stream operation
        if( pFile->eof( ) )
            this->setFinished( );
        return pRet;
    } // if
#endif
    // if we reach this point we have read all content of the file
    throw AnnotatedException( "Tried to read query past EoF" );
} // function

std::shared_ptr<TP_PAIRED_READS> PairedFileReader::execute( )
{
    auto pRet = std::make_shared<TP_PAIRED_READS>( );
    pRet->push_back( xF1.execute( ) );
    pRet->push_back( xF2.execute( ) );
    // forward the finished flags...
    if( xF1.isFinished( ) || xF2.isFinished( ) )
        this->setFinished( );
    if( ( *pRet )[ 0 ] == nullptr )
        return nullptr;
    if( ( *pRet )[ 1 ] == nullptr )
        return nullptr;
    return pRet;
} // function

#ifdef WITH_PYTHON

#ifdef BOOST_PYTHON
void exportFileReader( )
{
    // export the FileReader class
    exportModule<FileReader, std::string>( "FileReader" );

    boost::python::
        class_<TP_PAIRED_READS, boost::noncopyable, boost::python::bases<Container>, std::shared_ptr<TP_PAIRED_READS>>(
            "QueryVector" )
            /*
             * true = noproxy this means that the content of
             * the vector is already exposed by boost python.
             * if this is kept as false, Container would be
             * exposed a second time. the two Containers would
             * be different and not inter castable.
             */
            .def( boost::python::vector_indexing_suite<TP_PAIRED_READS, true>( ) );
    // export the PairedFileReader class
    exportModule<PairedFileReader, std::string, std::string>( "PairedFileReader" );
} // function
#else
void exportFileReader( py::module& rxPyModuleId )
{
    // export the FileReader class
    exportModule<FileReader, std::string>( rxPyModuleId, "FileReader" );

    py::bind_vector_ext<TP_PAIRED_READS, Container, std::shared_ptr<TP_PAIRED_READS>>( rxPyModuleId, "QueryVector" );

    // export the PairedFileReader class
    exportModule<PairedFileReader, std::string, std::string>( rxPyModuleId, "PairedFileReader" );
} // function
#endif
#endif