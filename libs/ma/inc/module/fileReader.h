/**
 * @file fileReader.h
 * @brief Reads queries from a file.
 * @author Markus Schmidt
 */
#ifndef FILE_READER_H
#define FILE_READER_H

#include "container/nucSeq.h"
#include "module/module.h"
#include "util/support.h"


#ifdef WITH_ZLIB
#include "zlib.h"
#endif


/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <fstream>
#include <sstream>
/// @endcond

namespace libMA
{

class FileStream
{
  public:
    DEBUG( size_t uiNumLinesRead = 0; ) // DEBUG

    FileStream( )
    {}

    FileStream( const FileStream& ) = delete;

    virtual bool eof( ) const
    {
        throw AnnotatedException( "Unimplemented" );
    }

    virtual bool is_open( ) const
    {
        throw AnnotatedException( "Unimplemented" );
    }
    virtual void close( )
    {
        throw AnnotatedException( "Unimplemented" );
    }
    virtual size_t tellg( )
    {
        throw AnnotatedException( "Unimplemented" );
    }
    virtual char peek( )
    {
        throw AnnotatedException( "Unimplemented" );
    }
    virtual void safeGetLine( std::string& t )
    {
        throw AnnotatedException( "Unimplemented" );
    }
}; // class

class StdFileStream : public FileStream
{
  private:
    std::ifstream xStream;

  public:
    StdFileStream( fs::path sFilename ) : xStream( sFilename )
    {} // constructor

    bool eof( ) const
    {
        return !xStream.good( ) || xStream.eof( );
    } // method

    ~StdFileStream( )
    {
        DEBUG( std::cout << "StdFileStream read " << uiNumLinesRead << " lines in total." << std::endl;
               if( !eof( ) ) std::cerr << "WARNING: Did abort before end of File." << std::endl; ) // DEBUG
    } // deconstrucotr

    bool is_open( ) const
    {
        return xStream.is_open( );
    } // method

    void close( )
    {
        xStream.close( );
    } // method

    char peek( )
    {
        return xStream.peek( );
    } // method

    size_t tellg( )
    {
        return xStream.tellg( );
    } // method

    /**
     * code taken from
     * https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
     * suports windows linux and mac line endings
     */
    void safeGetLine( std::string& t )
    {
        t.clear( );
        DEBUG( uiNumLinesRead++; ) // DEBUG

        // The characters in the stream are read one-by-one using a std::streambuf.
        // That is faster than reading them one-by-one using the std::istream.
        // Code that uses streambuf this way must be guarded by a sentry object.
        // The sentry object performs various tasks,
        // such as thread synchronization and updating the stream state.

        std::istream::sentry se( xStream, true );
        std::streambuf* sb = xStream.rdbuf( );

        while( true )
        {
            int c = sb->sbumpc( );
            switch( c )
            {
                case '\n':
                    return;
                case '\r':
                    if( sb->sgetc( ) == '\n' )
                        sb->sbumpc( );
                    return;
                case std::streambuf::traits_type::eof( ):
                    // Also handle the case when the last line has no line ending
                    if( t.empty( ) )
                        xStream.setstate( std::ios::eofbit );
                    return;
                default:
                    t += (char)c;
            } // switch
        } // while
    } // method
}; // class

class StringStream : public FileStream
{
  private:
    std::stringstream xStream;

  public:
    StringStream( std::string& sString ) : xStream( sString )
    {} // constructor

    bool eof( ) const
    {
        return !xStream.good( ) || xStream.eof( );
    } // method

    ~StringStream( )
    {
        DEBUG( std::cout << "StringStream read " << uiNumLinesRead << " lines in total." << std::endl;
               if( !eof( ) ) std::cerr << "WARNING: Did abort before end of File." << std::endl; ) // DEBUG
    } // deconstrucotr

    bool is_open( ) const
    {
        return true;
    } // method

    void close( )
    {
        // we do not need to do anything
    } // method

    char peek( )
    {
        return xStream.peek( );
    } // method

    size_t tellg( )
    {
        return xStream.tellg( );
    } // method

    /**
     * code taken from
     * https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
     * suports windows linux and mac line endings
     */
    void safeGetLine( std::string& t ) // @todo duplicate code
    {
        t.clear( );
        DEBUG( uiNumLinesRead++; ) // DEBUG

        // The characters in the stream are read one-by-one using a std::streambuf.
        // That is faster than reading them one-by-one using the std::istream.
        // Code that uses streambuf this way must be guarded by a sentry object.
        // The sentry object performs various tasks,
        // such as thread synchronization and updating the stream state.

        std::istream::sentry se( xStream, true );
        std::streambuf* sb = xStream.rdbuf( );

        while( true )
        {
            int c = sb->sbumpc( );
            switch( c )
            {
                case '\n':
                    return;
                case '\r':
                    if( sb->sgetc( ) == '\n' )
                        sb->sbumpc( );
                    return;
                case std::streambuf::traits_type::eof( ):
                    // Also handle the case when the last line has no line ending
                    if( t.empty( ) )
                        xStream.setstate( std::ios::eofbit );
                    return;
                default:
                    t += (char)c;
            } // switch
        } // while
    } // method
}; // class

#ifdef WITH_ZLIB
class GzFileStream : public FileStream
{
  private:
    // gzFile is a pointer type
    gzFile pFile;
    int lastReadReturn = 1; // 1 == last read was ok; 0 == eof; -1 == error
    char cBuff;

  public:
    GzFileStream( fs::path sFilename )
    {
        // open the file in reading and binary mode
        pFile = gzopen( sFilename.string( ).c_str( ), "rb" );
        if( pFile != nullptr )
            // fill internal buffer
            lastReadReturn = gzread( pFile, &cBuff, 1 );
    } // constructor

    bool is_open( ) const
    {
        return pFile != nullptr;
    } // method

    bool eof( ) const
    {
        return lastReadReturn != 1;
    } // method

    ~GzFileStream( )
    {
        DEBUG( std::cout << "GzFileStream read " << uiNumLinesRead << " lines in total." << std::endl;
               if( !eof( ) ) std::cerr << "WARNING: Did abort before end of File." << std::endl; ) // DEBUG
    } // deconstrucotr

    void close( )
    {
        gzclose( pFile );
    } // method

    size_t tellg( )
    {
        return gzoffset( pFile );
    } // method

    char peek( )
    {
        return cBuff;
    } // method

    void safeGetLine( std::string& t )
    {
        t.clear( );
        DEBUG( uiNumLinesRead++; ) // DEBUG

        while( true )
        {
            // gzread returns the number of bytes read; if we do not get 1 here something went wrong
            if( lastReadReturn != 1 )
                break;
            // if we reached the end of the line return
            if( cBuff == '\r' )
                // read the \n that should be the next character
                lastReadReturn = gzread( pFile, &cBuff, 1 );
            if( cBuff == '\n' )
                break;
            t += cBuff;

            // read one char from the stream
            lastReadReturn = gzread( pFile, &cBuff, 1 );
        } // while

        // read one char from the stream
        lastReadReturn = gzread( pFile, &cBuff, 1 );
    } // method
}; // class
#endif

class Reader
{
  public:
    virtual size_t getCurrPosInFile( ) const = 0;
    virtual size_t getFileSize( ) const = 0;
    virtual size_t getCurrFileIndex( ) const = 0;
    virtual size_t getNumFiles( ) const = 0;
}; // class


class SingleFileReader : public Module<NucSeq, true>, public Reader
{
  public:
    virtual const std::string status( ) const
    {
        throw std::runtime_error( "This method must be overridden." );
    } // method
}; // class

/**
 * @brief Reads Queries from a file.
 * @details
 * Reads (multi-)fasta or fastaq format.
 */
class FileReader : public SingleFileReader
{
  public:
    const fs::path xFileName;
    std::shared_ptr<FileStream> pFile;
    size_t uiFileSize = 0;
    DEBUG( size_t uiNumLinesWithNs = 0; ) // DEBUG

  public:
    /**
     * @brief creates a new FileReader.
     */
    FileReader( fs::path sFileName ) : xFileName( sFileName )
    {
#ifdef WITH_ZLIB
        if( sFileName.extension( ).string( ) == ".gz" )
            pFile = std::make_shared<GzFileStream>( sFileName );
        else
#endif
            pFile = std::make_shared<StdFileStream>( sFileName );
        if( !pFile->is_open( ) )
        {
            throw AnnotatedException( "Unable to open file" + sFileName.string( ) );
        } // if
        std::ifstream xFileEnd( sFileName, std::ifstream::ate | std::ifstream::binary );
        uiFileSize = xFileEnd.tellg( );
        if( uiFileSize == 0 )
            this->setFinished( );
    } // constructor

    /**
     * @brief creates a new FileReader.
     */
    FileReader( const ParameterSetManager& rParameters, fs::path sFileName ) : FileReader( sFileName )
    {} // constructor

    /**
     * @brief creates a new FileReader.
     */
    FileReader( const ParameterSetManager& rParameters, std::string& sString, size_t uiStringSize ) : xFileName( )
    {
        pFile = std::make_shared<StringStream>( sString );
        uiFileSize = uiStringSize;
        if( uiFileSize == 0 )
            throw AnnotatedException( "Got empty query via text input." );
    } // constructor

    ~FileReader( )
    {
        pFile->close( );
    } // deconstructor

    std::shared_ptr<NucSeq> EXPORTED execute( );

    // @override
    virtual bool requiresLock( ) const
    {
        return true;
    } // function

    size_t getCurrPosInFile( ) const
    {
        if( pFile->eof( ) )
            return uiFileSize;
        return pFile->tellg( );
    } // function

    size_t getFileSize( ) const
    {
        // prevent floating point exception here (this is only used for progress bar...)
        if( uiFileSize == 0 )
            return 1;
        return uiFileSize;
    } // function

    size_t getCurrFileIndex( ) const
    {
        return 0;
    } // method

    size_t getNumFiles( ) const
    {
        return 1;
    } // method

    virtual const std::string status( ) const
    {
        return std::string( xFileName.string( ) )
            .append( ":" )
            .append( std::to_string( this->getCurrPosInFile( ) ) )
            .append( "/" )
            .append( std::to_string( this->getFileSize( ) ) );
    } // method
}; // class

class FileListReader : public SingleFileReader
{
  public:
    size_t uiFileIndex;
    std::vector<fs::path> vsFileNames;
    std::unique_ptr<FileReader> pFileReader;

  private:
    void openNextFileWhileNecessary( )
    {
        while( pFileReader->isFinished( ) && !this->isFinished( ) )
        {
            uiFileIndex++;
            if( uiFileIndex >= vsFileNames.size( ) )
                this->setFinished( );
            else
                pFileReader = std::make_unique<FileReader>( vsFileNames[ uiFileIndex ] );
        } // if
    } // method
  public:
    /**
     * @brief creates a new FileReader.
     */
    FileListReader( const std::vector<fs::path>& vsFileNames )
        : uiFileIndex( 0 ),
          vsFileNames( vsFileNames ),
          pFileReader( std::make_unique<FileReader>( vsFileNames[ uiFileIndex ] ) )
    {
        openNextFileWhileNecessary( );
    } // constructor

    FileListReader( const ParameterSetManager& rParameters, const std::vector<fs::path>& vsFileNames )
        : uiFileIndex( 0 ),
          vsFileNames( vsFileNames ),
          pFileReader( std::make_unique<FileReader>( vsFileNames[ uiFileIndex ] ) )
    {
        openNextFileWhileNecessary( );
    } // constructor

    std::shared_ptr<NucSeq> EXPORTED execute( )
    {
        auto pRet = pFileReader->execute( );
        openNextFileWhileNecessary( );
        return pRet;
    } // method

    // @override
    virtual bool requiresLock( ) const
    {
        return pFileReader->requiresLock( );
    } // function

    size_t getCurrPosInFile( ) const
    {
        return pFileReader->getCurrPosInFile( );
    } // function

    size_t getFileSize( ) const
    {
        return pFileReader->getFileSize( );
    } // function

    size_t getCurrFileIndex( ) const
    {
        return uiFileIndex;
    } // method

    size_t getNumFiles( ) const
    {
        return vsFileNames.size( );
    } // method

    virtual const std::string status( ) const
    {
        return std::string( std::to_string( this->getCurrFileIndex( ) ) )
            .append( "/" )
            .append( std::to_string( this->getNumFiles( ) ) )
            .append( "::" )
            .append( this->pFileReader->status( ) );
    } // method
}; // class

typedef ContainerVector<std::shared_ptr<NucSeq>> TP_PAIRED_READS;
/**
 * @brief Reads Queries from a file.
 * @details
 * Reads (multi-)fasta or fastaq format.
 */
class PairedFileReader : public Module<TP_PAIRED_READS, true>, public Reader
{
  public:
    std::shared_ptr<SingleFileReader> pF1;
    std::shared_ptr<SingleFileReader> pF2;

    /**
     * @brief creates a new FileReader.
     */
    PairedFileReader(
        const ParameterSetManager& rParameters, std::vector<fs::path> vsFileName1, std::vector<fs::path> vsFileName2 )
        : pF1( std::make_shared<FileListReader>( vsFileName1 ) ), pF2( std::make_shared<FileListReader>( vsFileName2 ) )
    {} // constructor
    /**
     * @brief creates a new FileReader.
     */
    PairedFileReader( const ParameterSetManager& rParameters, std::string sQuery, std::string sMate )
        : pF1( std::make_shared<FileReader>( rParameters, sQuery, sQuery.size( ) ) ),
          pF2( std::make_shared<FileReader>( rParameters, sMate, sMate.size( ) ) )
    {} // constructor

    std::shared_ptr<TP_PAIRED_READS> EXPORTED execute( );

    // @override
    virtual bool requiresLock( ) const
    {
        return pF1->requiresLock( ) || pF2->requiresLock( );
    } // function

    size_t getCurrPosInFile( ) const
    {
        return pF1->getCurrPosInFile( ) + pF2->getCurrPosInFile( );
    } // function

    size_t getFileSize( ) const
    {
        return pF1->getFileSize( ) + pF2->getFileSize( );
    } // function

    size_t getCurrFileIndex( ) const
    {
        return pF1->getCurrFileIndex( ) + pF2->getCurrFileIndex( );
    } // method

    size_t getNumFiles( ) const
    {
        return pF1->getNumFiles( ) + pF2->getNumFiles( );
    } // method

    virtual const std::string status( ) const
    {
        return std::string( this->pF1->status( ) ).append( ";" ).append( this->pF2->status( ) );
    } // method
}; // class

} // namespace libMA

#ifdef WITH_PYTHON
#ifdef WITH_BOOST
void exportFileReader( );
#else
void exportFileReader( py::module& rxPyModuleId );
#endif
#endif

#endif