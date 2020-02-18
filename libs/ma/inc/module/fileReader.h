/**
 * @file fileReader.h
 * @brief Reads queries from a file.
 * @author Markus Schmidt
 */
#ifndef FILE_READER_H
#define FILE_READER_H

#include "container/nucSeq.h"
#include "module/module.h"
#include "support.h"


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
        throw std::runtime_error( "This function should have been overridden" );
    }

    virtual bool is_open( ) const
    {
        throw std::runtime_error( "This function should have been overridden" );
    }
    virtual void close( )
    {
        throw std::runtime_error( "This function should have been overridden" );
    }
    virtual size_t tellg( )
    {
        throw std::runtime_error( "This function should have been overridden" );
    }
    virtual char peek( )
    {
        throw std::runtime_error( "This function should have been overridden" );
    }
    virtual void safeGetLine( std::string& t )
    {
        throw std::runtime_error( "This function should have been overridden" );
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
        DEBUG( if( !eof( ) ) std::cout << "StdFileStream read " << uiNumLinesRead << " lines in total." << std::endl;
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
        DEBUG( if( !eof( ) ) std::cout << "StringStream read " << uiNumLinesRead << " lines in total." << std::endl;
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
    unsigned char cBuff;

  public:
    GzFileStream( fs::path sFilename )
    {
        // open the file in reading and binary mode
        pFile = gzopen( sFilename.string( ).c_str( ), "rb" );
        if( pFile != nullptr )
            // fill internal buffer
            lastReadReturn = gzread( pFile, &cBuff, 1 );
        else
            lastReadReturn = -1; // if the file could net be opened set error immediately
    } // constructor

    bool eof( ) const
    {
        return lastReadReturn != 1;
    } // method

    ~GzFileStream( )
    {
        DEBUG( if( !eof( ) ) std::cout << "StdFileStream read " << uiNumLinesRead << " lines in total." << std::endl;
               if( !eof( ) ) std::cerr << "WARNING: Did abort before end of File." << std::endl; ) // DEBUG
    } // deconstructor

    bool is_open( ) const
    {
        return pFile != nullptr;
    } // method

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

        // read until EoF or EoL is reached
        while( true )
        {
            // gzread returns the number of bytes read; if we do not get 1 here something went wrong
            if( lastReadReturn != 1 )
                break;
            // if we reached the end of the line return
            if( cBuff == '\r' )
                // read the \n that should be the next character
                lastReadReturn = gzread( pFile, &cBuff, 1 );
            if( lastReadReturn != 1 || cBuff == '\n' )
                break;
            t += cBuff;

            // read one char from the stream
            lastReadReturn = gzread( pFile, &cBuff, 1 );
        } // while

        // read the EoL char from the stream (if the stream is ok)
        if( lastReadReturn == 1 )
            lastReadReturn = gzread( pFile, &cBuff, 1 );
    } // method
}; // class
#endif

class Reader
{
  public:
    virtual size_t getCurrPosInFile( ) = 0;
    virtual size_t getFileSize( ) = 0;
    virtual size_t getCurrFileIndex( ) = 0;
    virtual size_t getNumFiles( ) const = 0;
}; // class

class FileReader;

class SingleFileReader : public Module<NucSeq, true>, public Reader
{
  public:
    virtual const std::string status( )
    {
        throw std::runtime_error( "This method must be overridden." );
    } // method

    virtual void reset( )
    {
        throw std::runtime_error( "This method must be overridden." );
    } // method

    virtual void withFileReader( std::function<void( FileReader& )> func )
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
            throw std::runtime_error( "Unable to open file" + sFileName.string( ) );
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
            throw std::runtime_error( "Got empty query via text input." );
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

    size_t getCurrPosInFile( )
    {
        if( pFile->eof( ) )
            return uiFileSize;
        return pFile->tellg( );
    } // function

    size_t getFileSize( )
    {
        // prevent floating point exception here (this is only used for progress bar...)
        if( uiFileSize == 0 )
            return 1;
        return uiFileSize;
    } // function

    size_t getCurrFileIndex( )
    {
        return 0;
    } // method

    size_t getNumFiles( ) const
    {
        return 1;
    } // method

    virtual const std::string status( )
    {
        return std::string( xFileName.string( ) )
            .append( ":" )
            .append( std::to_string( this->getCurrPosInFile( ) ) )
            .append( "/" )
            .append( std::to_string( this->getFileSize( ) ) );
    } // method

    virtual void reset( )
    {
        pFile->close( );
#ifdef WITH_ZLIB
        if( xFileName.extension( ).string( ) == ".gz" )
            pFile = std::make_shared<GzFileStream>( xFileName );
        else
#endif
            pFile = std::make_shared<StdFileStream>( xFileName );
    } // method

    virtual void withFileReader( std::function<void( FileReader& )> func )
    {
        func( *this );
    } // method

}; // class

/**
 * @brief reads all files in a list
 * @details
 * Applies a round robin approach instead of reading the files sequentially.
 * That way several threads can work simultaneously if multiple input files are given.
 */
class FileListReader : public SingleFileReader
{
  public:
    std::mutex xMutex;
    std::condition_variable xCondition; // For wait, notify synchronization purposes.
    std::queue<std::shared_ptr<FileReader>> xFileReaders;
    size_t uiNumFiles = 0;
    bool bIsFinished = false;

    /** @brief gets std::unique_ptr<FileReader>
     * @details
     * This function is threadsave, it uses a queue to distribute the file readers among threads.
     * It makes sure that each file reader will only ever be used by one thread.
     * @todo move this to container
     */
    virtual void withFileReader( std::function<void( FileReader& )> func )
    {
        std::shared_ptr<FileReader> pReader;
        {
            std::unique_lock<std::mutex> xLock( xMutex );
            while( xFileReaders.empty( ) )
            {
                // if some other thread finished this module,
                if( bIsFinished )
                    return;
                xCondition.wait( xLock );
            } // while
            pReader = xFileReaders.front( );
            xFileReaders.pop( );
        } // scope for xLock
        // parallel execution of func allowed
        func( *pReader );
        {
            std::unique_lock<std::mutex> xLock( xMutex );
            if( !pReader->isFinished( ) )
                xFileReaders.push( pReader );

            if( xFileReaders.empty( ) )
            {
                bIsFinished = true;
                xCondition.notify_all( );
            } // if
            else
                xCondition.notify_one( );
        } // scope for xLock
    } // method

    /**
     * @brief creates a new FileReader.
     */
    FileListReader( const std::vector<fs::path>& vsFileNames ) : uiNumFiles( vsFileNames.size( ) )
    {
        for( auto sFileName : vsFileNames )
            xFileReaders.push( std::make_shared<FileReader>( sFileName ) );
    } // constructor

    FileListReader( const ParameterSetManager& rParameters, const std::vector<fs::path>& vsFileNames )
        : FileListReader( vsFileNames )
    {} // constructor

    FileListReader( const FileListReader& ) = delete; // no copies

    std::shared_ptr<NucSeq> EXPORTED execute( )
    {
        std::shared_ptr<NucSeq> pRet = std::make_shared<NucSeq>(); // @todo nullptr
        withFileReader( [&pRet]( FileReader& xReader ) { pRet = xReader.execute( ); } );
        return pRet;
    } // method

    std::vector<std::shared_ptr<NucSeq>> read_all( )
    {
        std::vector<std::shared_ptr<NucSeq>> vRet;

        while( !isFinished( ) )
            vRet.push_back( execute( ) );

        return vRet;
    } // method

    // @override
    virtual bool requiresLock( ) const
    {
        return false; // locking mechanism is provided internally
    } // function

    size_t getCurrPosInFile( )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        return xFileReaders.front( )->getCurrPosInFile( );
    } // function

    size_t getFileSize( )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        return xFileReaders.front( )->getFileSize( );
    } // function

    size_t getCurrFileIndex( )
    {
        std::unique_lock<std::mutex> xLock( xMutex );
        return uiNumFiles - xFileReaders.size( );
    } // method

    size_t getNumFiles( ) const
    {
        return uiNumFiles;
    } // method

    virtual const std::string status( )
    {
        std::string sOut = std::string( std::to_string( this->getCurrFileIndex( ) ) )
                               .append( "/" )
                               .append( std::to_string( this->getNumFiles( ) ) )
                               .append( "::" );
        std::unique_lock<std::mutex> xLock( xMutex );
        sOut.append( xFileReaders.front( )->status( ) );
        return sOut;
    } // method

    virtual bool isFinished( )
    {
        // need to obtain lock in order to know if reader is finished
        std::unique_lock<std::mutex> xLock( xMutex );
        return bIsFinished;
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
    const bool bRevCompMate;

    /**
     * @brief creates a new FileReader that reads from a list of files.
     */
    PairedFileReader(
        const ParameterSetManager& rParameters, std::vector<fs::path> vsFileName1, std::vector<fs::path> vsFileName2 )
        : pF1( std::make_shared<FileListReader>( vsFileName1 ) ),
          pF2( std::make_shared<FileListReader>( vsFileName2 ) ),
          bRevCompMate( rParameters.getSelected( )->xRevCompPairedReadMates->get( ) )
    {} // constructor

    /**
     * @brief creates a new FileReader with a single query
     */
    PairedFileReader( const ParameterSetManager& rParameters, std::string sQuery, std::string sMate )
        : pF1( std::make_shared<FileReader>( rParameters, sQuery, sQuery.size( ) ) ),
          pF2( std::make_shared<FileReader>( rParameters, sMate, sMate.size( ) ) ),
          bRevCompMate( rParameters.getSelected( )->xRevCompPairedReadMates->get( ) )
    {} // constructor

    void EXPORTED checkPaired( );

    std::shared_ptr<TP_PAIRED_READS> EXPORTED execute( );

    // @override
    virtual bool requiresLock( ) const
    {
        return pF1->requiresLock( ) || pF2->requiresLock( );
    } // function

    size_t getCurrPosInFile( )
    {
        return pF1->getCurrPosInFile( ) + pF2->getCurrPosInFile( );
    } // function

    size_t getFileSize( )
    {
        return pF1->getFileSize( ) + pF2->getFileSize( );
    } // function

    size_t getCurrFileIndex( )
    {
        return pF1->getCurrFileIndex( ) + pF2->getCurrFileIndex( );
    } // method

    size_t getNumFiles( ) const
    {
        return pF1->getNumFiles( ) + pF2->getNumFiles( );
    } // method

    virtual const std::string status( )
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