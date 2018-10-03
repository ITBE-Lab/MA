/**
 * @file fileReader.h
 * @brief Reads queries from a file.
 * @author Markus Schmidt
 */
#ifndef FILE_READER_H
#define FILE_READER_H

#include "container/nucSeq.h"
#include "module/module.h"


/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <fstream>
/// @endcond

namespace libMA
{
class Reader : public Module<NucSeq, true>
{
  public:
    virtual size_t getCurrPosInFile( ) const = 0;
    virtual size_t getFileSize( ) const = 0;
}; // class

/**
 * @brief Reads Queries from a file.
 * @details
 * Reads (multi-)fasta or fastaq format.
 */
class FileReader : public Reader
{
  public:
    std::shared_ptr<std::ifstream> pFile;
    size_t uiFileSize = 0;
    DEBUG( size_t uiNumLinesRead = 0; size_t uiNumLinesWithNs = 0; ) // DEBUG
    std::shared_ptr<NucSeq> EXPORTED execute( );
    /**
     * @brief creates a new FileReader.
     */
    FileReader( std::string sFileName ) : pFile( new std::ifstream( sFileName ) )
    {
        if( !pFile->is_open( ) )
        {
            throw AnnotatedException( "Unable to open file" + sFileName );
        } // if
        std::ifstream xFileEnd( sFileName, std::ifstream::ate | std::ifstream::binary );
        uiFileSize = xFileEnd.tellg( );
        if( uiFileSize == 0 )
            std::cerr << "Warning: empty file: " << sFileName << std::endl;
    } // constructor

    ~FileReader( )
    {
        DEBUG( std::cout << "read " << uiNumLinesRead << " lines in total." << std::endl;
               std::cout << "read " << uiNumLinesWithNs << " N's." << std::endl;
               if( !pFile->eof( ) ) std::cerr << "WARNING: Did abort before end of File." << std::endl; ) // DEBUG
        pFile->close( );
    } // deconstructor

    // @override
    virtual bool requiresLock( ) const
    {
        return true;
    } // function

    size_t getCurrPosInFile( ) const
    {
        if( !pFile->good( ) || pFile->eof( ) )
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
  private:
    /**
     * code taken from
     * https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
     * suports windows linux and mac line endings
     */
    inline void safeGetline( std::string& t )
    {
        t.clear( );
        DEBUG( uiNumLinesRead++; ) // DEBUG

        // The characters in the stream are read one-by-one using a std::streambuf.
        // That is faster than reading them one-by-one using the std::istream.
        // Code that uses streambuf this way must be guarded by a sentry object.
        // The sentry object performs various tasks,
        // such as thread synchronization and updating the stream state.

        std::istream::sentry se( *pFile, true );
        std::streambuf* sb = pFile->rdbuf( );

        for( ;; )
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
                        pFile->setstate( std::ios::eofbit );
                    return;
                default:
                    t += (char)c;
            }
        }
    } // method
}; // class

#if 0
/**
 * @brief Reads Queries from a file.
 * @details
 * Reads (multi-)fasta or fastaq format.
 */
class PairedFileReader : public Reader
{
  public:
    FileReader xF1;
    FileReader xF2;

    /**
     * @brief creates a new FileReader.
     */
    PairedFileReader( std::string sFileName1, std::string sFileName2 ) : xF1( sFileName1 ), xF2( sFileName2 )
    {
        if( xF1.getFileSize( ) != xF2.getFileSize( ) )
            std::cerr << "Paired alignment with differently sized files." << std::endl;
    } // constructor

    std::shared_ptr<Container> EXPORTED execute( std::shared_ptr<ContainerVector> vpInput );

    /**
     * @brief Used to check the input of execute.
     * @details
     * Returns:
     * - Nil
     */
    ContainerVector EXPORTED getInputType( ) const;

    /**
     * @brief Used to check the output of execute.
     * @details
     * Returns:
     * - ContainerVector(NucSeq)
     */
    std::shared_ptr<Container> EXPORTED getOutputType( ) const;

    // @override
    std::string getName( ) const
    {
        return "PairedFileReader";
    } // function

    // @override
    std::string getFullDesc( ) const
    {
        return std::string( "PairedFileReader" );
    } // function

    // @override
    bool outputsVolatile( ) const
    {
        return true;
    } // function

    // @override
    bool requiresLock( ) const
    {
        return true;
    } // function

    size_t getCurrPosInFile( ) const
    {
        return xF1.getCurrPosInFile( ) + xF2.getCurrPosInFile( );
    } // function

    size_t getFileSize( ) const
    {
        return xF1.getFileSize( ) + xF2.getFileSize( );
    } // function
}; // class
#endif

} // namespace libMA

#ifdef WITH_PYTHON
void exportFileReader( );
#endif

#endif