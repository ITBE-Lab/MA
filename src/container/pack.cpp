/** 
 * @file pack.cpp
 * @author Arne Kutzner
 */
#include "container/pack.h"
#include "container/nucSeq.h"

using namespace libMA;

//#ifdef FASTA_READER



    #define USE_STL

    #ifdef USE_STL
        typedef std::string TextSequence;
    #endif

    /**
     * @brief Structure that describes single FASTA Record.
     */
    class FastaDescriptor
    {
    public :
        std::string sName;    // name of the FASTA Record
        std::string sComment; // comment for the entry
        TextSequence qualifier; // qualifier (describes quality of a read.)

        /* Here we take a text sequence because FASTA files contain text.
        */
        std::unique_ptr<TextSequence> pSequenceRef; 

        FastaDescriptor() 
            : sName(), sComment(), qualifier(), pSequenceRef( new TextSequence() ) // std::make_shared<TextSequence>() )
        { } // default constructor
    }; // class

    /**
     * @brief StreamType can be GzipInputStream or GzipInputFileStream.
     */
    template<int BUFFER_SIZE, typename StreamType>
    class BufferedStreamer
    {
    protected :
        /* The stream used for all reading operations.
        * FIX ME: The position of the stream over here is a design flaw.
        */
        StreamType xBufferedStream;

    private :
        /* We prohibit the copy of Objects of this class.
        */
        BufferedStreamer<BUFFER_SIZE, StreamType>(const BufferedStreamer&);

        /* The buffer is placed on the heap in order to safe stack space.
        */
        std::array<char, BUFFER_SIZE> aBuffer;

        /* References to the first and last element in the buffer
        */
        size_t uiBegin = 0; 
        size_t uiEnd = 0;

        /* bIsEOF becomes true, if we reached the last block
        */
        bool bIsEOF = false;   

        inline void vRefillBuffer()
        {
            uiBegin = 0;   

            /* Fill the buffer by reading from the stream.
            */
            xBufferedStream.read( &aBuffer[0], BUFFER_SIZE );

            if ( xBufferedStream.bad() )
            {
                throw fasta_reader_exception( "Something went wrong during FASTA stream reading" );
            } // if

            /* uiEnd saves the number of characters that we could read in the stream operation
            */
            uiEnd = (size_t)xBufferedStream.gcount();

            if (uiEnd < BUFFER_SIZE)
            {
                bIsEOF = 1; 
            } // if
        } // inline method

    public :
        /* Constructor:
        * rxBufferedStreamConstructorArg is simply forwarded to the constructor of the buffered stream.
        */
        template<typename TP>
        BufferedStreamer<BUFFER_SIZE, StreamType>( TP &rxBufferedStreamConstructorArg /* const char* pcFileNameRef */ ) :
        xBufferedStream( rxBufferedStreamConstructorArg /* pcFileNameRef */ ),
            uiBegin( 0 ), 
            uiEnd( 0 ), 
            bIsEOF( false )
        {} // constructor

        template<typename TP1, typename TP2>
        BufferedStreamer<BUFFER_SIZE, StreamType>( TP1 &a, TP2 &b /* const char* pcFileNameRef */ ) :
        xBufferedStream( a, b /* pcFileNameRef */ ),
            uiBegin( 0 ), 
            uiEnd( 0 ), 
            bIsEOF( false )
        {} // constructor

        /* Virtual destructor. (interesting in the context of polymorphism.)
        */
        virtual ~BufferedStreamer()
        {} // destructor

        bool fail()
        {
            return xBufferedStream.fail();
        }

        /* keep me here, so that I stay automatically inlined ...
        */
        inline bool getSingleChar( char &cChar )                               
        {                                                                                                            
            if (uiBegin >= uiEnd) 
            {                                                   
                if ( bIsEOF == false )
                {
                    vRefillBuffer();

                    if (uiEnd == 0 )
                    {
                        return false;
                    } // if level 3
                } // if level 2
                else
                {
                    return false;
                } // else
            } // if level 1

            cChar = aBuffer[uiBegin++];

            return true;                                      
        } // inline method

        /* Reads from the stream until it sees the given delimiter.
        * Returns: true  : could find some string until delimiter
        *          false : otherwise
        */
        inline bool bGetUntilDelimiter( const bool bUntil_CR_or_LF, TextSequence &rSequence, char &rcLastCharacterRead ) 
        {                                                                                                                                       
            rcLastCharacterRead = '\0';  

            while ( true ) 
            {                                                                                                              
                size_t uiIterator;   

                if (uiBegin >= uiEnd) 
                {                                                                     
                    /* We consumed the complete buffer content
                    */
                    if ( bIsEOF == false )
                    {                                                                              
                        vRefillBuffer();
                        
                        if ( uiEnd == 0 ) 
                        {
                            /* We reached the end of file and couldn't read anything.
                            */
                            return false; 
                        }
                    }  // if
                    else
                    {
                        /* The buffer is empty and we have an EOF, so we exit the while loop and return
                        */
                        return false;  
                    } // else
                } // if                                                                                                                

                if ( bUntil_CR_or_LF ) 
                {                                                                                        
                    /* We read until we see a CR or LF
                    */
                    for ( uiIterator = uiBegin; uiIterator < uiEnd; ++uiIterator ) 
                    {
                        if ( aBuffer[uiIterator] == '\n' || aBuffer[uiIterator] == '\r' ) 
                        {
                            break;
                        } // if
                    } // for            
                } // if
                else 
                {                                                                                                       
                    /* We read until we see some space 
                    */
                    for( uiIterator = uiBegin; uiIterator < uiEnd; ++uiIterator )
                    {
                        if ( isspace( static_cast<unsigned char>( aBuffer[uiIterator] ) ) )  // CHECK: Is it possible to work with a reinterpret_cast over here?
                            break;
                    } // for
                }   
    #ifdef USE_STL
                rSequence.append( &aBuffer[uiBegin], uiIterator - uiBegin ); 
    #else
                rSequence.PlainSequence<char>::vAppend( &aBuffer[0] + uiBegin, uiIterator - uiBegin );
    #endif
                uiBegin = uiIterator + 1; 

                if( uiIterator < uiEnd ) 
                {
                    /* We found the requested delimiters or spaces.
                    * We store the last character that we got in cLastCharacterRead
                    */
                    rcLastCharacterRead = aBuffer[uiIterator];
                    
                    break;                                                                                                
                } // outer if                                                                                                                
            } // while                                                                                                                             

            return true;                                                                                                 
        } // method
    }; // class

    /**
     * @brief A stream reader for FASTA streams.
     * @details
     * A FastaStreamReader can read a FASTA-file and creates the appropriate FASTA record.
     * Requires as template argument a stream type:
     * E.g.: std::istream, std::ifstream, GzipInputFileStream, GzipInputStream
     * FIX ME: Rethink the design of the FASTA reader. The current version is not really well made.
     */
    template<typename StreamType>
    class FastaStreamReader : public BufferedStreamer<8192, StreamType>
    { 
    private:
        /* Indicates internally, whether the currently processed sequence is the first one.
        */
        bool bIsFirstSequence = true;

        /* Return values:
        >=0  length of the sequence (more sequences following)
            0 : read a sequence no more sequences following
            1 : There are more sequences following
        -1   reading failed
        -2   truncated quality string
        */                                                                                                           
        int readFastaRecord( FastaDescriptor &fastaRecord ) //, bool bIsFirstSequence )                                                                      
        {                                                                                                                                       
            char cChar;     

            if ( bIsFirstSequence ) 
            { 
                /* then jump to the next header line 
                */ 
                while ( true )
                {
                    if ( BufferedStreamer<8192, StreamType>::getSingleChar( cChar ) == false )
                    {
                        return -1;
                    } // if
                    
                    if ( cChar == '>' || cChar == '@' )
                    {
                        /* We did see the beginning of a FASTA record. 
                        */
                        break;
                    } // if
                } // while   
                bIsFirstSequence = false;
            } // if    

            /* We read the name. (String until the first space / new line)
            */
            TextSequence bufferSequence;
            if ( BufferedStreamer<8192, StreamType>::bGetUntilDelimiter( false, bufferSequence, cChar ) == false )
            {
                return -1;
            } // if

    #ifdef USE_STL
            fastaRecord.sName = bufferSequence; 
            bufferSequence.clear();
    #else
            fastaRecord.sName = bufferSequence.cString();
            bufferSequence.vClear();
    #endif

            /* We read the remaining part of the first line as comment
            */
            if ( cChar != '\n' && cChar != '\r' )
            {
                if ( BufferedStreamer<8192, StreamType>::bGetUntilDelimiter( true, bufferSequence, cChar ) == false )
                {
                    return -1;
                } // inner if
            } // if

    #ifdef USE_STL
            fastaRecord.sComment = bufferSequence; 
            bufferSequence.clear();
    #else
            /* TO DO: We have to trim the comment, because there can be LF inside in some file formats.
            */
            fastaRecord.sComment = bufferSequence.cString();
            bufferSequence.vClear();
    #endif
        
            /* We read the core sequence
            */
            while (    true ) 
            { 
                if ( BufferedStreamer<8192, StreamType>::getSingleChar( cChar ) == false )
                {
                    /* EOF there isn't anything more ...
                    */
                    return 0;
                } // if

                if ( (cChar == '>') || (cChar == '@') )
                {
                    /* We found the beginning of a next sequence
                    * '>' == FASTA, '@' == FASTQ
                    */
                    return 1;
                } // if

                if ( cChar == '+' )
                {
                    /* We found a qualifier ...
                    */
                    break;
                }
                
                if ( isgraph( cChar ) ) 
                { 
                    /* We found a printable non-space character and 
                    * append the single character to the sequence
                    */
    #ifdef USE_STL
                    fastaRecord.pSequenceRef->push_back( cChar ); 
    #else
                    fastaRecord.pSequenceRef->vAppend( cChar );
    #endif
                } // if                                                                                                                 
            } // while
            
            /* skip the rest of the '+' line 
            */
    #ifdef USE_STL
            bufferSequence.clear(); // <- bufferSequence.vClear();
    #else
            bufferSequence.vClear();
    #endif
            BufferedStreamer<8192, StreamType>::bGetUntilDelimiter( true, bufferSequence, cChar );

            /* TO OD; Fore a reserve for qualifier using the size of the sequence
            */
            while (    true ) 
            { 
                if ( BufferedStreamer<8192, StreamType>::getSingleChar( cChar ) == false )
                {
                    /* EOF there isn't anything more ...
                    */
                    break;
                } // if

                if ( cChar >= 33 && cChar <= 127 ) 
                {
    #ifdef USE_STL
                    fastaRecord.qualifier.push_back( cChar ); 
    #else
                    fastaRecord.qualifier.vAppend( cChar );
    #endif
                } // if
            } // while

    #ifdef USE_STL        
            if ( fastaRecord.pSequenceRef->length() != fastaRecord.qualifier.length() ) // <- if ( fastaRecord.pSequenceRef->uxGetSequenceSize() != fastaRecord.qualifier.uxGetSequenceSize() )
    #else
            if ( fastaRecord.pSequenceRef->uxGetSequenceSize() != fastaRecord.qualifier.uxGetSequenceSize() )
    #endif
            {
                return -2;
            }  // if  

            return 0;
        } //method

    public :
        /* Apply the function func to all FASTA sequences defined in the given FASTA stream.
        */
        template <typename FunctionType>
        void forAllSequencesDo( FunctionType &&function )
        {
            int iContinue = 1;

            while( iContinue > 0 )
            {
                FastaDescriptor xFastaRecord; // FIX ME: Analyze whether it is better to always reuse the same record and to clear instead.
                
                /* We read a single record.
                */
                iContinue = readFastaRecord( xFastaRecord );
                function( xFastaRecord );
            } // while
        } // generic function

        /* Sequence provider gets access to the private stuff! 
        */
        friend class FastaReader;

        /* Deprecated. Remove me.
        */
        FastaStreamReader( std::istream &xInputStream, FastaDescriptor &fastaRecord ) : 
            BufferedStreamer<8192, StreamType>( xInputStream ), 
            bIsFirstSequence( true )
        {
            /* We read a single record.
            */
            readFastaRecord( fastaRecord );
        } // constructor

        /* Constructor:
        * The current design always creates a separated stream for all reading. 
        * This stream is invisible to the outside and constructed as part of the buffer.
        * TP must be: 
        * std::istream for GzipInputStream (std::istream is the source for GzipInputStream)
        * std::istream for std::istream     
        */
        template<typename TP> 
        FastaStreamReader( TP &rxBufferedStreamConstructorArg ) : 
            BufferedStreamer<8192, StreamType>( rxBufferedStreamConstructorArg )
        {} // constructor

        template<typename TP1, typename TP2> 
        FastaStreamReader( TP1 &a, TP2 &b ) : 
            BufferedStreamer<8192, StreamType>( a, b )
        {} // constructor

        virtual ~FastaStreamReader()
        {
            /* Automatic call of the superclass destructor
            */
        } // destructor    
    }; // class

    /**
     * @brief a FASTA rader for file streams.
     * @details
     * Use this a FASTA rader for file streams. (GzipInputFileStream and std::ifstream)
     * The constructor checks, whether the file could be successfully opened.
     */
    template<typename StreamType>
    class FastaFileStreamReader : public FastaStreamReader<StreamType>
    {
    public :
        /* The stream used for all reading operations.
        */
        FastaFileStreamReader( const std::string &rsFileName ) :
            FastaStreamReader<StreamType>( rsFileName )
        {} // constructor
    }; // class FastaFileStreamReader

    /**
     * @brief Reads a single FASTA record from a file
     */
    class FastaReader : public FastaDescriptor
    {                                                                                    
    public :
        /* The default constructor.
        */
        FastaReader()
        { } 

        /* The central load function.
        */
        void vLoadFastaFile( const char *pcFileName )
        {
            /* If the FASTA reader experiences some problem, it will throw an exception.
            */
           auto xMode = std::ios::in | std::ios::binary;
            FastaStreamReader<std::ifstream> fastaReader( 
                pcFileName, xMode
            );
            
        } // method Anonymous

        /* The destruction will free the memory of the FASTA-Record
        */
        ~FastaReader()
        { } // destructor
    }; // class


    void Pack::vAppendFastaSequence( const FastaDescriptor &rxFastaDescriptor ) 
    {
        /* This is a bit inefficient. We could boost performance by allowing a move for the sequence
        */
        if(rxFastaDescriptor.pSequenceRef->length() == 0)
            return;
        //if(rxFastaDescriptor.pSequenceRef->length() > 50000)
        //    return;

        vAppendSequence( 
                rxFastaDescriptor.sName,        // Name of the embedded sequence
                rxFastaDescriptor.sComment,    // Comments for the sequence
                NucSeq( *rxFastaDescriptor.pSequenceRef )
            );
    } // method    

    /* Entry point, for the construction of packs.
    * pcPackPrefix is some prefix for the pack-files.
    * Reads all sequences on the file system and creates a sequence collection out of them.
    */
    void Pack::vAppendFASTA( const std::string &sFastaFilePath )
    {   
        /* Open a stream for FASTA File reading.
         * Open a FASTA Stream reader using the input stream. 
         */
        auto xMode = std::ios::in | std::ios::binary;
        FastaStreamReader<std::ifstream> xFastaReader( 
                sFastaFilePath, 
                xMode
            );

        /* We check, whether file opening worked well.
        */
        if ( xFastaReader.fail() )
        {    /* Something is wrong with respect to the input-file
            */
            throw fasta_reader_exception( "File open error." );
        } // if

        
        /* Apply the lambda expression to all records in the file
        */
        xFastaReader.forAllSequencesDo
        (
            [&] ( const FastaDescriptor &xFastaRecord ) 
            {    /* Add the current FASTA sequence to the pack.
                */
                //BOOST_LOG_TRIVIAL( trace ) << "Add sequence " << xFastaRecord.sName;
                vAppendFastaSequence( xFastaRecord );
            } // lambda
        ); // function call
    } // method

    void Pack::vAppendFastaFile( const char *pcFileName ) 
    {
        vAppendFASTA(pcFileName);
    } // method


#ifdef WITH_PYTHON
void exportPack()
{
    boost::python::class_<
            Pack, 
            boost::noncopyable,
            boost::python::bases<Container>,
            std::shared_ptr<Pack>
        >(
                "Pack",
                "unpacked_size_single_strand: the size of one strand\n"
                "\n"
                "Holds a packed sequence.\n"
                "Usefull for long sequences.\n"
            )
        .def(
                "unpacked_size", 
                &Pack::uiUnpackedSizeForwardPlusReverse,
                "arg1: self\n"
                "returns: the length of the sequence within the pack.\n"
            )
        .def(
                "append", 
                &Pack::vAppendSequence_boost,
                "arg1: self\n"
                "arg2: a NucSeq\n"
                "returns: nil\n"
                "\n"
                "Appends seq at the end of the pack.\n"
            )
#if 1
        .def(
                "append_fasta_file", 
                &Pack::vAppendFastaFile,
                "arg1: self\n"
                "arg2: the filename\n"
                "returns: nil\n"
                "\n"
                "Appends seq at the end of the pack.\n"
            )
#endif
        .def(
                "store", 
                &Pack::vStoreCollection,
                "arg1: self\n"
                "arg2: the folder and filename on disk\n"
                "returns: nil\n"
                "\n"
                "Stores this pack at the given location.\n"
            )
        .def(
                "exists", 
                &Pack::packExistsOnFileSystem,
                "arg1: self\n"
                "arg2: the folder and filename on disk\n"
                "returns: a bool indicating if the file exists\n"
                "\n"
                "Checks weather a pack exists at the given location.\n"
            )
        .staticmethod("exists")
        .def(
                "load", 
                &Pack::vLoadCollection,
                "arg1: self\n"
                "arg2: the folder and filename on disk\n"
                "returns: nil\n"
                "\n"
                "Loads a pack from the given location on disk.\n"
            )
        .def(
                "extract_from_to", 
                &Pack::vExtract,
                "arg1: self\n"
                "arg2: begin of extraction\n"
                "arg3: end of extraction\n"
                "returns: the extracted sequence as NucSeq\n"
                "\n"
                "Extracts a sequence from the pack.\n"
                "Indices are inclusive.\n"
            )
        .def(
                "extract_complete", 
                &Pack::vColletionAsNucSeq,
                "arg1: self\n"
                "returns: the extracted sequence as NucSeq\n"
                "\n"
                "Extracts the entire pack as sequence.\n"
            )
        .def(
                "extract_forward_strand", 
                &Pack::vColletionWithoutReverseStrandAsNucSeq,
                "arg1: self\n"
                "returns: the extracted sequence as NucSeq\n"
                "\n"
                "Extracts the forward strand of the pack as sequence.\n"
            )
        .def(
                "extract_reverse_strand", 
                &Pack::vColletionOnlyReverseStrandAsNucSeq,
                "arg1: self\n"
                "returns: the extracted sequence as NucSeq\n"
                "\n"
                "Extracts the reverse strand of the pack as sequence.\n"
            )
        .def(
                "is_bridging", 
                &Pack::bridgingSubsection_boost
            )
        .def(
                "printHoles", 
                &Pack::printHoles
            )
        .def(
                "start_of_sequence", 
                &Pack::startOfSequenceWithName
            )
        .def(
                "length_of_sequence", 
                &Pack::lengthOfSequenceWithName
            )
        .def(
                "start_of_sequence_id", 
                &Pack::startOfSequenceWithId
            )
        .def(
                "name_of_sequence", 
                &Pack::nameOfSequenceForPosition
            )
        .def_readonly(
                "unpacked_size_single_strand", 
                &Pack::uiUnpackedSizeForwardStrand
            )
        ;

    //tell boost python that pointers of these classes can be converted implicitly
    boost::python::implicitly_convertible<
            std::shared_ptr<Pack>,
            std::shared_ptr<Container>
        >();
}//function
#endif