#include <iomanip> // for std::setw

#include "util/export.h"
#include "util/parameter.h"

#undef snprintf // required under MSVC
// Documentation: https://github.com/nlohmann/json
#include "json.hpp" // FIXME: Move to different location

using namespace libMA;
using namespace libMS;

using nlohmann::json; // might be better: namespace json = nlohmann;

/* Manages Pack and FMD-Index of genomes */
class GenomeManager
{
  public:
    typedef decltype( makePledge<Pack>( std::string( ) ) ) PackPledgeType;
    typedef decltype( makePledge<FMIndex>( std::string( ) ) ) FMDIndexPledgeType;

  private:
    // Genome-prefix initially empty (equal to no genome selected)
    std::string sPackPrefix;
    std::string sGenomeName;

    // shared pointer to pack (pledge)
    PackPledgeType pxPackPledge;
    // shared pointer to FMD-index (pledge)
    FMDIndexPledgeType pxFMDIndexPledge;

  public:
    /* Constructor */
    GenomeManager( void ) : sPackPrefix( "" ), sGenomeName( "" ), pxPackPledge( nullptr ), pxFMDIndexPledge( nullptr )
    {} // constructor

    /* Getter for pack pledge */
    PackPledgeType getPackPledge( void )
    {
        if( this->pxPackPledge )
            return this->pxPackPledge;
        else
            // Avoid null-pointer problems and throw exception
            throw std::runtime_error( "Genome undetermined (Missing Pack)." );
    } // method

    /* Getter for FMD-Index pledge */
    FMDIndexPledgeType getFMDIndexPledge( void )
    {
        if( this->pxFMDIndexPledge )
            return this->pxFMDIndexPledge;
        else
            // Avoid null-pointer problems and throw exception
            throw std::runtime_error( "Genome undetermined (Missing FMD-Index)." );
    } // method

    /* Loads genome using info given in JSON-file
     * Returns empty string if loading went well
     */
    std::string loadGenome( const fs::path& rsJsonFilePath ) // folder containing json, pack and fmd-index
    {
        // read a JSON file
        try
        {
            // Read JSON with genome info
            std::ifstream xInputStream( rsJsonFilePath );
            json xJSON;
            xInputStream >> xJSON;

            // Check JSON for correctness
            if( xJSON[ "type" ] != "MA Genome" )
                return ( "JSON file does not contain valid MA genome information." );
            if( !( xJSON[ "version" ][ "major" ] == 1 && xJSON[ "version" ][ "minor" ] == 0 ) )
                return ( "Wrong version of MA genome.\n(Expected major:1 minor:0)" );

            // Extract genome data from JSON ...
            this->sGenomeName = xJSON[ "name" ];
            this->sPackPrefix = ( rsJsonFilePath.parent_path( ) / std::string( xJSON[ "prefix" ] ) ).string( );

            // Get pledges for pack and FMD-index
            this->pxPackPledge = makePledge<Pack>( sPackPrefix );
            this->pxFMDIndexPledge = makePledge<FMIndex>( sPackPrefix );
            // (DEBUG) auto xDummyGenome = this->pxPackPledge->get( );
            // (DEBUG) xDummyGenome->printHoles( );
        } // try
        catch( std::exception& rxException )
        {
            // JSON loading failed
            return std::string( rxException.what( ) );
        } // catch

        return ""; // loading went fine
    } // method

    /* Create JSON that informs about genome */
    json createGenomeJSON( const std::string& rsGenomeTitle, const std::string& rsGenomePrefix )
    {
        return {{"type", "MA Genome"},
                {"version", {{"major", 1}, {"minor", 0}}},
                {"name", rsGenomeTitle},
                {"prefix", rsGenomePrefix}};
    } // method

    // std::vector<json::json> findAllGenomesInFolder( const fs::path& sFolderPath )
    // {} // method

    /* Computation of Pack and FM-Index for a genome in FASTA file */
    void makeIndexAndPackForGenome( const fs::path& rsGenomeFolderPath, // Folder for genome storage
                                    const fs::path& rsFastaFilePath, // Path to FASTA-file that contains genome
                                    const std::string& rsGenomeTitle, // Name of genome
                                    std::function<void( const std::string& )> // Feedback to the caller
                                        fCallBack )
    {
        // Genome prefix is always FASTA prefix
        auto sGenomePrefix = rsFastaFilePath.stem( );
        auto sGenomeFullPathPrefix = fs::path( rsGenomeFolderPath ) /= sGenomePrefix;
        fCallBack( std::string( "Genome prefix:\n" ) + sGenomeFullPathPrefix.string( ) + "\n" );

        // Create and store pack
        fCallBack( "Create pack and write pack to file system.\n" );
        std::shared_ptr<Pack> pxPack( new Pack( ) );
        pxPack->vAppendFASTA( rsFastaFilePath.string( ) );
        pxPack->vStoreCollection( sGenomeFullPathPrefix.string( ) );

        // Create and store FMD index
        fCallBack( "Compute FMD-index...\nImportant note: This may take long time.\n" );
        FMIndex xFMDIndex( pxPack );
        fCallBack( "Write FMD-Index to file system.\n" );
        xFMDIndex.vStoreFMIndex( sGenomeFullPathPrefix.string( ).c_str( ) );

        // Create JSON genome info
        fCallBack( "Create JSON info.\n" );

        std::ofstream xOutStream( ( fs::path( rsGenomeFolderPath ) /= rsGenomeTitle ) += ".json" );
        xOutStream << std::setw( 4 ) << createGenomeJSON( rsGenomeTitle, sGenomePrefix.string( ) ) << std::endl;

        fCallBack( "All done!\n" );
    } // method

    /* Delivers name of selected genome */
    std::string getGenomeName( void )
    {
        return this->pxPackPledge ? this->sGenomeName : "No index selected";
    } // method

    bool isReady( void )
    {
        return this->pxPackPledge && this->pxFMDIndexPledge;
    } // method
}; // class


/* Manages location of queries (reads) */
class ReadsManager
{
  public:
    // Option 1: User entered FASTA filenames
    std::vector<fs::path> vsPrimaryQueryFullFileName;
    std::vector<fs::path> vsMateQueryFullFileName;


    // Option 2: User directly entered FASTA sequences
    std::function<std::string( void )> // Feedback to the caller
        fCallBackGetPrimaryQuery = []( ) {
            throw std::runtime_error( "fCallBackGetPrimaryQuery not set." );
            return "";
        };
    std::function<std::string( void )> // Feedback to the caller
        fCallBackGetMateyQuery = []( ) {
            throw std::runtime_error( "fCallBackGetMateQuery not set." );
            return "";
        };
#if 0
    /* For paired reads (Illumina etc.) */
    void setReadsFileSource( const std::string& rsPrimaryQuery, const std::string& rsMateQuery )
    {
        this->vsPrimaryQueryFullFileName.front() = rsPrimaryQuery;
        this->vsMateQueryFullFileName.front() = rsMateQuery;
    } // method

    void setReadsFileSource( const std::string& rsPrimaryQuery )
    {
        this->sPrimaryQueryFullFileName = rsPrimaryQuery;
        this->sMateQueryFullFileName.clear( );
    } // method
#endif

    /* Delivers filename of primary reads without suffix */
    fs::path getReadsFileNameStem( void )
    {
        if( vsPrimaryQueryFullFileName.empty( ) )
            return fs::path( "Unnamed" );
        return vsPrimaryQueryFullFileName[ 0 ].stem( );
    } // method

    /* Delivers folder of primary reads */
    fs::path getReadsFolderPath( void )
    {
        if( vsPrimaryQueryFullFileName.empty( ) )
            return fs::path( "./" );
        return vsPrimaryQueryFullFileName[ 0 ].parent_path( );
    } // method

    void reset( void )
    {
        this->vsPrimaryQueryFullFileName.clear( );
        this->vsMateQueryFullFileName.clear( );
    } // method

    bool hasPrimaryPath( )
    {
        return !( this->vsPrimaryQueryFullFileName.empty( ) );
    } // method

    bool hasMatePath( )
    {
        return !( this->vsMateQueryFullFileName.empty( ) );
    } // method
}; // class


/* Manages the output of alignments.
 * FIXME: Provisional design.
 */
class OutputManager
{
  public:
    std::map<std::string, std::string> xKindsOfOutput;
    ParameterSetManager* pxParameterSetManager; // required for configuration data
    ReadsManager* pxReadsManager; // required for filename construction

    /* Constructor */
    OutputManager( ParameterSetManager* pxParameterSetManager, ReadsManager* pxReadsManager )
        : pxParameterSetManager( pxParameterSetManager ), pxReadsManager( pxReadsManager )

    {
        xKindsOfOutput.emplace( "SAM File", "unused" );
        // xKindsOfOutput.emplace( "SQLite Database", "unused" );
    } // constructor

    /* Delivers date-time string for SAM-file naming */
    std::string dateTimeString( void )
    {
        time_t xTime = time( 0 ); // get time now
#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4996 ) // disable the un-save warning for localtime
#endif
        tm* pLocalTimeNow = localtime( &xTime );
#ifdef _MSC_VER
#pragma warning( pop )
#endif

        char aBuffer[ 80 ];
        strftime( aBuffer, 80, "(%F)-(%H-%M-%S)", pLocalTimeNow );
        return std::string( aBuffer );
    } // method

    /* Creates full path for SAM output */
    std::string SAMFullFileName( void )
    {
        if( this->pxParameterSetManager->pGeneralParameterSet->xSAMOutputTypeChoice->uiSelection == 2 )
            return this->pxParameterSetManager->pGeneralParameterSet->xSAMOutputFileName->get( ).string( );
        // SAM filename generation according to parameter settings.
        auto sFullFileName = ( this->pxParameterSetManager->pGeneralParameterSet->xSAMOutputTypeChoice->uiSelection == 0
                                   ? pxReadsManager->getReadsFolderPath( )
                                   : pxParameterSetManager->pGeneralParameterSet->xSAMOutputPath->get( ) );
        ( ( ( sFullFileName /= pxReadsManager->getReadsFileNameStem( ) ) += '-' ) += this->dateTimeString( ) ) +=
            ".sam";
        return sFullFileName.string( );
    } // method
}; // class


class ExecutionContext
{
  public:
    // Major context elements
    ParameterSetManager xParameterSetManager;
    GenomeManager xGenomeManager;
    ReadsManager xReadsManager;
    OutputManager xOutputManager; // depends on xParameterSetManager, xReadsManager

    /* Constructor */
    ExecutionContext( void ) : xOutputManager( &xParameterSetManager, &xReadsManager )
    {} // constructor

    /* Computes alignments for reads.
     * Throws an exception if something goes wrong.
     */
    void
    doAlign( std::function<bool( double dPercentageProgress, int iFileNum, int iNumFilesOverall )> fProgressCallBack =
                 []( double, int, int ) { return true; },
             std::function<void( const std::string& )> fCheckCallBack =
                 []( const std::string& rS ) { std::cout << rS << std::endl; } )
    {
#if 0
        size_t uiConcurency = xParameterSetManager.pGeneralParameterSet->getNumThreads();

        // For now, we build a computational graph for each call of doAlign
        // Possible Improvement: Cache the graph ...
        std::vector<std::shared_ptr<BasePledge>> aGraphSinks;

        // For progress computation used merely.
        std::shared_ptr<Reader> pxReader;

        const std::string sSAMFileName = xOutputManager.SAMFullFileName( );

        // Build a vector of aligner instances. (aGraphSinks is a vector of final pledges)
        // Each Pledge corresponds to a thread.
        // (The number of instances is decided by the degree of concurrency) .
        // The Pledge is later given to a worker for evaluation.
        if( xParameterSetManager.getSelected( )->usesPairedReads( ) )
        {
            // Paired reads.
            // Wish AK: TP_PAIRED_WRITER -> PairedWriterType
            std::shared_ptr<TP_PAIRED_WRITER> pxPairedWriter;
            pxPairedWriter.reset(
                new PairedFileWriter( xParameterSetManager, sSAMFileName, xGenomeManager.getPackPledge( )->get( ) ) );


            if( this->xReadsManager.hasPrimaryPath( ) != this->xReadsManager.hasMatePath( ) )
                throw std::runtime_error( "Cannot combine file and text input." );

            std::shared_ptr<PairedFileReader> pxPairedFileReader;
            if( this->xReadsManager.hasPrimaryPath( ) && this->xReadsManager.hasMatePath( ) )
                pxPairedFileReader = std::make_shared<PairedFileReader>( xParameterSetManager,
                                                                         xReadsManager.vsPrimaryQueryFullFileName,
                                                                         xReadsManager.vsMateQueryFullFileName );
            else
                pxPairedFileReader = std::make_shared<PairedFileReader>( xParameterSetManager,
                                                                         xReadsManager.fCallBackGetPrimaryQuery( ),
                                                                         xReadsManager.fCallBackGetMateyQuery( ) );
            pxReader = pxPairedFileReader;

            if( xParameterSetManager.getSelected( )->xPairedCheck->get( ) )
            {
                fCheckCallBack( "checking paired reads..." );
                pxPairedFileReader->checkPaired( );
                fCheckCallBack( "done checking paired reads" );
            } // if

            auto pxPairedQueriesPledge = promiseMe( pxPairedFileReader );
            aGraphSinks = setUpCompGraphPaired( xParameterSetManager,
                                                xGenomeManager.getPackPledge( ), // Pack
                                                xGenomeManager.getFMDIndexPledge( ), // FMD index
                                                pxPairedQueriesPledge, // (for paired reads we require two queries!)
                                                pxPairedWriter, // Output writer module(output of alignments)
                                                (unsigned int)uiConcurency ); // Number of threads
        } // if
        else
        {
            // Singular (Non-Paired) reads.
            // Wish AK: TP_WRITER -> WriterType
            std::shared_ptr<TP_WRITER> pxWriter;
            pxWriter.reset(
                new FileWriter( xParameterSetManager, sSAMFileName, xGenomeManager.getPackPledge( )->get( ) ) );
            std::string sPrimaryQuery;

            std::shared_ptr<SingleFileReader> pxFileReader;
            if( xReadsManager.hasPrimaryPath( ) )
            {
                pxFileReader =
                    std::make_shared<FileListReader>( xParameterSetManager, xReadsManager.vsPrimaryQueryFullFileName );
            }
            else
            {
                sPrimaryQuery = xReadsManager.fCallBackGetPrimaryQuery( );
                pxFileReader =
                    std::make_shared<FileReader>( xParameterSetManager, sPrimaryQuery, sPrimaryQuery.size( ) );
            } // else

            pxReader = pxFileReader;
            auto pxQueriesPledge = promiseMe( pxFileReader );
            aGraphSinks = setUpCompGraph( xParameterSetManager,
                                          xGenomeManager.getPackPledge( ), // Pack
                                          xGenomeManager.getFMDIndexPledge( ), // FMD index
                                          pxQueriesPledge, // Queries
                                          pxWriter, // Output writer module(output of alignments)
                                          (unsigned int)uiConcurency ); // Number of threads
        } // else

        // Compute the actual alignments.
        // Sets the progress bar after each finished alignment.
        BasePledge::simultaneousGet( aGraphSinks,
                                     [&]( ) {
                                         int iI1 = (int)pxReader->getCurrFileIndex( );
                                         int iI2 = (int)pxReader->getNumFiles( );
                                         return fProgressCallBack( ( (double)pxReader->getCurrPosInFile( ) * 100 ) /
                                                                       (double)pxReader->getFileSize( ),
                                                                   iI1, iI2 );
                                     } // lambda
        ); // function call

#endif
    } // method

    void doAlignCallbackLess( )
    {
        this->doAlign( );
    } // method
}; // class

// exported in export.cpp
