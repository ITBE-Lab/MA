#include "util/export.h" // MA aligner interface
#include "util/parameter.h"

using namespace libMA;

/* Manages Pack and FMD-Index of genomes */
class GenomeManager
{
    // Genome-prefix initially empty (equal to no genome selected)
    std::string sPackPrefix;
    std::string sGenomeName;

  public:
    // shared pointer to pack (pledge)
    decltype( makePledge<Pack>( std::string( ) ) ) pxPackPledge;
    // shared pointer to FMD-index (pledge)
    decltype( makePledge<FMIndex>( std::string( ) ) ) pxFMDIndexPledge;

    /* sFileName must be a MAreference file
     * FIX ME: Better exception management if something goes wrong.
     * FIX ME: In GUI avoid the call with empty file-name during initializing
     */
    void setGenome( const std::string& sFileName )
    {
        if( sFileName.empty( ) )
            return;

        std::ifstream xInfile( sFileName );
        std::getline( xInfile, this->sGenomeName );
        // sPackPrefix is the full path without the final suffix
        sPackPrefix = sFileName.substr( 0, sFileName.find_last_of( '.' ) );

        this->pxPackPledge = makePledge<Pack>( sPackPrefix );
        this->pxFMDIndexPledge = makePledge<FMIndex>( sPackPrefix );

        //-- auto xDummyGenome = this->pxPackPledge->get( );
        //-- xDummyGenome->printHoles( );
    } // method

    GenomeManager( void ) : pxPackPledge( nullptr ), pxFMDIndexPledge( nullptr )
    {} // constructor

    /* Delivers name of selected genome */
    std::string getGenomeName( void )
    {
        return pxPackPledge ? sGenomeName : "No genome selected";
    } // method

    bool isReady( void )
    {
        return pxPackPledge && pxFMDIndexPledge;
    } // method
}; // class


/* Manages location of queries (reads) */
class ReadsManager
{
  public:
    fs::path sPrimaryQueryFullFileName;
    fs::path sMateQueryFullFileName;

    /* For paired reads (Illumina etc.) */
    void setReadsFileSource( const std::string& rsPrimaryQuery, const std::string& rsMateQuery )
    {
        this->sPrimaryQueryFullFileName = rsPrimaryQuery;
        this->sMateQueryFullFileName = rsMateQuery;
    } // method

    void setReadsFileSource( const std::string& rsPrimaryQuery )
    {
        this->sPrimaryQueryFullFileName = rsPrimaryQuery;
        this->sMateQueryFullFileName.clear( );
    } // method

    /* Delivers folder of primary reads */
    fs::path getReadsFolderPath( void )
    {
        return sPrimaryQueryFullFileName.parent_path( );
    } // method

    /* Delivers filename of primary reads without suffix */
    fs::path getReadsFileNameStem( void )
    {
        return sPrimaryQueryFullFileName.stem( );
    } // method

    void reset( void )
    {
        this->sPrimaryQueryFullFileName.clear( );
        this->sMateQueryFullFileName.clear( );
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
        xKindsOfOutput.emplace( "SQLite Database", "unused" );
    } // constructor

    /* Delivers date-time string for SAM-file naming */
    std::string dateTimeString( void )
    {
        time_t xTime = time( 0 ); // get time now
        tm* pLocalTimeNow = localtime( &xTime );

        char aBuffer[ 80 ];
        strftime( aBuffer, 80, "(%F)-(%H-%M-%S)", pLocalTimeNow );
        return std::string( aBuffer );
    } // method

    /* Creates full path for SAM output */
    std::string SAMFullFileName( void )
    {
        // SAM filename generation according to parameter settings.
        auto sFullFileName = ( this->pxParameterSetManager->xGlobalParameterSet.bSAMOutputInReadsFolder.value
                                   ? pxReadsManager->getReadsFolderPath( )
                                   : pxParameterSetManager->xGlobalParameterSet.xSAMOutputPath.get( ) );
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
    void doAlign( std::function<void( double dPercentageProgress )> fProgressCallBack = []( double ) {} )
    {
        std::cout << "Do align started ..." << std::endl;
        // FIXME: Read correct value for uiConcurrency
        unsigned int uiConcurency = 1;

        // FIXME: This kind of call updateGlobalParameter is not nice ...
        // Update the global parameter setting using the selected parameter-set
        xParameterSetManager.getSelected( )->updateGlobalParameter( );

        // For now, we build a computational graph for each call of doAlign
        // Possible Improvement: Cache the graph ...
        std::vector<std::shared_ptr<BasePledge>> aGraphSinks;

        // For progress computation used merely.
        std::shared_ptr<Reader> pxReader;

        const std::string sSAMFileName = xOutputManager.SAMFullFileName( );
        std::cout << "SAM output to file: " << sSAMFileName << std::endl;

        // Build a vector of aligner instances. (aGraphSinks is a vector of final pledges)
        // Each Pledge corresponds to a thread.
        // (The number of instances is decided by the degree of concurrency) .
        // The Pledge is later given to a worker for evaluation.
        if( xParameterSetManager.getSelected( )->usesPairedReads( ) )
        {
            // Wish AK: TP_PAIRED_WRITER -> PairedWriterType
            std::shared_ptr<TP_PAIRED_WRITER> pxPairedWriter;
            pxPairedWriter.reset( new PairedFileWriter( sSAMFileName, xGenomeManager.pxPackPledge->get( ) ) );

            std::vector<std::string> vA{xReadsManager.sPrimaryQueryFullFileName.string( )};
            std::vector<std::string> vB{xReadsManager.sMateQueryFullFileName.string( )};
            auto pxPairedFileReader = std::make_shared<PairedFileReader>( vA, vB );
            pxReader = pxPairedFileReader;
            auto pxPairedQueriesPledge = promiseMe( pxPairedFileReader );
            aGraphSinks = setUpCompGraphPaired( xGenomeManager.pxPackPledge, // Pack
                                                xGenomeManager.pxFMDIndexPledge, // FMD index
                                                pxPairedQueriesPledge, // (for paired reads we require two queries!)
                                                pxPairedWriter, // Output writer module(output of alignments)
                                                uiConcurency ); // Number of threads
        } // if
        else
        {
            // Wish AK: TP_WRITER -> WriterType
            std::shared_ptr<TP_WRITER> pxWriter;
            pxWriter.reset( new FileWriter( sSAMFileName, xGenomeManager.pxPackPledge->get( ) ) );

            auto pxFileReader = std::make_shared<FileReader>( xReadsManager.sPrimaryQueryFullFileName.string( ) );
            pxReader = pxFileReader;
            auto pxQueriesPledge = promiseMe( pxFileReader );
            aGraphSinks = setUpCompGraph( xGenomeManager.pxPackPledge, // Pack
                                          xGenomeManager.pxFMDIndexPledge, // FMD index
                                          pxQueriesPledge, // Queries
                                          pxWriter, // Output writer module(output of alignments)
                                          uiConcurency ); // Number of threads
        } // else

        // Compute the actual alignments.
        // Sets the progress bar after each finished alignment.
        BasePledge::simultaneousGet(
            aGraphSinks,
            [&]( ) {
                fProgressCallBack( ( (double)pxReader->getCurrPosInFile( ) * 100 ) / (double)pxReader->getFileSize( ) );
                return true;
            } // lambda
        ); // function call

        // Destroy computational graph; release memory
        aGraphSinks.clear( );
    } // method
#if 0
	void onAlign( wxCommandEvent& /*e*/ )
    {
        if( pxGenomeSelectionChoice->GetSelection( ) == wxNOT_FOUND )
        {
            wxMessageBox( "You have to select a reference genome.", "Error", wxICON_ERROR );
            return;
        } // if

        bool bInputIsFromFile = !xQueryEdit->IsEditable( );
        if( paired && !xMateEdit->IsEditable( ) != bInputIsFromFile )
        {
            wxMessageBox( "Cannot do paired alignment with one input from file and one from text.", "Error",
                          wxICON_ERROR );
            return;
        } // if


        if( bConsole )
        {
#ifdef _MSC_VER
            AllocConsole( );
            freopen( "CONOUT$", "w", stdout );
            freopen( "CONOUT$", "w", stderr );
#endif
        }
        /* Build the pack for the selected genome */
        try
        {
            auto iSelction = pxGenomeSelectionChoice->GetSelection( );
            std::string s = vRefPrefix[ iSelction ];
            auto pPackPledge = makePledge<Pack>( s );

            auto pFMDIndexPledge = makePledge<FMIndex>( vRefPrefix[ iSelction ] );
            std::shared_ptr<Reader> pReader;


            if( bInputIsFromFile )
            {
                // sets the progress bar range
                if( paired )
                {
                    std::shared_ptr<TP_PAIRED_WRITER> pOut;
                    wxFileDialog saveFileDialog( this, _( "Save Alignment to file" ), "", "",
                                                 "SAMFiles(*.sam)|*.sam|All Files (*.*)|*.*", wxFD_SAVE );
                    if( saveFileDialog.ShowModal( ) == wxID_OK )
                        pOut.reset(
                            new PairedFileWriter( std::string( saveFileDialog.GetPath( ) ), pPackPledge->get( ) ) );
                    else
                        return;
                    auto pFileReader = std::make_shared<PairedFileReader>( sQueryIn, sMateIn );
                    auto pQueriesPledge = promiseMe( pFileReader );
                    pReader = pFileReader;
                    /* Build a vector of aligner instances. (The number of instances is decided by the degress og
                     * concurrency) aGraphSinks is a vector of final pledges. Each Pledge corresponds to a thread.
                     */
                    aGraphSinks = setUpCompGraphPaired( pPackPledge, // pack
                                                        pFMDIndexPledge, // FMD index
                                                        pQueriesPledge, // (for paired reads we require two queries!)
                                                        pOut, // Output writer module(output of alignments)
                                                        uiConcurency // num threads
                    );
                } // if
                else
                {
                    // bool bPrintInfos = false;
                    std::shared_ptr<TP_WRITER> pOut;
                    wxFileDialog saveFileDialog( this, _( "Save Alignment to file" ), "", "",
                                                 "SAMFiles(*.sam)|*.sam|All Files (*.*)|*.*", wxFD_SAVE );
                    if( saveFileDialog.ShowModal( ) == wxID_OK )
                        pOut.reset( new FileWriter( std::string( saveFileDialog.GetPath( ) ), pPackPledge->get( ) ) );
                    else
                        return;

                    auto pFileReader = std::make_shared<FileReader>( sQueryIn );
                    auto pQueriesPledge = promiseMe( pFileReader );
                    pReader = pFileReader;
                    /* Build a vector of aligner instances. (The number of instances is decided by the degress og
                     * concurrency) aGraphSinks is a vector of final pledges. Each Pledge corresponds to a thread.
                     */
                    aGraphSinks = setUpCompGraph( pPackPledge, // pack
                                                  pFMDIndexPledge, // FMD index
                                                  pQueriesPledge, // queries
                                                  pOut, // Output writer module(output of alignments)
                                                  uiConcurency // num threads
                    );
                } // else
            } // if
            // else @todo re-enable input from text box
            //{
            //    auto pQuery = std::make_shared<NucSeq>( std::string( xQueryEdit->GetValue( ) ) );
            //    aQueriesPledge.push_back( std::make_shared<Pledge>( std::make_shared<NucSeq>( ) ) );
            //    aQueriesPledge.back( )->set( pQuery );
            //
            //    if( paired )
            //    {
            //        auto pQuery2 = std::make_shared<NucSeq>( std::string( xMateEdit->GetValue( ) ) );
            //        aQueriesPledge.push_back( std::make_shared<Pledge>( std::make_shared<NucSeq>( ) ) );
            //        aQueriesPledge.back( )->set( pQuery2 );
            //    } // if
            //} // else


            try
            {
                // run the alignment
                if( pWorker != nullptr )
                    pWorker->join( );
                pWorker = std::make_shared<std::thread>( [&]( ) {
                    Disable( );
                    /**
                     * compute the actual alignments.
                     * sets the progress bar after each finished alignment.
                     */
                    BasePledge::simultaneousGet( aGraphSinks,
                                                 [&]( ) {
                                                     if( bInputIsFromFile )
                                                     { // sets the progress bar
                                                         size_t uiCurrProg = ( 1000 * pReader->getCurrPosInFile( ) ) /
                                                                             pReader->getFileSize( );
                                                         pProgressBar->SetValue( uiCurrProg );
                                                     }
                                                     return true;
                                                 } // lambda
                    );
                    // release memory
                    {
                        // delete the computational graph
                        aGraphSinks.clear( );
                        auto xCommandEvent = wxCommandEvent( );
                        onRefClear( xCommandEvent );
                    }
                    {
                        // clear the main query field
                        auto xCommandEvent = wxCommandEvent( );
                        onQueryClear( xCommandEvent );
                    }
                    {
                        // clear the mate query field
                        auto xCommandEvent = wxCommandEvent( );
                        onMateClear( xCommandEvent );
                    }
                    pProgressBar->SetValue( 0 );
                    wxMessageBox( "Task completed.", "Done", wxICON_INFORMATION );
                    Enable( );
                } // lambda
                );
            }
            catch( AnnotatedException e )
            {
                std::cerr << "failed setting up aligner\n" << std::endl;
                std::cerr << e.what( ) << std::endl;
                return;
            }
            catch( ... )
            {
                std::cerr << "failed setting up aligner\n" << std::endl;
                std::cerr << "got unkown exception\n" << std::endl;
                return;
            }
#if 0
			Task* frame = new Task(
				"test",
				[&]
				(Task* t, std::string sQueryIn, unsigned int uiToFileSelect, bool bPaired, std::string sMateIn, unsigned int uiThreads,
					unsigned int uiSeedSet, unsigned int uiDistribution, unsigned int uiMean, unsigned int uiSTD, unsigned int pairedPenalty,
					unsigned int uiReportN, bool bGlobal, unsigned int uiMaxAmbiguity, unsigned int uiNumSOC, std::string sPref, bool bConsole)
				{
                
						t->print("done\n");
				},
				sQueryIn,
				uiToFileSelect,
				paired,
				sMateIn,
				uiConcurency,
				(selected == 0 ? /*fast:*/ 0 : (selected == 1 ? /*accurate:*/ 1 : /*custom:*/ uiSeedSetSelect) ),
				uiDistributionSelect,
				uiMeanDistance,
				uiSTD,
				uiUnpairedPenalty,
				uiReport,
				bGlobalAlignment,
				(selected == 0 ? /*fast:*/ 10 : (selected == 1 ? /*accurate:*/ 100 : /*custom:*/ uiMaxAmbiguity)),
				(selected == 0 ? /*fast:*/ 2 :  (selected == 1 ? /*accurate:*/  10 : /*custom:*/ uiNumSOC)),
				aRefPrefix[pxGenomeSelectionChoice->GetSelection()],
				bConsole
			);
			frame->Show();
#endif
        }
        catch( std::exception& rxEception )
        {
            std::cout << "Caught: " << rxEception.what( ) << std::endl;
        }
    }
#endif
}; // class
