/////////////////////////////////////////////////////////////////////////////
// Name:        minimal.cpp
// Purpose:     Minimal wxWidgets sample
// Author:      Julian Smart
// Modified by:
// Created:     04/01/98
// Copyright:   (c) Julian Smart
// Licence:     wxWindows licence
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <fstream> // STL
#include <functional> // STL
#include <iostream> //STL
#include <string> // STL
#include <thread> // STL
#include <vector> // STL


// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWidgets headers)
#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "mwx.h"
/* Additional wxWidgets includes that are not part of wx.h already */
#include <wx/artprov.h>
#include <wx/bmpcbox.h>
#include <wx/event.h> // Bind
#include <wx/filepicker.h>
#include <wx/html/helpctrl.h>
#include <wx/notebook.h>
#include <wx/tglbtn.h>
#include <wx/wizard.h>
#include <wx/wrapsizer.h>

/* Due to MSVC libMa should be included after the wxWidgets includes */
#include "util/default_parameters.h"
#include "util/export.h"
using namespace libMA;

#include "execution-context.h"
#include "parameter.h"

// Global MA execution context
ExecutionContext xExecutionContext;

/* Include XPM bitmpas */
#include "IconMA.xpm" // defines xIconMA_XPM

/* Global execution context for ID management.
 * Defined in mwx.h
 * TODO: Remove me
 */
int mxwExecutionContext::iHighestID = wxID_HIGHEST + 30;

/* Define all event types */
wxDEFINE_EVENT( wxEVT_WORKER_MESSAGE, wxCommandEvent );
wxDEFINE_EVENT( wxEVT_WORKER_COMPLETED, wxCommandEvent );
wxDEFINE_EVENT( wxEVT_WORKER_UPDATE_GAUGE, wxCommandEvent );
wxDEFINE_EVENT( wxEVT_WORKER_ALIGNMENT_ERROR, wxCommandEvent );

/* Settings Dialog for aligner parameter management */
class mwxSettingsDialog : public mwxOK_Cancel_Dialog
{
    int iValue;

  public:
    /* Constructor */
    mwxSettingsDialog( wxWindow* pxHostWindow, // Host window of box context (responsible for destruction)
                       Presetting& rxParameterSet // Parameter used for the dialog
                       )
        : mwxOK_Cancel_Dialog(
              pxHostWindow, // host ( responsible for destruction)
              "Parameter Settings", // Title of Dialog
              [&]( wxWindow* pxHostWindow ) { // Content of Dialog
                  this->iValue = 0;
                  auto* pxNotebook = new mwxPropertyNotebook( pxHostWindow );

                  auto* pxScrolledStatixBoxesContext = pxNotebook->addPage( "Custom" );
                  ( new mwxPropertyPanel( pxScrolledStatixBoxesContext->addStaticBox( )->getConnector( ) ) )
                      ->append( "Maximal Ambiguity", &this->iValue )
                      .append( "Number of SoC", &this->iValue )
                      .append( rxParameterSet.xSeedingTechnique )
                      .append( rxParameterSet.xExampleCheckBox );

                  pxScrolledStatixBoxesContext = pxNotebook->addPage( "General" );
                  ( new mwxPropertyPanel( pxScrolledStatixBoxesContext->addStaticBox( )->getConnector( ) ) )
                      ->append( "Report", &iValue )
                      .append( rxParameterSet.xUsePairedReads )
                      .append( rxParameterSet.xExamplePath )
                      .appendCheckBox( "Global ALignment", &iValue )
                      .append( "Write to file if", {"input is from file", "always", "never"}, &iValue )
                      .append( "Concurrency", &iValue )
                      .append( "Max NM area", &iValue )
                      .appendCheckBox( "Show stdout/stderr", &iValue );

                  pxScrolledStatixBoxesContext = pxNotebook->addPage( "Paired Alignments" );
                  ( new mwxPropertyPanel( pxScrolledStatixBoxesContext->addStaticBox( )->getConnector( ) ) )
                      ->append( "Used distribution", {"normal", "uniform"}, &iValue )
                      .append( "Mean distance", &iValue )
                      .append( "Standard deviation", &iValue )
                      .append( "Unpaired penalty", &iValue );

                  pxScrolledStatixBoxesContext = pxNotebook->addPage( "DP Scoring" );
                  ( new mwxPropertyPanel( pxScrolledStatixBoxesContext->addStaticBox( )->getConnector( ) ) )
                      ->append( rxParameterSet.xMatch )
                      .append( rxParameterSet.xMisMatch )
                      .append( "Gap open penalty", &iValue )
                      .append( "Gap extend penalty", &iValue );

                  return pxNotebook;
              }, // lambda
              wxDefaultPosition )
    {} // constructor

    /* Destructor */
    virtual ~mwxSettingsDialog( )
    {
        std::cout << "Destruct mwxSettingsDialog" << std::endl;
    } // destructor
}; // class


/* Dialog for entering the paired read parameter */
class mwxPairedSettingsDialog : public mwxOK_Cancel_Dialog
{
    int iValue = 0;

  public:
    /* Constructor */
    mwxPairedSettingsDialog( wxWindow* pxHostWindow, // Host window of box context (responsible for destruction)
                             Presetting& rxParameterSet // Parameter used for the dialog
                             )
        : mwxOK_Cancel_Dialog(
              pxHostWindow, // host ( responsible for destruction)
              "Paired Reads Settings", // Title of Dialog
              [&]( wxWindow* pxHostWindow ) { // Content of Dialog
                  auto* pxPropertyPanel = new mwxPropertyPanel( pxHostWindow, NULL );
                  pxPropertyPanel->append( "Mean distance", &iValue )
                      .append( "Used distribution", {"normal", "uniform"}, &iValue )
                      .append( "Standard deviation", &iValue )
                      .append( "Unpaired penalty", &iValue );

                  return pxPropertyPanel;
              }, // lambda

              wxDefaultPosition,
              true ) // fit the size of the dialog to the size of its content
    {} // constructor

    /* Destructor */
    virtual ~mwxPairedSettingsDialog( )
    {
        std::cout << "Destruct mwxPairedSettingsDialog" << std::endl;
    } // destructor
}; // class


/* Dialog for entering the paired read parameter */
class mwxSAMSettingsDialog : public mwxOK_Cancel_Dialog
{
  public:
    /* Constructor */
    mwxSAMSettingsDialog( wxWindow* pxHostWindow, // Host window of box context (responsible for destruction)
                          GeneralParameter& rxGlobalParameterSet // Parameter used for the dialog
                          )
        : mwxOK_Cancel_Dialog(
              pxHostWindow, // host ( responsible for destruction)
              "SAM Settings", // Title of Dialog
              [&]( wxWindow* pxHostWindow ) { // Content of Dialog
                  auto* pxPropertyPanel = new mwxPropertyPanel( pxHostWindow, NULL );
                  pxPropertyPanel->append( rxGlobalParameterSet.bSAMOutputInReadsFolder )
                      .append( rxGlobalParameterSet.xSAMOutputPath );

                  return pxPropertyPanel;
              }, // lambda

              wxDefaultPosition,
              true ) // fit the size of the dialog to the size of its content
    {} // constructor
}; // class


WX_DEFINE_ARRAY_PTR( wxWizardPageSimple*, WizardPages );

class FMIndexCreationWizard : public wxWizard
{
    // Three buttons of wizards can be found via FindWindowById.
    // Ids are: wxID_CANCEL, wxID_FORWARD, wxID_BACKWARD
  protected:
    wxFilePickerCtrl* pxFilePickerIndexFasta;
    wxTextCtrl* pxTextCtrlIndexName;
    wxDirPickerCtrl* pxDirPickerIndexLocation;
    wxTextCtrl* pxTextCtrlConsole;
    wxWizardPageSimple* pxWizardPage2;

    void hideButtonByID( wxWindowID wxId = wxID_CANCEL )
    {
        wxButton* pxButton = static_cast<wxButton*>( this->FindWindowById( wxId ) );
        if( pxButton )
            pxButton->Hide( );
    } // method

    void disableButtonById( wxWindowID wxId = wxID_CANCEL )
    {
        wxButton* pxButton = static_cast<wxButton*>( this->FindWindowById( wxId ) );
        if( pxButton )
            pxButton->Disable( );
    } // method

    void enableButtonById( wxWindowID wxId = wxID_CANCEL )
    {
        wxButton* pxButton = static_cast<wxButton*>( this->FindWindowById( wxId ) );
        if( pxButton )
            pxButton->Enable( );
    } // method

    /* Handler for printing from worker */
    void onWorkerPrint( wxCommandEvent& rxEvent )
    {
        pxTextCtrlConsole->AppendText( rxEvent.GetString( ) );
    } // method

    /* Handler for completion by worker */
    void onWorkerCompleted( wxCommandEvent& rxEvent )
    {
        // enable forward button
        enableButtonById( wxID_FORWARD );
    } // method

    /* Queues event of the worker comprising text */
    void queueStringMessage( const wxEventTypeTag<wxCommandEvent>& rxEventType, const std::string& sMessage )
    {
        wxCommandEvent* pxEvent = new wxCommandEvent( rxEventType, wxID_ANY );
        pxEvent->SetString( wxString( sMessage ) );
        this->QueueEvent( pxEvent ); // handler takes ownership
    } // method

    std::unique_ptr<std::thread> pWorker;

    /* Handler called, if 'Next'-button pushed.
     * We need a worker, so that the event-loop can continue ...
     */
    void onPageShown( void )
    {
        // If we see the second page, we automatically start the index creation
        if( GetCurrentPage( ) == pxWizardPage2 )
        {
            // Disable Finish Button and hide Cancel Button
            disableButtonById( wxID_FORWARD );
            hideButtonByID( wxID_CANCEL );

            // FM-Index creation happens in a separated thread
            pWorker = std::make_unique<std::thread>( [this]( ) {
                auto sGenomePrefix = fs::path( sGenomeFolderPath( ) ) /= sFastaFilePath( ).stem( );
                queueStringMessage( wxEVT_WORKER_MESSAGE,
                                    std::string( "Genome prefix:\n" ) + sGenomePrefix.string( ) + "\n" );

                // Create and store pack
                queueStringMessage( wxEVT_WORKER_MESSAGE, "Create pack and write pack to file system.\n" );
                std::shared_ptr<Pack> pxPack( new Pack( ) );
                pxPack->vAppendFASTA( sFastaFilePath( ).string( ) );
                pxPack->vStoreCollection( sGenomePrefix.string( ) );

                // Create and store FMD index
                queueStringMessage( wxEVT_WORKER_MESSAGE,
                                    "Compute FMD-index...\nImportant note: This may take long time.\n" );
                FMIndex xFMDIndex( pxPack );
                queueStringMessage( wxEVT_WORKER_MESSAGE, "Write FMD-Index to file system.\n" );
                xFMDIndex.vStoreFMIndex( sGenomePrefix.string( ).c_str( ) );

                // FIXME: Write genome name to description file.
                queueStringMessage( wxEVT_WORKER_MESSAGE, "All done!\n" );

                // Enable 'Finish' button
                queueStringMessage( wxEVT_WORKER_COMPLETED, "" );
            } // lambda
            );
        } // if
    } // method

    template <class VALIDATION_TARGET_CLASS> class WizardPage1Validator : public wxValidator
    {
      public:
        const VALIDATION_TARGET_CLASS* pxValidationTarget;

        /* Constructor */
        WizardPage1Validator( const VALIDATION_TARGET_CLASS* pxValidationTarget )
            : pxValidationTarget( pxValidationTarget )
        {}

        virtual wxObject* Clone( void ) const
        {
            return new WizardPage1Validator( this->pxValidationTarget );
        } // method

        virtual bool Validate( wxWindow* parent )
        {
            return pxValidationTarget->validate( );
        } // method

        /* Called to transfer data to the window */
        virtual bool TransferToWindow( )
        {
            return true;
        } // method

        /* Called to transfer data from the window */
        virtual bool TransferFromWindow( )
        {
            return true;
        } // method
    }; // class

  public:
    WizardPages xWizardPages;

    /* Returns FASTA file-path of dialog */
    fs::path sFastaFilePath( void ) const
    {
        return ( fs::path( std::string( pxFilePickerIndexFasta->GetPath( ).c_str( ) ) ) );
    } // method

    /* Returns folder for genome storage */
    fs::path sGenomeFolderPath( void ) const
    {
        return ( fs::path( std::string( pxDirPickerIndexLocation->GetPath( ).c_str( ) ) ) );
    } // method

    /* Input validator (only connected to the first page) */
    bool validate( void ) const
    {
        if( !fs::exists( sFastaFilePath( ) ) )
        {
            wxMessageBox( "Please select a FASTA file for index generation.", "Missing Input", wxICON_ERROR );
            return false;
        } // if
        if( this->pxTextCtrlIndexName->GetValue( ).empty( ) )
        {
            wxMessageBox( "Please enter an appropriate name for your genome.", "Missing Input", wxICON_ERROR );
            return false;
        } // if
        if( sGenomeFolderPath( ).empty( ) )
        {
            wxMessageBox( "Please select a folder for genome storage.", "Missing Input", wxICON_ERROR );
            return false;
        } // if

        return true;
    } // method

    /* Constructor */
    FMIndexCreationWizard( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxEmptyString,
                           const wxBitmap& bitmap = wxNullBitmap, const wxPoint& pos = wxDefaultPosition,
                           long style = wxDEFAULT_DIALOG_STYLE )
        : wxWizard( parent, id, title, bitmap, pos, style )
    {
        SetBackgroundColour( wxColor( 255, 255, 255 ) );
        SetSizeHints( wxDefaultSize, wxDefaultSize );
        CenterOnParent( );

        // Page 1 of Wizard
        wxWizardPageSimple* pxWizardPage1 = new wxWizardPageSimple( this );

        xWizardPages.Add( pxWizardPage1 );

        auto pxBoxSizerPage1 = new wxBoxSizer( wxVERTICAL );
        pxBoxSizerPage1->Add( new wxStaticText( pxWizardPage1, wxID_ANY,
                                                wxT( "1. Select a Fasta file for index creation" ), wxDefaultPosition,
                                                wxDefaultSize, 0 ),
                              0, wxALL | wxEXPAND, 5 );
        pxBoxSizerPage1->Add( pxFilePickerIndexFasta = new wxFilePickerCtrl(
                                  pxWizardPage1, wxID_ANY, wxEmptyString, wxT( "Select a file" ), wxT( "*.*" ),
                                  wxDefaultPosition, wxDefaultSize, wxFLP_DEFAULT_STYLE ),
                              0, wxALL | wxEXPAND, 5 );

        pxBoxSizerPage1->Add( new wxStaticText( pxWizardPage1, wxID_ANY,
                                                wxT( "2. Give your index a name (e.g. Human-hg38)" ), wxDefaultPosition,
                                                wxDefaultSize, 0 ),
                              0, wxALL, 5 );
        pxBoxSizerPage1->Add( pxTextCtrlIndexName = new wxTextCtrl( pxWizardPage1, wxID_ANY, wxEmptyString,
                                                                    wxDefaultPosition, wxDefaultSize, 0 ),
                              0, wxALL | wxEXPAND, 5 );
        pxBoxSizerPage1->Add( new wxStaticText( pxWizardPage1, wxID_ANY, wxT( "3. Select a folder for Index storage" ),
                                                wxDefaultPosition, wxDefaultSize, 0 ),
                              0, wxALL, 5 );
        pxBoxSizerPage1->Add( pxDirPickerIndexLocation =
                                  new wxDirPickerCtrl( pxWizardPage1, wxID_ANY, wxEmptyString, wxT( "Select a folder" ),
                                                       wxDefaultPosition, wxDefaultSize, wxDIRP_DEFAULT_STYLE ),
                              0, wxALL | wxEXPAND, 5 );


        pxBoxSizerPage1->Add( new wxStaticText( pxWizardPage1, wxID_ANY,
                                                wxT( "4. Click \"Next\" for starting the index creation." ),
                                                wxDefaultPosition, wxDefaultSize, 0 ),
                              0, wxALL, 5 );

        pxWizardPage1->SetSizer( pxBoxSizerPage1 );
        pxWizardPage1->Layout( );
        // pxWizardPage1->Center(wxBOTH);
        pxBoxSizerPage1->Fit( pxWizardPage1 );

        // Next button of page 1 is connected to validator for input verification
        pxWizardPage1->SetValidator( WizardPage1Validator<FMIndexCreationWizard>( this ) );

        // Page 2 of Wizard
        pxWizardPage2 = new wxWizardPageSimple( this );
        xWizardPages.Add( pxWizardPage2 );

        wxBoxSizer* pxBoxSizerPage2;
        pxBoxSizerPage2 = new wxBoxSizer( wxVERTICAL );
        pxBoxSizerPage2->Add( new wxStaticText( pxWizardPage2, wxID_ANY,
                                                wxT( "Index creation is now underway.\n"
                                                     "Depending on your hardware and the size of\n"
                                                     "the genome, this can take a long time." ),
                                                wxDefaultPosition, wxDefaultSize, 0 ),
                              0, wxALL, 5 );
        // Read-only multi-line window
        pxBoxSizerPage2->Add( pxTextCtrlConsole =
                                  new wxTextCtrl( pxWizardPage2, wxID_ANY, wxEmptyString, wxDefaultPosition,
                                                  wxDefaultSize, wxTE_MULTILINE | wxTE_READONLY ),
                              1, wxALL | wxEXPAND, 5 );

        pxWizardPage2->SetSizer( pxBoxSizerPage2 );
        pxWizardPage2->Layout( );
        // pxWizardPage2->Center(wxBOTH);
        pxBoxSizerPage2->Fit( pxWizardPage2 );

        for( unsigned int i = 1; i < xWizardPages.GetCount( ); i++ )
        {
            xWizardPages.Item( i )->SetPrev( xWizardPages.Item( i - 1 ) );
            xWizardPages.Item( i - 1 )->SetNext( xWizardPages.Item( i ) );
        } // for

        wxEvtHandler::Bind( wxEVT_WIZARD_PAGE_SHOWN, [&]( wxCommandEvent& rxEvent ) { onPageShown( ); } );

        // Bind worker events
        Bind( wxEVT_WORKER_MESSAGE, std::bind( &FMIndexCreationWizard::onWorkerPrint, this, std::placeholders::_1 ) );
        Bind( wxEVT_WORKER_COMPLETED,
              std::bind( &FMIndexCreationWizard::onWorkerCompleted, this, std::placeholders::_1 ) );

        // Switch of the 'back' button
        hideButtonByID( wxID_BACKWARD );

        Layout( );
        Centre( wxBOTH );
    } // constructor

    /* Destructor */
    ~FMIndexCreationWizard( void )
    {
        std::cout << "Destructor FMIndexCreationWizard" << std::endl;
        // Before destroying the worker we must join the worker! (must of thread-class)
        if( pWorker )
            pWorker->join( );
        this->xWizardPages.Clear( );
    }
}; // class


// Define a new application type, each program should derive a class from wxApp
class MyApp : public wxApp
{
  public:
    // override base class virtuals
    // ----------------------------

    // this one is called on application startup and is a good place for the app
    // initialization (doing it here and not in the ctor allows to have an error
    // return: if OnInit() returns false, the application terminates)
    virtual bool OnInit( ); // In 3.1: virtual bool OnInit( ) wxOVERRIDE;
};

// Create a new application object: this macro will allow wxWidgets to create
// the application object during program execution (it's better than using a
// static object for many reasons) and also implements the accessor function
// wxGetApp() which will return the reference of the right type (i.e. MyApp and
// not wxApp)

wxIMPLEMENT_APP( MyApp );

/* Some global parameter of MA-Gui */
unsigned int uiConcurency = std::thread::hardware_concurrency( );
unsigned int uiToFileSelect = 0;
bool bConsole = false;

// Define a new frame type: this is going to be our main frame
// Example for multi-threaded wxWidgtes Application: https://forums.wxwidgets.org/viewtopic.php?t=8392
// https://wiki.wxwidgets.org/Inter-Thread_and_Inter-Process_communication


class AlignFrame : public wxDialog
{
  public:
    wxTextCtrl* pxTextCtrl;
    wxGauge* pxGauge;
    wxButton* pxStopOKButton;
    double dPreviousPercentage = -1;
    std::unique_ptr<std::thread> pWorker;

    /* Callback called by worker.
     */
    void onCallBack( double dPercentage )
    {
        dPercentage = round( dPercentage );
        if( dPercentage > dPreviousPercentage )
        {
            dPreviousPercentage = dPercentage;
            auto* pxEvent = new wxCommandEvent( wxEVT_WORKER_UPDATE_GAUGE, wxID_ANY );
            pxEvent->SetInt( static_cast<int>( dPercentage ) );
            this->QueueEvent( pxEvent ); // owner-ship transfered to handler
        }
    } // method

    /* Handler for printing from worker */
    void onWorkerGaugeUpdate( wxCommandEvent& rxEvent )
    {
        pxGauge->SetValue( static_cast<int>( rxEvent.GetInt( ) ) );
    } // method

    /* Handler for printing from worker */
    void onWorkerPrint( wxCommandEvent& rxEvent )
    {
        pxTextCtrl->AppendText( rxEvent.GetString( ) );
    } // method

    /* Handler for completion by worker */
    void onWorkerCompleted( wxCommandEvent& rxEvent )
    {
        pxStopOKButton->Enable( );
    } // method

    /* Handler for display of error messages */
    void onWorkerError( wxCommandEvent& rxEvent )
    {
        wxMessageBox( "Aligner failed due to exception:\n" + rxEvent.GetString( ), "Error", wxICON_ERROR );
        // Handler enables the OK button
        pxTextCtrl->AppendText( "Alignment failed\n" );
        pxStopOKButton->Enable( );
    } // method

    /* Queues event of the worker comprising text */
    void queueStringMessage( const wxEventTypeTag<wxCommandEvent>& rxEventType, const std::string& sMessage )
    {
        wxCommandEvent* pxEvent = new wxCommandEvent( rxEventType, wxID_ANY );
        pxEvent->SetString( wxString( sMessage ) );
        this->QueueEvent( pxEvent );
    } // method

    // The wxWidgets ShowModal is not thread-safe.
    int ShowModal( void )
    {
        // create worker thread, which communicates via events
        pWorker = std::make_unique<std::thread>( [&]( ) {
            queueStringMessage( wxEVT_WORKER_MESSAGE, "Do alignment\n" );
            // If sErrorMessage is empty, the alignment is fine succeeded
            std::string sErrorMessage;
            bool bExexutionSucceeded = false;

            try
            {
                xExecutionContext.doAlign( std::bind( &AlignFrame::onCallBack, this, std::placeholders::_1 ) );
                bExexutionSucceeded = true;
            } // try
            catch( std::exception& xException )
            {
                queueStringMessage( wxEVT_WORKER_ALIGNMENT_ERROR, std::string( xException.what( ) ) + "\n" );
            } // catch
            catch( ... )
            {
                queueStringMessage( wxEVT_WORKER_ALIGNMENT_ERROR, "Aligner failed due to unknown reason.\n" );
            } // catch ellipsis

            if( bExexutionSucceeded )
            { // Enable the OK button only if the alignment succeeded
                queueStringMessage( wxEVT_WORKER_MESSAGE, "Task completed\n" );
                queueStringMessage( wxEVT_WORKER_COMPLETED, "" ); // triggers activation of OK button
            } // if
        } // lambda
        );

        // Continue in local event loop until button is clicked.
        // This loop checks if the Stop button has been pressed.
        int iModalReturnValue = wxDialog::ShowModal( );

        return iModalReturnValue;
    } // method

    /* Task Frame; created by clicking the start button */
    AlignFrame( wxWindow* parent, const wxString& title )
        : wxDialog( parent, wxID_ANY, title, wxDefaultPosition, wxDefaultSize, wxDEFAULT_FRAME_STYLE )
    {
        // Set frame icon and background color
        SetIcons( wxIcon( xIconMA_XPM ) );
        SetBackgroundColour( wxColor( 255, 255, 255 ) );
        SetSizeHints( wxDefaultSize, wxDefaultSize );
        CenterOnParent( );

        // Add textCTrl, gauge (progress bar), button
        auto* pxSizer = new wxBoxSizer( wxVERTICAL );
        pxSizer->Add( new wxStaticText( this, wxID_ANY, wxT( "Aligner Log Output" ) ), 0, wxALL, 5 );
        pxSizer->Add( pxTextCtrl = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize,
                                                   wxTE_MULTILINE | wxTE_READONLY ),
                      1, wxBOTTOM | wxEXPAND | wxLEFT | wxRIGHT, 5 );
        pxSizer->Add( new wxStaticText( this, wxID_ANY, wxT( "Progress" ) ), 0, wxALL, 5 );

        this->pxGauge = new wxGauge( this, wxID_ANY, 100, wxDefaultPosition, wxDefaultSize, wxGA_HORIZONTAL );
        this->pxGauge->SetValue( 0 );
        pxSizer->Add( this->pxGauge, 0, wxBOTTOM | wxEXPAND | wxLEFT | wxRIGHT, 5 );

        pxSizer->Add( pxStopOKButton = new wxButton( this, wxID_OK, wxT( "OK" ) ), 0, wxALIGN_CENTER | wxALL, 5 );

        SetSizer( pxSizer );
        Layout( );
        pxStopOKButton->Disable( );
        Centre( wxBOTH );

        // Bind worker events
        Bind( wxEVT_WORKER_MESSAGE, std::bind( &AlignFrame::onWorkerPrint, this, std::placeholders::_1 ) );
        Bind( wxEVT_WORKER_COMPLETED, std::bind( &AlignFrame::onWorkerCompleted, this, std::placeholders::_1 ) );
        Bind( wxEVT_WORKER_UPDATE_GAUGE, std::bind( &AlignFrame::onWorkerGaugeUpdate, this, std::placeholders::_1 ) );
        Bind( wxEVT_WORKER_ALIGNMENT_ERROR, std::bind( &AlignFrame::onWorkerError, this, std::placeholders::_1 ) );
    } // constructor

    ~AlignFrame( void )
    {
        if( pWorker )
            pWorker->join( );
        std::cout << "Worker joined... " << std::endl;
    } // destructor
}; // wxFrame


class MA_MainFrame : public wxFrame
{
  private:
    mwxMenuBar* xMenuBar; // Deallocation by host-frame

    bool inPairedMode( void )
    {
        return xExecutionContext.xParameterSetManager.getSelected( )->usesPairedReads( );
    } // method

    /* Handler for menu Item 'About' */
    void onAbout( wxCommandEvent& WXUNUSED( event ) )
    {
        wxMessageBox( wxT( "MA - The Modular Aligner\nCreated by Markus Schmidt and Arne "
                           "Kutzner\nVersion: 1.0" ),
                      wxT( "About" ), wxICON_INFORMATION );
    } // method

    /* Handler for the gear button */
    template <class SETTINGS_DIALOG_CLASS> void doSettingsDialog( wxCommandEvent& WXUNUSED( event ) )
    {
        Presetting* pSelectedParamSet = xExecutionContext.xParameterSetManager.getSelected( );

        // Create a copy of the selected parameter-set for editing
        Presetting xParameterSet( *pSelectedParamSet, "" );
        auto xSettingsDialog = new SETTINGS_DIALOG_CLASS( this, xParameterSet );
        auto iReturnedChoice = xSettingsDialog->ShowModal( );
        if( iReturnedChoice == wxID_OK )
        {
            // save the changed parameter into the selected parameter-set
            pSelectedParamSet->mirror( xParameterSet );
        } // if

        xSettingsDialog->Destroy( );

        // The parameter-set might have changed to paired ...
        this->updateLayout( );
    } // method

    /* Handler for gear button  of aligner settings */
    void onSettingsGearButton( wxCommandEvent& event )
    {
        this->doSettingsDialog<mwxSettingsDialog>( event );
    } // method

    /* Handler for gear button of outputs settings */
    void onOutputGearButton( wxCommandEvent& WXUNUSED( event ) )
    {
        GeneralParameter xGlobalParameterSet( xExecutionContext.xParameterSetManager.xGlobalParameterSet );
        mwxSAMSettingsDialog xSettingsDialog( nullptr, xGlobalParameterSet );
        if( xSettingsDialog.ShowModal( ) == wxID_OK )
        {
            xExecutionContext.xParameterSetManager.xGlobalParameterSet.mirror( xGlobalParameterSet );
        } // if
    } // method

    void onPairedGearButton( wxCommandEvent& event )
    {
        this->doSettingsDialog<mwxPairedSettingsDialog>( event );
    } // method

    // Paired reads related elements
    wxStaticText* xMatesStaticText;
    wxTextCtrl* xMatesTextCtrl;
    mwxBoxSizer* xMatesControlElements;

    /* Updates the layout according to the selected parameter-set */
    void updateLayout( void )
    {
        bool bInPairedMode = this->inPairedMode( );
        // Show or hide elements related to paired reads
        this->xMatesStaticText->Show( bInPairedMode );
        this->xMatesTextCtrl->Show( bInPairedMode );
        this->xMatesControlElements->Show( bInPairedMode );
        this->Layout( );
    } // method

    /* Handler for parameter selection combo-box */
    void onParameterComboBox( wxCommandEvent&( event ) )
    {
        wxComboBox* pxComboBox = dynamic_cast<wxComboBox*>( event.GetEventObject( ) );
        std::cout << pxComboBox->GetSelection( ) << std::endl; // gives number of selection
        std::string sSelected( pxComboBox->GetStringSelection( ) );

        xExecutionContext.xParameterSetManager.setSelected( sSelected );
        this->updateLayout( );
    } // method

    wxTextCtrl* xGenomeNameTextCtrl;

    /* Handler for genome selection button */
    void onGenomeSelection( const std::string& sFileName )
    {
        xExecutionContext.xGenomeManager.setGenome( sFileName );
        xGenomeNameTextCtrl->SetValue( xExecutionContext.xGenomeManager.getGenomeName( ) );
    } // method

    wxTextCtrl* xQueryTextCtrl;

    /* Handler for the start button */
    void onStart( wxCommandEvent& WXUNUSED( event ) )
    {
        // First ask execution environment if all data are available
        // auto xAlignFrame = new AlignFrame( "Alignment" );

        AlignFrame xAlignFrame( this, "Alignment" );
        xAlignFrame.InitDialog( );
        if( xAlignFrame.ShowModal( ) == wxID_OK )
        {
        }
    } // method

    void onCreateIndexWizard( wxCommandEvent& WXUNUSED( event ) )
    {
        FMIndexCreationWizard xWizard( this, wxID_ANY, "FM-Index Generation", wxBITMAP_PNG_FROM_DATA( WizardLabel ),
                                       wxDefaultPosition, wxDEFAULT_DIALOG_STYLE );

        if( xWizard.RunWizard( xWizard.xWizardPages.Item( 0 ) ) )
        {
        } // if
    } // method

    /* Handler query select/clear button pair */
    void onQuerySelection( const std::string& sFileName, wxTextCtrl* xTargetTextCtrl )
    {
        // FIXME: Add code for paired inputs
        xExecutionContext.xReadsManager.sPrimaryQueryFullFileName = sFileName;
        if( sFileName.empty( ) )
            xTargetTextCtrl->SetValue( "Type your query here or select a FASTA-file ..." );
        else
            xTargetTextCtrl->SetValue( sFileName );

        xTargetTextCtrl->SetEditable( sFileName.empty( ) );
    } // method

  public:
    MA_MainFrame( const wxString& title )
        : wxFrame( NULL,
                   wxID_ANY,
                   title,
                   wxDefaultPosition,
                   wxDefaultSize,
                   wxDEFAULT_FRAME_STYLE & ~( wxMAXIMIZE_BOX | wxFRAME_FLOAT_ON_PARENT | wxSTAY_ON_TOP ) )
    {
        wxImage::AddHandler( new wxPNGHandler );

        // set icon, size and background color
        this->SetIcons( wxIcon( xIconMA_XPM ) );
        this->SetSize( wxSize( 800 + wxSYS_FRAMESIZE_X, 800 + wxSYS_FRAMESIZE_Y + 11 ) );
        this->SetBackgroundColour( wxColor( 255, 255, 255 ) );

        // Create and initialize menu-bar
        // (deallocation automatically by wxWidgets by frame object)
        xMenuBar = new mwxMenuBar;
        this->xMenuBar->push_back( "&File" )
            .menu( ) // File Menu
            .append( "O&pen\tAlt-O", // Open
                     [&]( wxCommandEvent& ) { this->Close( true ); } )
            .append( "E&xit\tAlt-X", // Exit
                     [&]( wxCommandEvent& ) { this->Close( true ); } );

        this->xMenuBar->push_back( "&Genome" )
            .menu( ) // Genome Menu
            .append( "C&reate Index\tF2", // Open
                     std::bind( &MA_MainFrame::onCreateIndexWizard, this, std::placeholders::_1 ) );

        this->xMenuBar->push_back( "&Help" )
            .menu( ) // Help Menu
            .append( "&About\tF1", std::bind( &MA_MainFrame::onAbout, this, std::placeholders::_1 ) );

        // Set menu bar for current frame. Thereby the frame becomes responsible for deallocation of the menu bar.
        SetMenuBar( xMenuBar );

        auto* pxFrameSizer = new mwxBoxSizer( mwxConnector( this, nullptr ), wxVERTICAL );
        pxFrameSizer //
            ->addBoxSizer( // horizontal BoxSizer (top of frame)
                wxHORIZONTAL,
                wxSizerFlags( 0 ).Expand( ).Border( wxALL, 10 ),
                [this]( mwxBoxSizer& pxBoxSizer ) // boxSizer content
                {
                    // MA LOGO
                    pxBoxSizer.Add( new wxStaticBitmap( pxBoxSizer.pxConnector.pxWindow, wxID_ANY,
                                                        wxBITMAP_PNG_FROM_DATA( MALogo ), wxDefaultPosition,
                                                        wxDefaultSize, 0 ),
                                    0, wxALL, 5 );

                    // Aligner settings
                    pxBoxSizer.addStaticBoxSizer( // vertical BoxSizer
                        "Aligner Settings",
                        wxVERTICAL,
                        wxSizerFlags( 0 ).Expand( ).Border( wxRIGHT, 5 ).Align( wxCENTER ),
                        [this]( mwxStaticBoxSizer& pxBoxSizer ) // BoxSizer content
                        {
                            // Parameter selection combo
                            pxBoxSizer.Add(
                                new mwxMapDrivenComboBox( //
                                    pxBoxSizer.xConnector.pxWindow,
                                    xExecutionContext.xParameterSetManager.xParametersSets,
                                    std::bind( &MA_MainFrame::onParameterComboBox, this, std::placeholders::_1 ) ),
                                0, wxTOP | wxBOTTOM | wxEXPAND, 5 );
                            pxBoxSizer.Add( 0, 0, 0, wxBOTTOM, 5 ); // Vertical spacer

                            // Gear button for aligner parameter management
                            pxBoxSizer //
                                .Add( new mwxBitmapButton //
                                      ( pxBoxSizer.xConnector.pxWindow,
                                        wxBITMAP_PNG_FROM_DATA( GearButtonNormal ),
                                        wxBITMAP_PNG_FROM_DATA( GearButtonSelected ),
                                        std::bind( &MA_MainFrame::onSettingsGearButton, this, std::placeholders::_1 ) ),
                                      0, wxALIGN_CENTER | wxTOP | wxRIGHT | wxLEFT, 5 );
                        } ); // addBoxSizer

                    // Output Selection
                    pxBoxSizer.addStaticBoxSizer( // vertical BoxSizer
                        "Output Destination",
                        wxVERTICAL,
                        wxSizerFlags( 0 ).Expand( ).Border( wxRIGHT, 5 ),
                        [this]( mwxStaticBoxSizer& pxBoxSizer ) // BoxSizer content
                        {
                            pxBoxSizer.Add( new mwxMapDrivenComboBox( //
                                                pxBoxSizer.xConnector.pxWindow,
                                                xExecutionContext.xOutputManager.xKindsOfOutput ),
                                            0, wxTOP | wxBOTTOM | wxEXPAND, 5 );
                            pxBoxSizer.Add( 0, 0, 0, wxBOTTOM, 5 ); // Vertical space

                            // Gear button for output parameter management
                            pxBoxSizer //
                                .Add( new mwxBitmapButton //
                                      ( pxBoxSizer.xConnector.pxWindow,
                                        wxBITMAP_PNG_FROM_DATA( GearButtonNormal ),
                                        wxBITMAP_PNG_FROM_DATA( GearButtonSelected ),
                                        std::bind( &MA_MainFrame::onOutputGearButton, this, std::placeholders::_1 ) ),
                                      0, wxALIGN_CENTER | wxTOP | wxRIGHT | wxLEFT, 5 );
                        } ); // addBoxSizer

                    // Quickstart-Text
                    pxBoxSizer.addStaticBoxSizer(
                        "Quickstart", wxVERTICAL,
                        wxSizerFlags( 2 ).Expand( ), //.Border( wxALL, 5 ),
                        [this]( mwxStaticBoxSizer& xStaticBoxSizer ) // BoxSizer content
                        {
                            xStaticBoxSizer.Add(
                                new wxTextCtrl(
                                    xStaticBoxSizer.xConnector.pxWindow, wxID_ANY,
                                    wxT( "1. Select the sequencing technology of your reads. (Illumina etc.)\n"
                                         "2. Choose your output format" ),
                                    wxDefaultPosition, wxDefaultSize,
                                    wxTE_AUTO_URL | wxTE_RICH | wxHSCROLL | wxTE_MULTILINE | wxTE_READONLY |
                                        wxNO_BORDER ),
                                1, wxEXPAND, 5 );
                        } );
                } ) // add horizontal BoxSizer

            ->addBoxSizer(
                wxHORIZONTAL,
                wxSizerFlags( 0 ).Expand( ).Border( wxLEFT | wxRIGHT, 10 ),
                [this]( mwxBoxSizer& pxBoxSizer ) // boxSizer content
                {
                    // Text plus combo
                    pxBoxSizer.addBoxSizer( // vertical BoxSizer
                        wxVERTICAL,
                        wxSizerFlags( 1 ).Expand( ),
                        [this]( mwxBoxSizer& pxBoxSizer ) // BoxSizer content
                        {
                            pxBoxSizer.Add( new wxStaticText( pxBoxSizer.pxConnector.pxWindow, wxID_ANY,
                                                              wxT( "Reference genome" ), wxDefaultPosition,
                                                              wxDefaultSize, 0 ),
                                            0, wxLEFT | wxTOP, 5 );
                            pxBoxSizer.Add( this->xGenomeNameTextCtrl = new wxTextCtrl(
                                                pxBoxSizer.pxConnector.pxWindow, wxID_ANY, wxEmptyString,
                                                wxDefaultPosition, wxDefaultSize, wxTE_READONLY ),
                                            0, wxEXPAND | wxALL, 5 );
                        } ); // add vertical BoxSizer

                    // File Selector for genome selection
                    pxBoxSizer.Add( new mwxFileSelectDeleteButtonSizer(
                                        pxBoxSizer.pxConnector, "Genome selection", "Select Reference Genome",
                                        "Reference files (*.maReference)|*.maReference|All Files|*",
                                        std::bind( &MA_MainFrame::onGenomeSelection, this, std::placeholders::_1 ),
                                        false ), // no clear button
                                    wxSizerFlags( 0 ) );
                    this->onGenomeSelection( "" );
                } // lambda
                ) // add horizontal BoxSizer
            ; // add horizontal BoxSizer

        // Bottom part of frame
        pxFrameSizer->addBoxSizer(
            wxHORIZONTAL,
            wxSizerFlags( 1 ).Expand( ).Border( wxLEFT | wxRIGHT | wxBOTTOM, 10 ),
            [this]( mwxBoxSizer& pxBoxSizer ) {
                // Query Input Area (bsizer9)
                pxBoxSizer.addBoxSizer(
                    wxVERTICAL,
                    wxSizerFlags( 1 ).Expand( ),
                    [this]( mwxBoxSizer& pxBoxSizer ) { // BoxSizer content
                        pxBoxSizer.Add( new wxStaticText( pxBoxSizer.pxConnector.pxWindow, wxID_ANY,
                                                          wxT( "Query reads" ), wxDefaultPosition, wxDefaultSize, 0 ),
                                        0, wxTOP | wxLEFT, 5 );
                        this->xQueryTextCtrl =
                            new wxTextCtrl( pxBoxSizer.pxConnector.pxWindow, wxID_ANY, wxEmptyString, wxDefaultPosition,
                                            wxDefaultSize, wxHSCROLL | wxTE_MULTILINE );
                        pxBoxSizer.Add( this->xQueryTextCtrl, 1, wxALL | wxEXPAND, 5 );

                        this->xMatesStaticText =
                            new wxStaticText( pxBoxSizer.pxConnector.pxWindow, wxID_ANY, wxT( "Paired reads mates" ),
                                              wxDefaultPosition, wxDefaultSize, 0 );
                        pxBoxSizer.Add( this->xMatesStaticText, 0, wxTOP | wxLEFT, 5 );

                        this->xMatesTextCtrl = new wxTextCtrl( pxBoxSizer.pxConnector.pxWindow, wxID_ANY, wxEmptyString,
                                                               wxDefaultPosition, wxDefaultSize, 0 );
                        pxBoxSizer.Add( this->xMatesTextCtrl, 1, wxALL | wxEXPAND, 5 );

                        xMatesStaticText->Show( false );
                        xMatesTextCtrl->Show( false );
                    } ); // add vertical BoxSizer

                // File Selection and Start Button
                pxBoxSizer.addBoxSizer( // vertical BoxSizer (inside bottom of frame) bSizer14
                    wxVERTICAL,
                    wxSizerFlags( 0 ).Expand( ),
                    [this]( mwxBoxSizer& pxBoxSizer ) // BoxSizer content
                    {
                        // File Selection for Reads
                        pxBoxSizer.Add( new mwxFileSelectDeleteButtonSizer(
                                            pxBoxSizer.pxConnector, "Reads Selection", "Select Query File",
                                            "FASTA(Q) Files(*.fasta;*fastaq)|*.fasta;*.fastaq|All Files (*.*)|*.*",
                                            std::bind( &MA_MainFrame::onQuerySelection, this, std::placeholders::_1,
                                                       this->xQueryTextCtrl ) ),
                                        wxSizerFlags( 1 ) );

                        // Paired Read Mates File Selection
                        pxBoxSizer.addBoxSizer( // vertical BoxSizer bSizer12
                            wxVERTICAL,
                            wxSizerFlags( 1 ),
                            [this]( mwxBoxSizer& pxBoxSizer ) // BoxSizer content
                            {
                                pxBoxSizer.addBoxSizer( // vertical BoxSizer bSizer12
                                    wxVERTICAL,
                                    wxSizerFlags( 1 ),
                                    [this]( mwxBoxSizer& pxBoxSizer ) // Section with paired control elements
                                    {
                                        this->xMatesControlElements = &pxBoxSizer;
                                        pxBoxSizer.Add( new mwxFileSelectDeleteButtonSizer(
                                                            pxBoxSizer.pxConnector, "Mates Selection",
                                                            "Open Query Mate file",
                                                            "FASTA(Q) Files(*.fasta;*fastaq)|*.fasta;*.fastaq|All "
                                                            "Files (*.*)|*.*",
                                                            std::bind( &MA_MainFrame::onQuerySelection, this,
                                                                       std::placeholders::_1, this->xMatesTextCtrl ) ),
                                                        wxSizerFlags( 0 ) );
                                        pxBoxSizer.Add( 0, 0, 0, wxEXPAND | wxTOP, 5 ); // Vertical spacer
                                        pxBoxSizer.Add( new wxStaticText( pxBoxSizer.pxConnector.pxWindow, wxID_ANY,
                                                                          wxT( "Paired Settings" ), wxDefaultPosition,
                                                                          wxDefaultSize, 0 ),
                                                        0, wxALL | wxLEFT, 5 );
                                        // Gear button for paired settings
                                        pxBoxSizer.Add(
                                            new mwxBitmapButton( pxBoxSizer.pxConnector.pxWindow,
                                                                 wxBITMAP_PNG_FROM_DATA( GearPairedButtonNormal ),
                                                                 wxBITMAP_PNG_FROM_DATA( GearPairedButtonSelected ),
                                                                 std::bind( &MA_MainFrame::onPairedGearButton, this,
                                                                            std::placeholders::_1 ) ),
                                            0, wxALIGN_CENTER_HORIZONTAL, 5 );
                                        pxBoxSizer.Show( false );
                                    } );

                                // Start button
                                pxBoxSizer.Add( new mwxBitmapButton(
                                                    pxBoxSizer.pxConnector.pxWindow,
                                                    wxBITMAP_PNG_FROM_DATA( StartButton ),
                                                    wxBITMAP_PNG_FROM_DATA( StartButtonSelected ),
                                                    std::bind( &MA_MainFrame::onStart, this, std::placeholders::_1 ) ),
                                                0, wxALIGN_BOTTOM | wxALIGN_CENTER_HORIZONTAL, 5 );
                            } ); // addBoxSizer
                    } ); // add vertical BoxSizer
            } // lambda
        ); //  add horizontal BoxSizer

        SetSizer( pxFrameSizer );
        Layout( );
        Center( wxBOTH );
    } // constructor

    void OnQuit( wxCommandEvent& WXUNUSED( event ) )
    {
        // true is to force the frame to close
        Close( true );
    } // method
}; // class


// 'Main program' equivalent: the program execution "starts" here
bool MyApp::OnInit( )
{
#ifdef _MSC_VER
    AllocConsole( );
    freopen( "conin$", "r", stdin );
    freopen( "conout$", "w", stdout );
    freopen( "conout$", "w", stderr );
#endif
    // call the base class initialization method, currently it only parses a
    // few common command-line options but it could be do more in the future
    if( !wxApp::OnInit( ) )
        return false;

    // create the main application window
    MA_MainFrame* xFrame = new MA_MainFrame( "MA - The Modular Aligner" );

    // and show it (the frames, unlike simple controls, are not shown when
    // created initially)
    xFrame->Show( true );

    // success: wxApp::OnRun() will be called which will enter the main message
    // loop and the application will run. If we returned false here, the
    // application would exit immediately.
    return true;
} // method

// ----------------------------------------------------------------------------
// event tables and other macros for wxWidgets
// ----------------------------------------------------------------------------

// the event tables connect the wxWidgets events with the functions (event
// handlers) which process them. It can be also done at run-time, but for the
// simple menu events like this the static method is much simpler.
// clang-format off
//wxBEGIN_EVENT_TABLE( AlignFrame, wxFrame ) 
//     EVT_CLOSE( AlignFrame::OnClose ) 
//wxEND_EVENT_TABLE( )
// clang-format on
