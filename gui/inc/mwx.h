/////////////////////////////////////////////////////////////////////////////
// Name:        mwx.h
// Purpose:     MA wxWidgets Extensions
// Author:      Arne Kutzner
// Created:     01/01/2019
// Copyright:   (c) Arne Kutzner, Markus Schmidt
// Licence:		MIT
/////////////////////////////////////////////////////////////////////////////
#pragma once

/* Chosen Icons are from:
 * https://material.io/tools/icons/?icon=cancel&style=outline
 * https://fontawesome.com/icons?d=gallery
 */

/* for all others, include the necessary headers (this file is usually all you
 * need because it includes almost all "standard" wxWidgets headers)
 */
#ifndef WX_PRECOMP
#include "wx/wx.h"
#include <wx/filepicker.h>
#include <wx/notebook.h>
#include <wx/textctrl.h>
#endif

#include <memory>
#include <vector>

#include "util/parameter.h"
// Include all bitmaps (Note: All bitmaps are static)
#include "png_bitmaps.c"

/* Global execution context */
class mxwExecutionContext
{
  public:
    static int iHighestID;

    /* Deliver a unique ID)*/
    static int uniqueID( void )
    {
        return iHighestID++;
    } // method
}; // class


/* Menu Item */
class mwxMenuItem
{
  public:
    mwxMenuItem( const std::string& ItemText ) : ItemText( ItemText )
    {} // default constructor

    std::string ItemText;
}; // class


/* Managed variant of wxMenuBar */
class mwxMenu : public wxMenu
{
  public:
    mwxMenu( const std::string& ItemText ) : wxMenu( )
    {} // constructor

    virtual ~mwxMenu( )
    {
        // std::cout << "In ~mwxMenu()" << std::endl;
    } // destructor

    /* All menu-items of the menu */
    std::vector<mwxMenuItem> vChilds;

    /* Appends an item to the menu. */
    mwxMenu& append( const std::string& sMenuItemText, // text of menu-item
                     std::function<void( wxCommandEvent& )> handler = []( wxCommandEvent& ) {}, // handler
                     const std::string& sHelpText = "" ) // helper text for menu item
    {
        auto iItemID = mxwExecutionContext::uniqueID( );
        vChilds.emplace_back( mwxMenuItem( sMenuItemText ) );
        this->Append( iItemID, sMenuItemText.c_str( ), sHelpText.c_str( ) );
        wxEvtHandler::Bind( wxEVT_MENU, handler, iItemID );
        return *this;
    } // method
}; // class


/* Managed variant of wxMenuBar.
 * All menus must be created on the heap because all menus attached
 * to a menu-bar or to another menu will be deleted by their
 * parent when it is deleted.
 */
class mwxMenuBar : public wxMenuBar
{
  public:
    mwxMenuBar( ) : wxMenuBar( )
    {} // constructor

    virtual ~mwxMenuBar( )
    {
        // std::cout << "In ~mwxMenuBar()" << std::endl;
    } // destructor

    /* All menus of the menu-bar */
    std::vector<mwxMenu*> vChilds;

    /* Deliver reference to menu item */
    mwxMenu& menu( int iMenuId = -1 )
    {
        return iMenuId < 0 ? *( vChilds.back( ) ) : *( vChilds.at( iMenuId ) );
    } // method

    /* Appends a new( empty ) menu */
    mwxMenuBar& push_back( const std::string& sMenuText )
    {
        vChilds.emplace_back( new mwxMenu( sMenuText ) );
        this->Append( vChilds.back( ), sMenuText );
        return *this;
    } // method
}; // class


/* Connects widgets via windows and sizers */
class mwxConnector
{
  public:
    wxWindow* const pxWindow;
    wxSizer* const pxSizer;

    mwxConnector( wxWindow* pxWindow, wxSizer* pxSizer ) : pxWindow( pxWindow ), pxSizer( pxSizer )
    {}
}; // class


/* Extended  wxStaticBoxSizer */
class mwxStaticBoxSizer : public wxStaticBoxSizer
{
  public:
    const mwxConnector xConnector; // connector for the current box sizer

    /* constructor */
    mwxStaticBoxSizer( const mwxConnector& xConnector, // connector to host
                       const std::string& sBoxText, // label text of box
                       int iOrientation ) // box orientation
        : wxStaticBoxSizer( new wxStaticBox( xConnector.pxWindow, wxID_ANY, sBoxText ), iOrientation ),
          xConnector( this->GetStaticBox( ), nullptr )
    {} // constructor

    /* Add via lambda */
    mwxStaticBoxSizer* add( const wxSizerFlags& flags, // sizer flags for fresh content
                            std::function<wxWindow*( const mwxConnector& )>
                                fMakeContent )
    {
        this->Add( fMakeContent( this->xConnector ), flags );
        return this;
    } // method
}; // class


/* Managed wxBoxSizer */
class mwxBoxSizer : public wxBoxSizer
{
  public:
    const mwxConnector pxConnector; // connector for the current box sizer

    /* constructor */
    mwxBoxSizer( const mwxConnector& pxConnector, // connector for BoxSizer
                 int iOrientation )
        : wxBoxSizer( iOrientation ), pxConnector( pxConnector.pxWindow, this )
    {} // method

    /* Add via lambda */
    mwxBoxSizer* add( const wxSizerFlags& flags, // sizer flags for fresh content
                      std::function<wxWindow*( const mwxConnector& )>
                          fMakeContent )
    {
        this->Add( fMakeContent( this->pxConnector ), flags );
        return this;
    } // method

    /* Add via lambda */
    mwxBoxSizer* add( const wxSizerFlags& flags, // sizer flags for fresh content
                      std::function<wxSizer*( const mwxConnector& )>
                          fMakeContent )
    {
        this->Add( fMakeContent( this->pxConnector ), flags );
        return this;
    } // method

    mwxBoxSizer* addBoxSizer( int iOrientation, // orientation of the new bxoSizer
                              const wxSizerFlags& flags, // sizer flags for fresh content
                              std::function<void( mwxBoxSizer& )>
                                  fMakeContent )
    {
        auto* pxBoxSizer = new mwxBoxSizer( pxConnector, iOrientation );
        this->Add( pxBoxSizer, flags );
        fMakeContent( *pxBoxSizer );
        return this;
    } // method

    mwxBoxSizer* addStaticBoxSizer( const std::string& sBoxLabel, // label of static sizer
                                    int iOrientation, // orientation of the new bxoSizer

                                    const wxSizerFlags& flags, // sizer flags for fresh content
                                    std::function<void( mwxStaticBoxSizer& )>
                                        fMakeContent )
    {
        auto* pxStaticBoxSizer = new mwxStaticBoxSizer( pxConnector, sBoxLabel, iOrientation );
        this->Add( pxStaticBoxSizer, flags );
        fMakeContent( *pxStaticBoxSizer );
        return this;
    } // method
}; // class


/* Managed variant of bitmap button */
class mwxBitmapButton : public wxBitmapButton
{
  public:
    wxSize calculateSize( const wxBitmap& rxBitmapNormal, const int iWidthOverwrite )
    {
        // All bitmaps are expected to have an height of 35
        int iWidth = std::max( 48, rxBitmapNormal.GetWidth( ) ) + 12;
        int iHeight = std::max( 35, rxBitmapNormal.GetHeight( ) ) + 12;
        return wxSize( std::max( iWidth, iWidthOverwrite ), iHeight );
    } // method
    /* Constructor */
    mwxBitmapButton( wxWindow* pxHostWindow, // host window of box context (responsible for destruction)
                     const wxBitmap& rxBitmapNormal, // Bitmap unpressed
                     std::function<void( wxCommandEvent& )> handler = []( wxCommandEvent& ) {}, // handler
                     const wxPoint rxPos = wxDefaultPosition,
                     const int iWidthOverwrite = 0 )
        : wxBitmapButton( pxHostWindow, wxID_ANY, rxBitmapNormal, rxPos,
                          calculateSize( rxBitmapNormal, iWidthOverwrite ) ) //, wxNO_BORDER )
    {
        wxEvtHandler::Bind( wxEVT_BUTTON, handler, wxID_ANY );
    } // constructor
}; // class


/* Classic OK, Cancel Dialog having another window as content.
 * Returns wxID_OK if OK button was pressed, wxID_CANCEL for cancel button.
 */
class mwxOK_Cancel_Dialog : public wxDialog
{
  public:
    /* Constructor */
    mwxOK_Cancel_Dialog( wxWindow* pxHostWindow, // host window of box context (responsible for destruction)
                         const wxString& sTitle, // title of the dialog
                         std::function<wxWindow*( wxWindow* )>
                             fMakeContent, // called for creating the content of the dialog
                         const wxPoint& xPos = wxDefaultPosition, // position of the dialog
                         const std::string& rsKOButtonLabel = "OK", // Label of OK Button
                         bool bDoFit = false, // adapt the size of the dialog to the size of its content
                         std::function<void( wxCommandEvent& )> OK_handler =
                             []( wxCommandEvent& ) {} ) // handler of OK button

        : wxDialog( pxHostWindow, wxID_ANY, sTitle, xPos, wxSize( 500, 400 ),
                    ( wxDEFAULT_FRAME_STYLE | wxTAB_TRAVERSAL ) /* ^ wxRESIZE_BORDER */ ^ wxMAXIMIZE_BOX ^
                        wxMINIMIZE_BOX ^ ( bDoFit ? wxRESIZE_BORDER : 0 ) )
    {
        SetSizeHints( wxDefaultSize, wxDefaultSize );
        SetBackgroundColour( wxColor( 255, 255, 255 ) );

        wxBoxSizer* pxTopSizer = new wxBoxSizer( wxVERTICAL );

        // Add inner content of dialog
        // std::cout << "Before mwxOK_Cancel_Dialog " << std::endl;
        pxTopSizer->Add( fMakeContent( this ), wxSizerFlags( 1 ).Align( wxALIGN_CENTER ).Expand( ).Border( wxALL, 5 ) );
        // std::cout << "After mwxOK_Cancel_Dialog " << std::endl;

        wxBoxSizer* pxButtomSizer = new wxBoxSizer( wxHORIZONTAL );

        // Create two buttons that are horizontally unstretchable,
        // with an all-around border with a width of 10 and implicit top alignment
        // auto iOK_Button_ID = mxwExecutionContext::uniqueID( );
        pxButtomSizer->Add( new wxButton( this, wxID_OK, rsKOButtonLabel ),
                            wxSizerFlags( 0 ).Align( wxALIGN_CENTER ).Border( wxALL, 5 ) );

        pxButtomSizer->Add( new wxButton( this, wxID_CANCEL, "Cancel" ),
                            wxSizerFlags( 0 ).Align( wxALIGN_CENTER ).Border( wxALL, 5 ) );

        // Create a sizer with no border and centered horizontally
        pxTopSizer->Add( pxButtomSizer, wxSizerFlags( 0 ).Center( ) );

        SetSizer( pxTopSizer );
        Layout( );
        if( bDoFit )
            // 'Fit' fits the dialog to the size of its children
            // Adapt the size of the dialog to the size of its content
            Fit( );
        Centre( wxBOTH );
    } // constructor
}; // class


/* Pair of button for FileDialog and clear-button */
class mwxFileSelectDeleteButtonSizer : public mwxBoxSizer
{
  public:
    const mwxConnector pxConnector; // connector for the current box sizer
    const bool bMultiFileSelect;

  private:
    const std::string sFileDialogTitle; // "Open Reference file"
    const std::string sFilePattern; // "Reference files (*.maReference)|*.maReference|All Files|*"
    std::function<void( const std::vector<fs::path>& )> fHandler;
    std::function<void( )> fClearHandler = [this]( void ) { this->fHandler( std::vector<fs::path>{} ); };

    void onFolderButton( wxCommandEvent& WXUNUSED( event ) )
    {
        wxFileDialog openFileDialog( pxConnector.pxWindow, sFileDialogTitle, "", "", sFilePattern,
                                     wxFD_OPEN | wxFD_FILE_MUST_EXIST |
                                         ( this->bMultiFileSelect ? wxFD_MULTIPLE : 0 ) );
        if( openFileDialog.ShowModal( ) == wxID_OK )
        {
            // Call handler to inform about new selection
            std::vector<fs::path> vxRet;
            wxArrayString xFilenames;
            openFileDialog.GetPaths( xFilenames );
            for( auto& rsPath : xFilenames )
                vxRet.push_back( fs::path( std::string( rsPath.c_str( ) ) ) );
            fHandler( vxRet );
        } // if
    } // method

    /* Event handler for clear button */
    void onClearButton( wxCommandEvent& WXUNUSED( event ) )
    {
        // Clear Button results in call of handler with empty string
        fClearHandler( );
    } // method

  public:
    /* Primary Constructor */
    mwxFileSelectDeleteButtonSizer( const mwxConnector& pxConnector, // connector for BoxSizer
                                    const std::string& rsText, // text of the button group
                                    const std::string& rsFileDialogTitle,
                                    const std::string& rsFilePattern,
                                    const bool bMultiFileSelect = false,
                                    std::function<void( const std::vector<fs::path>& )> fHandler =
                                        []( const std::vector<fs::path>& ) {},
                                    const unsigned int uiIconId = 0 )
        : mwxBoxSizer( pxConnector, wxVERTICAL ),
          pxConnector( pxConnector ),
          bMultiFileSelect( bMultiFileSelect ),
          sFileDialogTitle( rsFileDialogTitle ),
          sFilePattern( rsFilePattern ),
          fHandler( fHandler )
    {
        // const unsigned char* pIconBitmap = uiIconId == 0 ? DeleteButtonNormal_png : HammerIcon_png;
        Add( new wxStaticText( pxConnector.pxWindow, wxID_ANY, rsText, wxDefaultPosition, wxDefaultSize, 0 ), 0,
             wxALL | wxLEFT, 5 );

        addBoxSizer( // horizontal BoxSizer (File Selection and Deletion)
            wxHORIZONTAL,
            wxSizerFlags( 0 ).ReserveSpaceEvenIfHidden( ),
            [&]( mwxBoxSizer& pxBoxSizer ) // BoxSizer content
            {
                pxBoxSizer.Add( new mwxBitmapButton( pxBoxSizer.pxConnector.pxWindow,
                                                     wxBITMAP_PNG_FROM_DATA( FolderButtonNormal ),
                                                     std::bind( &mwxFileSelectDeleteButtonSizer::onFolderButton, this,
                                                                std::placeholders::_1 ) ),
                                0, wxRESERVE_SPACE_EVEN_IF_HIDDEN | wxRIGHT | wxLEFT, 5 ); // GTK Bug Border 0

                mwxBitmapButton* pxClearButton;
                pxBoxSizer.Add( pxClearButton =
                                    new mwxBitmapButton( pxBoxSizer.pxConnector.pxWindow,
                                                         uiIconId == 0 ? wxBITMAP_PNG_FROM_DATA( DeleteButtonNormal )
                                                                       : wxBITMAP_PNG_FROM_DATA( BuildIcon ),
                                                         std::bind( &mwxFileSelectDeleteButtonSizer::onClearButton,
                                                                    this, std::placeholders::_1 ) ),
                                0, wxRESERVE_SPACE_EVEN_IF_HIDDEN | wxRIGHT | wxLEFT, 5 ); // GTK Bug Border 0
            } ); // addBoxSizer

        fHandler( std::vector<fs::path>{} );
    } // constructor

    mwxFileSelectDeleteButtonSizer* setClearHandler( std::function<void( )> fClearHandler )
    {
        this->fClearHandler = fClearHandler;
        return this;
    } // method


}; // class


/* Map driven combo box.
 * Creator has to guarantee that the map stays unchanged during objects lifetime.
 */
class mwxMapDrivenComboBox : public wxComboBox
{

  public:
    wxArrayString xChoices; // Different choices of combo-box

    // Constructor
    template <typename TYPE>
    mwxMapDrivenComboBox(
        wxWindow* pxHost,
        const std::map<std::string, TYPE>& rxMap,
        std::function<void( wxCommandEvent& )> handler = []( wxCommandEvent& ) {}, // handler
        std::function<std::string( typename std::map<std::string, TYPE>::const_iterator )> fNameFromKeyValuePair =
            []( typename std::map<std::string, TYPE>::const_iterator xIt ) { return xIt->first; } )
        : wxComboBox( )
    {
        for( typename std::map<std::string, TYPE>::const_iterator xIt = rxMap.begin( ); xIt != rxMap.end( ); ++xIt )
        {
            xChoices.Add( fNameFromKeyValuePair( xIt ).c_str( ) ); // rxKeyValue.second->sName.c_str( )
        } // for
        /*
         * Define order so that 'Custom' is always last and 'Default' always first.
         * Other strings are sorted alphabetically.
         *
         * Widgets documentation: CompareFunction is defined as a function taking two const wxString& parameters and
         * returning an int value less than, equal to or greater than 0 if the first string is less than, equal to or
         * greater than the second one.
         *
         */
        xChoices.Sort( []( const wxString& rsA, const wxString& rsB ) {
            if( rsA == rsB ) // strings are equal
                return 0;
            if( rsA == "Default" || rsB == "Custom" ) // sA is smaller than sB
                return -1;
            if( rsA == "Custom" || rsB == "Default" ) // sA is larger than sB
                return 1;
            if( rsA < rsB ) // determine which string is larger lexigraphically
                return -1;
            return 1;
        } );

        this->Create( pxHost, wxID_ANY, xChoices[ 0 ], wxDefaultPosition, wxDefaultSize, xChoices, wxCB_READONLY );
        // Bind handler that process changes
        this->wxEvtHandler::Bind( wxEVT_COMBOBOX, handler );
    } // constructor

    void setSelected( const std::string& rsSelection )
    {
        // try out all choices of combo box
        for( size_t uiI = 0; uiI < xChoices.Count( ); uiI++ )
            if( std::string( xChoices[ uiI ].c_str( ) ) == rsSelection )
            {
                this->SetSelection( uiI );
                // return as soon as we find a matching choice
                return;
            } // if
        // if none match throw exception
        throw std::runtime_error( "Could not find " + rsSelection + " in combobox." );
    } // method
}; // class


/* Managed static box context */
class mwxStaticBoxContext : public wxPanel
{
  public:
    wxStaticBox* pxStaticBox;
    wxStaticBoxSizer* pxStaticBoxSizer;
    wxSizer* pxSizerHost;

    mwxStaticBoxContext( wxWindow* pxWindowHost, //, // host window of box context
                         wxSizer* pxSizerHost, // sizer of the host window
                         const std::string& sBoxLabel ) // label of the static box
        : wxPanel( pxWindowHost, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL ),
          pxSizerHost( pxSizerHost )
    {
        // deallocation controlled by this (wxPanel)
        pxStaticBox = new wxStaticBox( this, wxID_ANY, sBoxLabel.c_str( ) );
        // deallocation controlled by this (wxPanel)
        pxStaticBoxSizer = new wxStaticBoxSizer( this->pxStaticBox, wxVERTICAL );

        this->SetSizer( this->pxStaticBoxSizer );
        pxSizerHost->Add( this, 1, wxEXPAND | wxALL, 5 );
    } // constructor

    /* Delivers the windows for connections with other widgets */
    wxStaticBox* getWindow( void )
    {
        return this->pxStaticBoxSizer->GetStaticBox( );
    } // method

    /* Delivers the sizer for connections with other widgets */
    wxStaticBoxSizer* getSizer( void )
    {
        return pxStaticBoxSizer;
    } // method

    /* Delivers a connector for the current static box panel */
    mwxConnector getConnector( void )
    {
        return mwxConnector( pxStaticBoxSizer->GetStaticBox( ), pxStaticBoxSizer );
    } // method

    void layout( void )
    {
        this->Layout( );
        this->pxStaticBoxSizer->Fit( this );
    } // method

    virtual ~mwxStaticBoxContext( )
    {
        // std::cout << "In ~mwxStaticBoxContext()" << std::endl;
    } // destructor
}; // class


/* Scrolled Window that comprises several static box panels.
 * The topmost panel avoids problem with GTK.
 */
class mwxScrolledStatixBoxesContext : public wxPanel
{
  public:
    wxBoxSizer* pxBoxSizer;
    wxScrolledWindow* pxScrolledWindow; // Inner scrolled Window
    wxFlexGridSizer* pxFlexGridSizer;
    std::vector<mwxStaticBoxContext*> vpStaticBoxes; // Vector of pointer to boxes

    /* Constructor*/
    mwxScrolledStatixBoxesContext( wxWindow* pxHostWindow, // host window of box context (responsible for destruction)
                                   wxSizer* pxHostSizer ) // sizer of the host window
        : wxPanel( pxHostWindow, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                   wxTAB_TRAVERSAL ), // base class constructor call
          pxBoxSizer( new wxBoxSizer( wxVERTICAL ) ),
          pxScrolledWindow(
              // wxSize( 100, 100 ) is required for avoiding GTK errors.
              new wxScrolledWindow( this, wxID_ANY, wxDefaultPosition, wxSize( 100, 100 ), wxHSCROLL | wxVSCROLL ) ),
          pxFlexGridSizer( new wxFlexGridSizer( 0, 1, 0, 0 ) )
    {
        pxFlexGridSizer->SetFlexibleDirection( wxBOTH );
        pxFlexGridSizer->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );

        pxScrolledWindow->SetSizer( pxFlexGridSizer );
        pxScrolledWindow->SetScrollRate( 5, 5 ); // enables scrollbar
        pxScrolledWindow->Layout( );
        pxBoxSizer->Add( pxScrolledWindow, 1, wxEXPAND | wxALL, 5 );

        SetSizer( pxBoxSizer );
        if( pxHostSizer != NULL )
            pxHostSizer->Add( this, 1, wxEXPAND | wxALL, 5 );
    } // constructor

    /* Delegating constructor */
    mwxScrolledStatixBoxesContext( const mwxConnector& xConnector )
        : mwxScrolledStatixBoxesContext( xConnector.pxWindow, xConnector.pxSizer )
    {} // constructor

    /* Creates a new static box and returns a connector */
    mwxStaticBoxContext* addStaticBox( )
    {
        auto* pxStaticBoxContext = new mwxStaticBoxContext( pxScrolledWindow, pxFlexGridSizer, "Settings" );
        vpStaticBoxes.push_back( pxStaticBoxContext );
        return pxStaticBoxContext;
    } // class

    void layout( void )
    {
        Layout( );
        pxBoxSizer->Fit( this );
    } // method
}; // class


/* Property notebook consisting of several property panels (pages) */
class mwxPropertyNotebook : public wxNotebook
{
  public:
    std::vector<mwxScrolledStatixBoxesContext*> vPages; // each page is a scrolled static boxes context

    /* Constructor */
    mwxPropertyNotebook( wxWindow* pxHostWindow ) // host window of the notebook
        : wxNotebook( pxHostWindow, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 )
    {} // constructor

    mwxScrolledStatixBoxesContext* addPage( const std::string& sPageTitle )
    {
        // Create static box context for fresh page
        vPages.push_back( new mwxScrolledStatixBoxesContext( this, nullptr ) );
        this->AddPage( vPages.back( ), sPageTitle.c_str( ), vPages.size( ) == 1 ? true : false );
        return vPages.back( );
    } // method
}; // class

/* Managed wxWidgets extension for property management.
 * FIXME: Extend pxFlexLayout in order to get deallocation problem solved.
 */
class mwxPropertyPanel : public wxPanel
{
    /* Inner class describing text property */
    class mwxProperty
    {
      public:
        mwxPropertyPanel* pxHost; // Parent of property (manages destruction of text fields)
        wxStaticText* pxStaticText;

        /* Constructor */
        mwxProperty( mwxPropertyPanel* pxHost, // parent of property
                     const std::string& sName ) // name of property
            : pxHost( pxHost )
        {
            /* Parent of host will manage destruction of wxStaticText */
            pxHost->pxFlexLayout->Add( pxStaticText = new wxStaticText( pxHost, wxID_ANY, sName.c_str( ),
                                                                        wxDefaultPosition, wxDefaultSize, 0 ),
                                       0, wxALL | wxALIGN_CENTER_VERTICAL, 5 );
        } // constructor

        virtual ~mwxProperty( )
        {} // destructor

        /* Update visual properties */
        virtual void updateEnabledDisabledDispatch( void )
        {} // method

        void updateEnabledDisabled( bool bEnabled, wxControl* pxControl )
        {
            if( bEnabled )
            {
                pxStaticText->Enable( );
                pxControl->Enable( );
            } // if
            else
            {
                pxStaticText->Disable( );
                pxControl->Disable( );
            } // else
        } // method
    }; // inner class

    /* Inner class describing text property */
    template <typename VALUE_TYPE> class mwxTextProperty : public mwxProperty
    {
      public:
        wxTextCtrl* pxTextCtrl; // Text control holding the actual value
        std::shared_ptr<AlignerParameter<VALUE_TYPE>> pxParameter = nullptr; // pointer to parameter object

        /* Update the value of the property */
        void update( void )
        {
            pxTextCtrl->SetValue( this->pxParameter->asText( ) );
        } // method

        /* Handler is called whenever the user leaves the input field (the user ends property editing) */
        void focusHandler( wxFocusEvent& rxEvent )
        {
            /* Wrong input (out of range etc.) are reported via exceptions */
            try
            {
                pxParameter->setByText( std::string( this->pxTextCtrl->GetValue( ).c_str( ) ) );
            } // try
            catch( std::exception& xException )
            {
                wxMessageBox( xException.what( ), wxT( "Wrong Input" ), wxICON_ERROR );
            } // catch
            update( );
            /* From wxWidgets documentation:
             * The focus event handlers should almost invariably call
             * wxEvent::Skip() on their event argument to allow the default handling to take place.
             */
            rxEvent.Skip( );

            pxHost->updateEnabledDisabled( );
        } // method

        /* Constructor with AlignerParameter */
        mwxTextProperty( mwxPropertyPanel* pxHost, // parent of property
                         std::shared_ptr<AlignerParameter<VALUE_TYPE>>
                             pxParameter // parameter reference
                         )
            : mwxProperty( pxHost, pxParameter->sName ),
              pxTextCtrl( new wxTextCtrl( pxHost, wxID_ANY, pxParameter->asText( ), wxDefaultPosition, wxDefaultSize,
                                          wxTE_PROCESS_ENTER ) ),
              pxParameter( pxParameter )
        {
            pxHost->pxFlexLayout->Add( pxTextCtrl, 0, wxALL, 5 );

            /* The description of the parameter is used as ToolTip-text */
            this->pxTextCtrl->SetToolTip( pxParameter->sDescription );
            this->pxStaticText->SetToolTip( pxParameter->sDescription );

            /* Bind "focus leave" to the handler */
            this->pxTextCtrl->wxEvtHandler::Bind( wxEVT_KILL_FOCUS,
                                                  [this]( wxFocusEvent& rxEvent ) { this->focusHandler( rxEvent ); } );
            /* Make the enter key to behave like the tab key. */
            this->pxTextCtrl->wxEvtHandler::Bind( wxEVT_COMMAND_TEXT_ENTER,
                                                  [&]( wxCommandEvent& rxEvent ) { this->pxTextCtrl->Navigate( ); } );
        } // constructor

        virtual void updateEnabledDisabledDispatch( void )
        {
            if( this->pxParameter )
                updateEnabledDisabled( this->pxParameter->fEnabled( ), pxTextCtrl );
        } // method
    }; // inner class

    /* Inner class describing text property */
    class mwxFolderProperty : public mwxProperty
    {
      public:
        mwxPropertyPanel* pxHost; // keep host for folder dialog creation
        std::shared_ptr<AlignerParameter<fs::path>> pxParameter = nullptr; // pointer to parameter object
        wxDirPickerCtrl* pxDirPickerCtrl; // Folder picker

        /* Handler is called whenever the user leaves the input field (the user ends property editing) */
        void handler( wxCommandEvent& rxEvent )
        {
            fs::path sPath( std::string( this->pxDirPickerCtrl->GetPath( ).c_str( ) ) );

            // Check folder for existence
            if( fs::exists( sPath ) )
                pxParameter->set( sPath );
            else
                wxMessageBox( wxString( std::string( "Folder " ) + sPath.string( ) + " does not exist." ), "Error",
                              wxICON_ERROR );

            pxHost->updateEnabledDisabled( );
        } // method

        /* Constructor for AlignerParameter */
        mwxFolderProperty( mwxPropertyPanel* pxHost, // parent of property
                           std::shared_ptr<AlignerParameter<fs::path>>
                               pxParameter // parameter reference
                           )
            : mwxProperty( pxHost, pxParameter->sName ),
              pxHost( pxHost ),
              pxParameter( pxParameter ),
              pxDirPickerCtrl( new wxDirPickerCtrl( pxHost, wxID_ANY, pxParameter->get( ).string( ),
                                                    wxT( "Select a folder" ), wxDefaultPosition, wxDefaultSize,
                                                    wxDIRP_DEFAULT_STYLE ) )
        {
            pxHost->pxFlexLayout->Add( pxDirPickerCtrl, 0, wxALL, 5 );

            this->pxDirPickerCtrl->SetToolTip( pxParameter->sDescription );
            this->pxStaticText->SetToolTip( pxParameter->sDescription );

            // Bind "focus leave" to the handler
            this->pxDirPickerCtrl->wxEvtHandler::Bind(
                wxEVT_DIRPICKER_CHANGED, [this]( wxCommandEvent& rxEvent ) { this->handler( rxEvent ); } );
        } // constructor

        virtual void updateEnabledDisabledDispatch( void )
        {
            if( this->pxParameter )
                updateEnabledDisabled( this->pxParameter->fEnabled( ), pxDirPickerCtrl );
        } // method
    }; // inner class

    /* Inner class describing text property */
    class mwxFileProperty : public mwxProperty
    {
      public:
        mwxPropertyPanel* pxHost; // keep host for file dialog creation
        std::shared_ptr<AlignerParameter<fs::path>> pxParameter = nullptr; // pointer to parameter object
        wxFilePickerCtrl* pxFilePickerCtrl; // File picker

        /* Handler is called whenever the user leaves the input field (the user ends property editing) */
        void handler( wxCommandEvent& rxEvent )
        {
            fs::path sFile( std::string( this->pxFilePickerCtrl->GetFileName( ).GetFullPath( ).c_str( ) ) );

            // @todo Markus: here should be a check on wether the extension fo the selected/written file matches...
            pxParameter->set( sFile );

            pxHost->updateEnabledDisabled( );
        } // method

        /* Constructor for AlignerParameter */
        mwxFileProperty( mwxPropertyPanel* pxHost, // parent of property
                         std::shared_ptr<AlignerParameter<fs::path>>
                             pxParameter, // parameter reference
                         const std::string& rsSuffixWildcardString // allowed file suffixes
                         )
            : mwxProperty( pxHost, pxParameter->sName ),
              pxHost( pxHost ),
              pxParameter( pxParameter ),
              pxFilePickerCtrl( new wxFilePickerCtrl( pxHost, wxID_ANY, pxParameter->get( ).string( ),
                                                      wxT( "Select a file" ), rsSuffixWildcardString, wxDefaultPosition,
                                                      wxDefaultSize, wxDIRP_DEFAULT_STYLE ) )
        {
            pxHost->pxFlexLayout->Add( pxFilePickerCtrl, 0, wxALL, 5 );

            this->pxFilePickerCtrl->SetToolTip( pxParameter->sDescription );
            this->pxStaticText->SetToolTip( pxParameter->sDescription );

            // Bind "focus leave" to the handler
            this->pxFilePickerCtrl->wxEvtHandler::Bind(
                wxEVT_FILEPICKER_CHANGED, [this]( wxCommandEvent& rxEvent ) { this->handler( rxEvent ); } );
        } // constructor

        virtual void updateEnabledDisabledDispatch( void )
        {
            if( this->pxParameter )
                updateEnabledDisabled( this->pxParameter->fEnabled( ), pxFilePickerCtrl );
        } // method
    }; // inner class


    /* Inner class describing choice property */
    class mwxChoiceProperty : public mwxProperty
    {
      public:
        std::shared_ptr<AlignerParameter<AlignerParameterBase::ChoicesType>> pxParameter =
            nullptr; // pointer to parameter object
        wxComboBox* pxCombo; // Text control holding the actual value

        /* Update the value of the property */
        void update( void )
        {
            pxCombo->SetSelection( pxParameter->uiSelection );
        } // method

        /* Handler is called whenever the user leaves the input field (the user ends property editing) */
        void selectionHandler( wxCommandEvent& rxEvent )
        {
            // Wrong input (out of range etc.) are reported via exceptions
            try
            {
                pxParameter->set( pxCombo->GetSelection( ) );
            } // try
            catch( std::exception& xException )
            {
                wxMessageBox( xException.what( ), wxT( "Wrong Input" ), wxICON_ERROR );
            } // catch
            update( );
            // From wxWidgets documentation:
            // The focus event handlers should almost invariably call
            // wxEvent::Skip() on their event argument to allow the default handling to take place.
            // rxEvent.Skip( );

            pxHost->updateEnabledDisabled( );
        } // method

        /* Constructor with AlignerParameter */
        mwxChoiceProperty( mwxPropertyPanel* pxHost, // parent of property
                           std::shared_ptr<AlignerParameter<AlignerParameterBase::ChoicesType>>
                               pxParameter // parameter reference
                           )
            : mwxProperty( pxHost, pxParameter->sName ), pxParameter( pxParameter )
        {
            /* Prepare required string array of wxWidgets and create the choice widget */
            wxArrayString aChoices;
            for( auto& sChoiceText : pxParameter->vChoices )
                aChoices.Add( sChoiceText.second.c_str( ) );
            this->pxCombo = new wxComboBox( pxHost, wxID_ANY, aChoices[ pxParameter->uiSelection ], wxDefaultPosition,
                                            wxDefaultSize, aChoices, wxCB_READONLY );
            pxHost->pxFlexLayout->Add( pxCombo, 0, wxALL, 5 );

            // Description of the parameter is used as ToolTip-text
            this->pxCombo->SetToolTip( pxParameter->sDescription );
            this->pxStaticText->SetToolTip( pxParameter->sDescription );

            // Bind handler that process changes
            this->pxCombo->wxEvtHandler::Bind(
                wxEVT_COMBOBOX, [this]( wxCommandEvent& rxEvent ) { this->selectionHandler( rxEvent ); } );
        } // constructor

        virtual void updateEnabledDisabledDispatch( void )
        {
            if( this->pxParameter )
                updateEnabledDisabled( this->pxParameter->fEnabled( ), pxCombo );
        } // method
    }; // inner class


    /* Inner class describing choice property */
    class mwxCheckBoxProperty : public mwxProperty
    {
      public:
        std::shared_ptr<AlignerParameter<bool>> pxParameter = nullptr; // pointer to parameter object
        wxCheckBox* pxCheckBox; // Text control holding the actual value

        /* Handler is called whenever the user leaves the input field (the user ends property editing) */
        void commandHandler( wxCommandEvent& rxEvent )
        {
            // Copy boolean value to parameter
            pxParameter->set( pxCheckBox->GetValue( ) );
            pxHost->updateEnabledDisabled( );
        } // method

        mwxCheckBoxProperty( mwxPropertyPanel* pxHost, // parent of property
                             std::shared_ptr<AlignerParameter<bool>>
                                 pxParameter // parameter reference
                             )
            : // pxHost( pxHost ),
              mwxProperty( pxHost, pxParameter->sName ),
              pxParameter( pxParameter ),
              pxCheckBox( new wxCheckBox( pxHost, wxID_ANY, wxT( "" ), wxDefaultPosition, wxDefaultSize, 0 ) )
        {
            // Add check-box to layout and set initial value
            pxHost->pxFlexLayout->Add( pxCheckBox, 0, wxALL, 5 );
            this->pxCheckBox->SetValue( pxParameter->get( ) );

            // Description of the parameter is used as ToolTip-text
            this->pxCheckBox->SetToolTip( pxParameter->sDescription );
            this->pxStaticText->SetToolTip( pxParameter->sDescription );

            // Bind handler that process changes
            this->pxCheckBox->wxEvtHandler::Bind(
                wxEVT_CHECKBOX, [this]( wxCommandEvent& rxEvent ) { this->commandHandler( rxEvent ); } );
        } // constructor

        virtual void updateEnabledDisabledDispatch( void )
        {
            if( this->pxParameter )
                updateEnabledDisabled( this->pxParameter->fEnabled( ), pxCheckBox );
        } // method
    }; // inner class

  public:
    wxWindow* pxHostWindow; // Host window
    wxFlexGridSizer* pxFlexLayout; // Layout of property editor
    std::vector<std::shared_ptr<mwxProperty>> vProperties; // Vector of all properties

    /* Constructor requires parent widgets */
    mwxPropertyPanel( wxWindow* pxHostWindow, // host window of box context (responsible for destruction)
                      wxSizer* pxHostSizer ) // sizer of the host window
        : wxPanel( pxHostWindow, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                   wxTAB_TRAVERSAL ), // base class constructor call
          pxHostWindow( pxHostWindow ),
          pxFlexLayout( new wxFlexGridSizer( 0, 2, 0, 0 ) )
    {
        this->pxFlexLayout->SetFlexibleDirection( wxBOTH ); // according to wxFormBuild
        this->pxFlexLayout->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED ); // according to wxFormBuild
        this->SetSizer( this->pxFlexLayout );
        if( pxHostSizer != NULL )
            pxHostSizer->Add( this, 1, wxEXPAND | wxALL, 5 );
    } // constructor

    /* Delegating constructor */
    mwxPropertyPanel( const mwxConnector& xConnector ) : mwxPropertyPanel( xConnector.pxWindow, xConnector.pxSizer )
    {} // constructor

    /* Add value property */
    mwxPropertyPanel& append( std::shared_ptr<AlignerParameter<int>>& pParameter )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxTextProperty<int>>( this, pParameter ) );

        return *this;
    } // method

    /* Add value property */
    mwxPropertyPanel& append( std::shared_ptr<AlignerParameter<bool>>& pParameter )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxCheckBoxProperty>( this, pParameter ) );

        return *this;
    } // method

    /* Add value property */
    mwxPropertyPanel& append( std::shared_ptr<AlignerParameter<double>>& pParameter )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxTextProperty<double>>( this, pParameter ) );

        return *this;
    } // method

    /* Add value property */
    mwxPropertyPanel& append( std::shared_ptr<AlignerParameter<short>>& pParameter )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxTextProperty<short>>( this, pParameter ) );

        return *this;
    } // method

    /* Add value property */
    mwxPropertyPanel& append( std::shared_ptr<AlignerParameter<uint64_t>>& pParameter )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxTextProperty<uint64_t>>( this, pParameter ) );

        return *this;
    } // method

    /* Add choices property */
    mwxPropertyPanel& append( std::shared_ptr<AlignerParameter<AlignerParameterBase::ChoicesType>>& pParameter )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxChoiceProperty>( this, pParameter ) );

        return *this;
    } // method

    /*
     * Add path/file (file-system) property
     * if a path shall be selected: rsSuffixWildcardString must be empty
     * if a file shall be selected: rsSuffixWildcardString must contain wildcard string as e.g.:
     *      "XYZ files (*.xyz)|*.xyz"
     */
    mwxPropertyPanel& append( std::shared_ptr<AlignerParameter<fs::path>>& pParameter,
                              const std::string& rsSuffixWildcardString = "" )
    { /* Add Property to vector*/
        if( rsSuffixWildcardString.empty( ) )
            vProperties.emplace_back( std::make_shared<mwxFolderProperty>( this, pParameter ) );
        else
            vProperties.emplace_back( std::make_shared<mwxFileProperty>( this, pParameter, rsSuffixWildcardString ) );

        return *this;
    } // method

    /* Add ... */
    mwxPropertyPanel& append( std::shared_ptr<AlignerParameterBase>& pParameter )
    { /* Add Property to vector*/
        {
            std::shared_ptr<AlignerParameter<int>> pCasted =
                std::dynamic_pointer_cast<AlignerParameter<int>>( pParameter );
            if( pCasted != nullptr )
                return append( pCasted );
        } // scope
        {
            auto pCasted = std::dynamic_pointer_cast<AlignerParameter<double>>( pParameter );
            if( pCasted != nullptr )
                return append( pCasted );
        } // scope
        {
            auto pCasted = std::dynamic_pointer_cast<AlignerParameter<bool>>( pParameter );
            if( pCasted != nullptr )
                return append( pCasted );
        } // scope
        {
            auto pCasted = std::dynamic_pointer_cast<AlignerParameter<short>>( pParameter );
            if( pCasted != nullptr )
                return append( pCasted );
        } // scope
        {
            auto pCasted = std::dynamic_pointer_cast<AlignerParameter<uint64_t>>( pParameter );
            if( pCasted != nullptr )
                return append( pCasted );
        } // scope
        {
            auto pCasted = std::dynamic_pointer_cast<AlignerParameter<AlignerParameterBase::ChoicesType>>( pParameter );
            if( pCasted != nullptr )
                return append( pCasted );
        } // scope
        {
            auto pCasted = std::dynamic_pointer_cast<AlignerParameter<fs::path>>( pParameter );
            if( pCasted != nullptr )
                return append( pCasted );
        } // scope

        throw std::runtime_error( std::string( "Could not append parameter:" )
                                      .append( pParameter->sName )
                                      .append( " to the settings menu. Maybe it's type is not overloaded?" ) );
    } // method

    void layout( )
    {
        // this->SetSizer( this->pxFlexLayout );
        this->Layout( );
        this->pxFlexLayout->Fit( this );
        // pxHostSizer->Add( this, 1, wxEXPAND | wxALL, 5 );
    } // method

    /* Updates enabled, disabled according to the corresponding parameter function */
    void updateEnabledDisabled( void )
    {
        for( auto pxProperty : vProperties )
        {
            pxProperty->updateEnabledDisabledDispatch( );
        } // for
    } // method
}; // class
