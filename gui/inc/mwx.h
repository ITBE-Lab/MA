#pragma once

/* for all others, include the necessary headers (this file is usually all you
 * need because it includes almost all "standard" wxWidgets headers)
 */
#ifndef WX_PRECOMP
#include "wx/wx.h"
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
        std::cout << "In ~mwxMenu()" << std::endl;
    } // destructor

    /* All menu-items of the menu */
    std::vector<mwxMenuItem> vChilds;

    /* Appends an item to the menu. */
    mwxMenu& append(
        const std::string& sMenuItemText, // text of menu-item
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
        std::cout << "In ~mwxMenuBar()" << std::endl;
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
    /* Constructor */
    mwxBitmapButton(
        wxWindow* pxHostWindow, // host window of box context (responsible for destruction)
        const wxBitmap& rxBitmapNormal, // Bitmap unpressed
        const wxBitmap& rxBitmapPressed, // Bitmap for pressed button
        std::function<void( wxCommandEvent& )> handler = []( wxCommandEvent& ) {}, // handler
        const wxPoint rxPos = wxDefaultPosition,
        const wxSize rxSize = wxDefaultSize )
        : wxBitmapButton( pxHostWindow, wxID_ANY, rxBitmapNormal, rxPos, rxSize, wxNO_BORDER )
    {
        this->SetBitmapPressed( rxBitmapPressed );
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
    mwxOK_Cancel_Dialog(
        wxWindow* pxHostWindow, // host window of box context (responsible for destruction)
        const wxString& sTitle, // title of the dialog
        std::function<wxWindow*( wxWindow* )>
            fMakeContent, // called for creating the content of the dialog
        const wxPoint& xPos = wxDefaultPosition, // position of the dialog
        bool bDoFit = false, // adapt the size of the dialog to the size of its content
        std::function<void( wxCommandEvent& )> OK_handler = []( wxCommandEvent& ) {} ) // handler of OK button

        : wxDialog( pxHostWindow, wxID_ANY, sTitle, xPos, wxSize( 500, 400 ),
                    ( wxDEFAULT_FRAME_STYLE | wxTAB_TRAVERSAL ) /* ^ wxRESIZE_BORDER */ ^
                        wxMAXIMIZE_BOX ^ wxMINIMIZE_BOX ^ ( bDoFit ? wxRESIZE_BORDER : 0 ) )
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
        pxButtomSizer->Add( new wxButton( this, wxID_OK, "OK" ),
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
  private:
    const std::string sFileDialogTitle; // "Open Reference file"
    const std::string sFilePattern; // "Reference files (*.maReference)|*.maReference|All Files|*"
    std::function<void( const std::string& )> fHandler;

    void onFolderButton( wxCommandEvent& WXUNUSED( event ) )
    {
        wxFileDialog openFileDialog( pxConnector.pxWindow, sFileDialogTitle, "", "", sFilePattern,
                                     wxFD_OPEN | wxFD_FILE_MUST_EXIST );
        if( openFileDialog.ShowModal( ) == wxID_OK )
        {
            // Call handler to inform about new selection
            fHandler( std::string( openFileDialog.GetPath( ).c_str( ) ) );
        } // if
    } // method

    /* Event handler for clear button */
    void onClearButton( wxCommandEvent& WXUNUSED( event ) )
    {
        // Clear Button results in call of handler with empty string
        fHandler( std::string( "" ) );
    } // method

  public:
    mwxFileSelectDeleteButtonSizer(
        const mwxConnector& pxConnector, // connector for BoxSizer
        const std::string& rsText, // text of the button group
        const std::string& rsFileDialogTitle,
        const std::string& rsFilePattern,
        std::function<void( const std::string& )> fHandler = []( const std::string& ) {}, // called, if OK is selected
        bool bWithClearButton = true )
        : mwxBoxSizer( pxConnector, wxVERTICAL ),
          pxConnector( pxConnector ),
          sFileDialogTitle( rsFileDialogTitle ),
          sFilePattern( rsFilePattern ),
          fHandler( fHandler )
    {
        Add( new wxStaticText( pxConnector.pxWindow, wxID_ANY, rsText, wxDefaultPosition, wxDefaultSize, 0 ), 0,
             wxALL | wxLEFT, 5 );
        addBoxSizer( // horizontal BoxSizer (File Selection and Deletion)
            wxHORIZONTAL,
            wxSizerFlags( 0 ).ReserveSpaceEvenIfHidden( ),
            [&]( mwxBoxSizer& pxBoxSizer ) // BoxSizer content
            {
                pxBoxSizer.Add( new mwxBitmapButton( pxBoxSizer.pxConnector.pxWindow,
                                                     wxBITMAP_PNG_FROM_DATA( FolderButtonNormal ),
                                                     wxBITMAP_PNG_FROM_DATA( FolderButtonSelected ),
                                                     std::bind( &mwxFileSelectDeleteButtonSizer::onFolderButton, this,
                                                                std::placeholders::_1 ) ),
                                0, wxRESERVE_SPACE_EVEN_IF_HIDDEN | wxRIGHT | wxLEFT, 5 ); // GTK Bug Border 0

                mwxBitmapButton* pxClearButton;
                pxBoxSizer.Add( pxClearButton =
                                    new mwxBitmapButton( pxBoxSizer.pxConnector.pxWindow,
                                                         wxBITMAP_PNG_FROM_DATA( DeleteButtonNormal ),
                                                         wxBITMAP_PNG_FROM_DATA( DeleteButtonSelected ),
                                                         std::bind( &mwxFileSelectDeleteButtonSizer::onClearButton,
                                                                    this, std::placeholders::_1 ) ),
                                0, wxRESERVE_SPACE_EVEN_IF_HIDDEN | wxRIGHT | wxLEFT, 5 ); // GTK Bug Border 0
                pxClearButton->Show( bWithClearButton );
            } ); // addBoxSizer

        fHandler( std::string( "" ) );
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
        std::function<void( wxCommandEvent& )> handler = []( wxCommandEvent& ) {} ) // handler
        : wxComboBox( )
    {
        for( auto& rxKeyValue : rxMap )
        {
            xChoices.Add( rxKeyValue.first.c_str( ) );
        } // for

        this->Create( pxHost, wxID_ANY, xChoices[ 0 ], wxDefaultPosition, wxDefaultSize, xChoices, wxCB_READONLY );
        // Bind handler that process changes
        this->wxEvtHandler::Bind( wxEVT_COMBOBOX, handler );
    } // constructor
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
        std::cout << "In ~mwxStaticBoxContext()" << std::endl;
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

        /* Constructor */
        mwxProperty( mwxPropertyPanel* pxHost, // parent of property
                     const std::string& sName ) // name of property
            : pxHost( pxHost )
        {
            /* Parent of host will manage destruction of wxStaticText */
            pxHost->pxFlexLayout->Add(
                new wxStaticText( pxHost, wxID_ANY, sName.c_str( ), wxDefaultPosition, wxDefaultSize, 0 ), 0,
                wxALL | wxALIGN_CENTER_VERTICAL, 5 );
        } // constructor

        virtual ~mwxProperty( )
        {}
    }; // inner class

    /* Inner class describing text property */
    template <typename VALUE_TYPE> class mwxTextProperty : public mwxProperty
    {
      public:
        wxTextCtrl* pxTextCtrl; // Text control holding the actual value
        AlignerParameter<VALUE_TYPE>* pxParameter = nullptr; // pointer to parameter object

        /* Update the value of the property */
        void update( void )
        {
            pxTextCtrl->SetValue( std::to_string( this->pxParameter->value ) );
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
        } // method

        /* Constructor */
        mwxTextProperty( mwxPropertyPanel* pxHost, // parent of property
                         const std::string& sName, // property name
                         int* pValue // address of the managed value
                         )
            : mwxProperty( pxHost, sName ),
              pxTextCtrl( new wxTextCtrl( pxHost, wxID_ANY, std::to_string( *pValue ), wxDefaultPosition, wxDefaultSize,
                                          wxTE_PROCESS_ENTER ) )
        {
            pxHost->pxFlexLayout->Add( pxTextCtrl, 0, wxALL, 5 );

            /* Make the enter key to behave like the tab key. */
            pxTextCtrl->wxEvtHandler::Bind( wxEVT_COMMAND_TEXT_ENTER,
                                            [&]( wxCommandEvent& rxEvent ) { this->pxTextCtrl->Navigate( ); } );
        } // constructor

        /* Constructor with AlignerParameter */
        mwxTextProperty( mwxPropertyPanel* pxHost, // parent of property
                         AlignerParameter<VALUE_TYPE>* pxParameter // parameter reference
                         )
            : mwxProperty( pxHost, pxParameter->sName ),
              pxTextCtrl( new wxTextCtrl( pxHost, wxID_ANY, std::to_string( pxParameter->value ), wxDefaultPosition,
                                          wxDefaultSize, wxTE_PROCESS_ENTER ) ),
              pxParameter( pxParameter )
        {
            pxHost->pxFlexLayout->Add( pxTextCtrl, 0, wxALL, 5 );

            /* The description of the parameter is used as ToolTip-text */
            pxTextCtrl->SetToolTip( pxParameter->sDescription );

            /* Bind "focus leave" to the handler */
            pxTextCtrl->wxEvtHandler::Bind( wxEVT_KILL_FOCUS,
                                            [this]( wxFocusEvent& rxEvent ) { this->focusHandler( rxEvent ); } );
            /* Make the enter key to behave like the tab key. */
            pxTextCtrl->wxEvtHandler::Bind( wxEVT_COMMAND_TEXT_ENTER,
                                            [&]( wxCommandEvent& rxEvent ) { this->pxTextCtrl->Navigate( ); } );
        } // constructor
    }; // inner class

    /* Inner class describing text property */
    class mwxFolderProperty : public mwxProperty
    {
      public:
        mwxPropertyPanel* pxHost; // keep host for folder dialog creation
        AlignerParameter<fs::path>* pxParameter = nullptr; // pointer to parameter object
        wxTextCtrl* pxTextCtrl; // Text control holding the actual value
        mwxBitmapButton* pxButton; // Button for folder selection


        void onButton( wxCommandEvent& WXUNUSED( event ) )
        {
            wxDirDialog xDirDialog( pxHost, "Select folder", "" );
            if( xDirDialog.ShowModal( ) == wxID_OK )
            {
                pxTextCtrl->SetValue( std::string( xDirDialog.GetPath( ) ) );
                this->pxParameter->set( std::string( xDirDialog.GetPath( ) ) );
            } // if
        } // method

        /* Update the value of the property */
        void update( void )
        {
            pxTextCtrl->SetValue( this->pxParameter->get( ).string( ) );
        } // method

        /* Handler is called whenever the user leaves the input field (the user ends property editing) */
        void focusHandler( wxFocusEvent& rxEvent )
        {
            // Wrong inputs (out of range etc.) are reported via exceptions
            std::string sPath( this->pxTextCtrl->GetValue( ).c_str( ) );

            // Check folder in textCtrl for existence
            if( fs::exists( fs::path( sPath ) ) )
                pxParameter->set( sPath );
            else
                wxMessageBox( wxString( std::string( "Folder " ) + sPath + " does not exist." ), "Error",
                              wxICON_ERROR );
            this->update( );
            // From wxWidgets documentation:
            // The focus event handlers should almost invariably call
            // wxEvent::Skip() on their event argument to allow the default handling to take place.
            rxEvent.Skip( );
        } // method

        /* Constructor for AlignerParameter */
        mwxFolderProperty( mwxPropertyPanel* pxHost, // parent of property
                           AlignerParameter<fs::path>* pxParameter // parameter reference
                           )
            : mwxProperty( pxHost, pxParameter->sName ),
              pxHost( pxHost ),
              pxParameter( pxParameter ),
              pxTextCtrl( new wxTextCtrl( pxHost, wxID_ANY, this->pxParameter->get( ).string( ), wxDefaultPosition,
                                          wxSize( 200, -1 ), wxTE_PROCESS_ENTER ) ),
              pxButton( new mwxBitmapButton( pxHost, wxBITMAP_PNG_FROM_DATA( FolderButtonSmall ),
                                             wxBITMAP_PNG_FROM_DATA( FolderButtonSmallSelected ),
                                             std::bind( &mwxFolderProperty::onButton, this, std::placeholders::_1 ) ) )

        {
            // Place textCtrl plus button in horizontal sizer
            auto* pxSizer = new wxBoxSizer( wxHORIZONTAL );
            pxSizer->Add( pxTextCtrl, 0, wxALL | wxALIGN_CENTER_VERTICAL, 5 );
            pxSizer->Add( pxButton, 0, wxALL | wxALIGN_CENTER_VERTICAL, 5 );
            pxHost->pxFlexLayout->Add( pxSizer );

            // // Bind button to folder selection
            // this->pxButton->wxEvtHandler::Bind(
            //     wxEVT_BUTTON, std::bind( &mwxFolderProperty::onButton, this, std::placeholders::_1 ), wxID_ANY );

            // The description of the parameter is used as ToolTip-text
            this->pxTextCtrl->SetToolTip( pxParameter->sDescription );
            this->pxButton->SetToolTip( "Folder selection" );

            // Bind "focus leave" to the handler
            this->pxTextCtrl->wxEvtHandler::Bind( wxEVT_KILL_FOCUS,
                                                  [this]( wxFocusEvent& rxEvent ) { this->focusHandler( rxEvent ); } );
            // Make the enter key to behave like the tab key
            this->pxTextCtrl->wxEvtHandler::Bind( wxEVT_COMMAND_TEXT_ENTER,
                                                  [&]( wxCommandEvent& rxEvent ) { this->pxTextCtrl->Navigate( ); } );
        } // constructor
    }; // inner class


    /* Inner class describing choice property */
    class mwxChoiceProperty : public mwxProperty
    {
      public:
        AlignerParameter<AlignerParameterBase::ChoicesType>* pxParameter = nullptr; // pointer to parameter object
        wxComboBox* pxCombo; // Text control holding the actual value

        /* Update the value of the property */
        void update( void )
        {
            pxCombo->SetSelection( pxParameter->uiSelection );
        } // method

        /* Handler is called whenever the user leaves the input field (the user ends property editing) */
        void focusHandler( wxFocusEvent& rxEvent )
        {
            // Wrong input (out of range etc.) are reported via exceptions
            try
            {
                std::cout << "Selected is: " << pxCombo->GetSelection( ) << std::endl; // GetValue
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
            rxEvent.Skip( );
        } // method

        /* Constructor */
        mwxChoiceProperty( mwxPropertyPanel* pxHost, // parent of property
                           const std::string& sName, // property name
                           const std::vector<std::string>& vsChoices, // property choices
                           int* pValue // address of the managed value corresponding to choice
                           )
            : mwxProperty( pxHost, sName )
        {
            // Prepare required string array of wxWidgets and create the choice widget
            wxArrayString xChoices;
            for( auto& sChoiceText : vsChoices )
                xChoices.Add( sChoiceText.c_str( ) );
            pxCombo = new wxComboBox( pxHost, wxID_ANY, xChoices[ 0 ], wxDefaultPosition, wxDefaultSize, xChoices );
            pxHost->pxFlexLayout->Add( pxCombo, 0, wxALL, 5 );
        } // constructor

        /* Constructor with AlignerParameter */
        mwxChoiceProperty( mwxPropertyPanel* pxHost, // parent of property
                           AlignerParameter<AlignerParameterBase::ChoicesType>* pxParameter // parameter reference
                           )
            : mwxProperty( pxHost, pxParameter->sName ), pxParameter( pxParameter )
        {
            /* Prepare required string array of wxWidgets and create the choice widget */
            wxArrayString aChoices;
            for( auto& sChoiceText : pxParameter->vChoices )
                aChoices.Add( sChoiceText.second.c_str( ) );
            pxCombo = new wxComboBox( pxHost, wxID_ANY, aChoices[ pxParameter->uiSelection ], wxDefaultPosition,
                                      wxDefaultSize, aChoices, wxCB_READONLY );
            pxHost->pxFlexLayout->Add( pxCombo, 0, wxALL, 5 );

            // Description of the parameter is used as ToolTip-text
            pxCombo->SetToolTip( pxParameter->sDescription );

            // Bind handler that process changes
            pxCombo->wxEvtHandler::Bind( wxEVT_KILL_FOCUS,
                                         [this]( wxFocusEvent& rxEvent ) { this->focusHandler( rxEvent ); } );
        } // constructor
    }; // inner class


    /* Inner class describing choice property */
    class mwxCheckBoxProperty : public mwxProperty
    {
      public:
        AlignerParameter<bool>* pxParameter = nullptr; // pointer to parameter object
        wxCheckBox* pxCheckBox; // Text control holding the actual value

        /* Handler is called whenever the user leaves the input field (the user ends property editing) */
        void commandHandler( wxCommandEvent& rxEvent )
        {
            // Copy boolean value to parameter
            pxParameter->set( pxCheckBox->GetValue( ) );
        } // method

        mwxCheckBoxProperty( mwxPropertyPanel* pxHost, // parent of property
                             const std::string& sName, // property name
                             int* pValue // address of the managed value corresponding to choice
                             )
            : mwxProperty( pxHost, sName ), // pxHost( pxHost ),
              pxCheckBox( new wxCheckBox( pxHost, wxID_ANY, wxT( "" ), wxDefaultPosition, wxDefaultSize, 0 ) )
        {
            pxHost->pxFlexLayout->Add( pxCheckBox, 0, wxALL, 5 );
        } // constructor

        mwxCheckBoxProperty( mwxPropertyPanel* pxHost, // parent of property
                             AlignerParameter<bool>* pxParameter // parameter reference
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

            // Bind handler that process changes
            pxCheckBox->wxEvtHandler::Bind( wxEVT_CHECKBOX,
                                            [this]( wxCommandEvent& rxEvent ) { this->commandHandler( rxEvent ); } );
        } // constructor
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
    mwxPropertyPanel& append( const std::string& sName, // property name
                              int* pValue // address of the managed value
    )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxTextProperty<int>>( this, sName, pValue ) );

        return *this;
    } // method

    /* Add value property */
    mwxPropertyPanel& append( AlignerParameter<int>& rxParameter )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxTextProperty<int>>( this, &rxParameter ) );

        return *this;
    } // method

    /* Add value property */
    mwxPropertyPanel& append( AlignerParameter<bool>& rxParameter )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxCheckBoxProperty>( this, &rxParameter ) );

        return *this;
    } // method

    /* Add choices property */
    mwxPropertyPanel& append( AlignerParameter<AlignerParameterBase::ChoicesType>& rxParameter )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxChoiceProperty>( this, &rxParameter ) );

        return *this;
    } // method

    /* Add path (file-system) property */
    mwxPropertyPanel& append( AlignerParameter<fs::path>& rxParameter )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxFolderProperty>( this, &rxParameter ) );

        return *this;
    } // method

    /* Add choice property */
    mwxPropertyPanel& append( const std::string& sName, // property name
                              const std::vector<std::string>& vsChoices, // choices
                              int* pValue // address of the managed value
    )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxChoiceProperty>( this, sName, vsChoices, pValue ) );
        return *this;
    } // method

    /* Add check-box property */
    mwxPropertyPanel& appendCheckBox( const std::string& sName, // property name
                                      int* pValue // address of the managed value
    )
    { /* Add Property to vector*/
        vProperties.emplace_back( std::make_shared<mwxCheckBoxProperty>( this, sName, pValue ) );
        return *this;
    } // method

    void layout( )
    {
        // this->SetSizer( this->pxFlexLayout );
        this->Layout( );
        this->pxFlexLayout->Fit( this );
        // pxHostSizer->Add( this, 1, wxEXPAND | wxALL, 5 );
    } // method
}; // class
