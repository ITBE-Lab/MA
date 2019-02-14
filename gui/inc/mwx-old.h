#pragma once

/* for all others, include the necessary headers (this file is usually all you
 * need because it includes almost all "standard" wxWidgets headers)
 */
#ifndef WX_PRECOMP
#include "wx/wx.h"
#include <wx/notebook.h>
#endif

#include <memory>
#include <vector>

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
 * to a menubar or to another menu will be deleted by their
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
        std::function<wxWindow*( mwxOK_Cancel_Dialog* )> // caller that makes the content of the dialog
            fMakeContent, 
        const wxPoint& xPos = wxDefaultPosition, //  position
        std::function<void( wxCommandEvent& )> OK_handler = []( wxCommandEvent& ) {} ) // position of the dialog
        : wxDialog( pxHostWindow, wxID_ANY, sTitle, xPos, wxSize( 500, 400 ),
                    ( wxDEFAULT_FRAME_STYLE | wxTAB_TRAVERSAL | wxSTAY_ON_TOP ) /* ^ wxRESIZE_BORDER */ ^ wxMAXIMIZE_BOX ^
                        wxMINIMIZE_BOX )
    {
        this->SetSizeHints( wxDefaultSize, wxDefaultSize );

        wxBoxSizer* pxTopSizer = new wxBoxSizer( wxVERTICAL );

        // Add inner content of dialog
        pxTopSizer->Add( fMakeContent( this ), wxSizerFlags( 1 ).Align( wxALIGN_CENTER ).Expand( ).Border( wxALL, 5 ) );

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

        // 'Fit' fits the dialog to the size of its children
        // SetSizerAndFit( pxTopSizer ); // use the sizer for layout and set size and hints

        // wxEvtHandler::Bind( wxEVT_BUTTON, OK_handler, iOK_Button_ID );
        this->SetSizer( pxTopSizer );
        this->Layout( );
        this->Centre( wxBOTH );
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
        auto* pxStaticBoxContext = new mwxStaticBoxContext( pxScrolledWindow, pxFlexGridSizer, "Basic Settings" );
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
        virtual ~mwxProperty( )
        {}
    }; // inner class

    /* Inner class describing text property */
    class mwxTextProperty : public mwxProperty
    {
      public:
        mwxPropertyPanel* pxHost; // Parent of property (manages destruction of text fields)
        wxTextCtrl* pxTextCtrl; // Text control holding the actual value

        mwxTextProperty( mwxPropertyPanel* pxHost, // parent of property
                         const std::string& sName, // property name
                         int* pValue // address of the managed value
                         )
            : pxHost( pxHost ),
              pxTextCtrl(
                  new wxTextCtrl( pxHost, wxID_ANY, std::to_string( *pValue ), wxDefaultPosition, wxDefaultSize, 0 ) )
        {
            /* Parent of host will manage destruction of wxStaticText and wxTextCtrl */
            pxHost->pxFlexLayout->Add(
                new wxStaticText( pxHost, wxID_ANY, sName.c_str( ), wxDefaultPosition, wxDefaultSize, 0 ), 0, wxALL,
                5 );

            pxHost->pxFlexLayout->Add( pxTextCtrl, 0, wxALL, 5 );
            pxTextCtrl->SetToolTip( "TextCtrl hint ..." );
        } // constructor
    }; // inner class

    /* Inner class describing choice property */
    class mwxChoiceProperty : public mwxProperty
    {
      public:
        mwxPropertyPanel* pxHost; // Parent of property (manages destruction of text fields)
        wxChoice* pxChoice; // Text control holding the actual value

        mwxChoiceProperty( mwxPropertyPanel* pxHost, // parent of property
                           const std::string& sName, // property name
                           const std::vector<std::string>& vsChoices, // property choices
                           int* pValue // address of the managed value corresponding to choice
                           )
            : pxHost( pxHost )
        {
            /* Parent of host will manage destruction of wxStaticText and wxTextCtrl */
            pxHost->pxFlexLayout->Add(
                new wxStaticText( pxHost, wxID_ANY, sName.c_str( ), wxDefaultPosition, wxDefaultSize, 0 ), 0, wxALL,
                5 );

            /* Prepare required string array of wxWidgets and create the choice widget */
            wxArrayString xChoices;
            for( auto& sChoiceText : vsChoices )
                xChoices.Add( sChoiceText.c_str( ) );
            pxChoice = new wxChoice( pxHost, wxID_ANY, wxDefaultPosition, wxDefaultSize, xChoices, 0 );
            pxHost->pxFlexLayout->Add( pxChoice, 0, wxALL, 5 );
        } // constructor
    }; // inner class

    /* Inner class describing choice property */
    class mwxCheckBoxProperty : public mwxProperty
    {
      public:
        mwxPropertyPanel* pxHost; // Parent of property (manages destruction of text fields)
        wxCheckBox* pxCheckBox; // Text control holding the actual value

        mwxCheckBoxProperty( mwxPropertyPanel* pxHost, // parent of property
                             const std::string& sName, // property name
                             int* pValue // address of the managed value corresponding to choice
                             )
            : pxHost( pxHost ),
              pxCheckBox( new wxCheckBox( pxHost, wxID_ANY, wxT( "" ), wxDefaultPosition, wxDefaultSize, 0 ) )
        {
            /* Parent of host will manage destruction of wxStaticText and wxTextCtrl */
            pxHost->pxFlexLayout->Add(
                new wxStaticText( pxHost, wxID_ANY, sName.c_str( ), wxDefaultPosition, wxDefaultSize, 0 ), 0, wxALL,
                5 );

            pxHost->pxFlexLayout->Add( pxCheckBox, 0, wxALL, 5 );
            pxCheckBox->SetToolTip( "Checkbox hint ..." );
        } // constructor
    }; // inner class

  public:
    wxWindow* pxHostWindow; // Host window
    wxSizer* pxHostSizer; // Host sizer
    wxFlexGridSizer* pxFlexLayout; // Layout of property editor
    std::vector<std::shared_ptr<mwxProperty>> vProperties; // Vector of all properties

    /* Constructor requires parent widgets */
    mwxPropertyPanel( wxWindow* pxHostWindow, // host window of box context (responsible for destruction)
                      wxSizer* pxHostSizer ) // sizer of the host window
        : wxPanel( pxHostWindow, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                   wxTAB_TRAVERSAL ), // base class constructor call
          pxHostWindow( pxHostWindow ),
          pxHostSizer( pxHostSizer ),
          pxFlexLayout( new wxFlexGridSizer( 0, 2, 0, 0 ) )
    {
        this->pxFlexLayout->SetFlexibleDirection( wxBOTH ); // acoording to wxFormBuild
        this->pxFlexLayout->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED ); // acoording to wxFormBuild
        this->SetSizer( this->pxFlexLayout );
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
        vProperties.emplace_back( std::make_shared<mwxTextProperty>( this, sName, pValue ) );

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

    /* Add checkbox property */
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
