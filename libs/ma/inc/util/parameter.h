#pragma once

#include "util/default_parameters.h"

#if defined( __GNUC__ ) && ( __GNUC__ < 8 )
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/* Generic conversion of string to different value types.
 * (Specialized for some types)
 */
template <typename VALUE_TYPE> VALUE_TYPE genericStringToValue( const std::string& sString );

template <> std::string genericStringToValue<std::string>( const std::string& sString )
{
    return std::string( sString );
} // function

template <> int genericStringToValue<int>( const std::string& sString )
{
    return stoi( sString );
} // function

template <typename VALUE_TYPE> VALUE_TYPE genericStringToValue( const std::string& sString )
{
    std::stringstream xLineStream( sString );
    /* Throw an exception if something goes wrong with conversion */
    xLineStream.exceptions( std::ios::failbit );
    VALUE_TYPE value;
    xLineStream >> value;
    return value;
} // function


/* Predicate that checks for positive value */
static void checkPositiveValue( const int& iValue )
{
    if( iValue < 0 )
        throw std::out_of_range( "Positive values are allowed only." );
} // static method


/* Base class of all aligner parameter classes */
class AlignerParameterBase
{
  public:
    typedef std::vector<std::pair<std::string, std::string>> ChoicesType;

    const std::string sName; // Name of parameter
    const std::string sDescription; // Description of parameter

    AlignerParameterBase( const std::string& sName, const std::string& sDescription )
        : sName( sName ), sDescription( sDescription )
    {} // constructor
}; // class


/* Generic aligner parameter class.
 * Suitable for primitive types like int, double etc.
 */
template <typename VALUE_TYPE> class AlignerParameter : public AlignerParameterBase
{
  public:
    VALUE_TYPE value; // Current Value of parameter
    // Universal predicate
    static void predicateAlwaysOK( const VALUE_TYPE& )
    {} // static method
    std::function<void( const VALUE_TYPE& )> fPredicate;

    /* Constructor */
    AlignerParameter( const std::string& sName, const std::string& sDescription,
                      const VALUE_TYPE value, // initial value
                      std::function<void( const VALUE_TYPE& )> fPredicate = predicateAlwaysOK )
        : AlignerParameterBase( sName, sDescription ), value( value ), fPredicate( fPredicate )
    {} // constructor

    /* Throws an exception if something goes wrong */
    void set( VALUE_TYPE newValue )
    {
        /* Check the value. ( Error is indicated by exception )*/
        fPredicate( newValue );

        /* In the case of an exception this line is skipped !*/
        this->value = newValue;
    } // method

    /* Throws an exception if something goes wrong */
    void setByText( const std::string& sValueAsString )
    {
        VALUE_TYPE newValue = genericStringToValue<VALUE_TYPE>( sValueAsString );
        this->set( newValue );
    } // method

    VALUE_TYPE get( void )
    {
        return this->value;
    } // method

    void mirror( const AlignerParameter<VALUE_TYPE>& rxOther )
    {
        this->value = rxOther.value;
    } // method
}; // class


/* Parameter class for choices (list of different textual alternatives)
 */
template <> class AlignerParameter<AlignerParameterBase::ChoicesType> : public AlignerParameterBase
{
  public:
    AlignerParameterBase::ChoicesType vChoices; // Possible text values of the parameter
    unsigned int uiSelection; // integral value of the selected choice (0 for first choice, 1 for second etc.)

    AlignerParameter( const std::string& sName, const std::string& sDescription,
                      const AlignerParameterBase::ChoicesType& rvChoices, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameterBase( sName, sDescription ), vChoices( rvChoices ), uiSelection( uiSelection )
    {} // constructor

    /* Throws an exception if something goes wrong */
    void set( unsigned int uiNewSelection )
    {
        this->uiSelection = uiNewSelection;
    } // method

    /* Get the corresponding internal string for selection */
    std::string get( void )
    {
        std::cout << "Internal Setting:" << vChoices.at( uiSelection ).first << std::endl;
        return vChoices.at( uiSelection ).first;
    } // method

    void mirror( const AlignerParameter<AlignerParameterBase::ChoicesType>& rxOther )
    {
        this->uiSelection = rxOther.uiSelection;
    } // method
}; // class


/* Parameter class for file system paths */
template <> class AlignerParameter<fs::path> : public AlignerParameterBase
{
  public:
    fs::path xPath; // Current path

    AlignerParameter( const std::string& sName, const std::string& sDescription,
                      const fs::path& rxPath, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameterBase( sName, sDescription ), xPath( rxPath )
    {} // constructor

    /* Set path to argument.
     * (Throws an exception if something goes wrong.)
     */
    void set( const fs::path& rxPath )
    {
        this->xPath = rxPath;
    } // method

    void set( const std::string& rsPath )
    {
        this->xPath = rsPath;
    } // method

    /* Get the corresponding internal string for selection */
    fs::path get( void )
    {
        return this->xPath;
    } // method

    /* Used in the case of objects copies (object mirroring) */
    void mirror( const AlignerParameter<fs::path>& rxOther )
    {
        this->xPath = rxOther.xPath;
    } // method
}; // class


/* Class for management of aligner parameter sets.
 * Idea: By using Reflections this code could be improved. (e.g. https://www.rttr.org/)
 */
class Presetting
{
  public:
    /* IMPORTANT: If you add a new parameter, don't forget to add it to the mirror method */
    AlignerParameter<int> xMatch; // score for a DP match (used in SoC width computation)
    AlignerParameter<int> xMisMatch; // score for a DP match (used in SoC width computation)
    AlignerParameter<AlignerParameterBase::ChoicesType> xSeedingTechnique; // Seeding Technique
    AlignerParameter<bool> xUsePairedReads; // If true, work with paired reads

    AlignerParameter<bool> xExampleCheckBox;
    AlignerParameter<fs::path> xExamplePath;

    /* Delete copy constructor */
    Presetting( const Presetting& rxOtherSet ) = delete;
    Presetting( Presetting&& other ) = default; // required for getting emplace with maps working
    Presetting& operator=( Presetting&& other ) = default; // required for getting emplace with maps working

    /* Constructor */
    Presetting( )
        : // sName( sName ), //
          xMatch( "Match Score",
                  "Match score used in the context of Dynamic Programming\nand for SoC width computation.", 2,
                  checkPositiveValue ),
          xMisMatch( "Mismatch Penalty", "Penalty for a Dynamic Programming mismatch", 4, checkPositiveValue ),
          xSeedingTechnique( "Seeding Technique", "Technique used for the initial seeding.",
                             {{"maxSpan", "Maximally Spanning"}, {"SMEMs", "SMEMs"}} ),
          xUsePairedReads( "Use Paired Reads", "If your reads occur as paired reads, activate this flag.", false ),

          xExampleCheckBox( "Example Checkbox", "This is an example of a checkbox.", true ),
          xExamplePath( "Example Path", "This is an example of a file path.", "C:/" )
    {} // constructor

    /* Named copy Constructor */
    Presetting( const Presetting& rxOtherSet, const std::string& sName ) : Presetting( )
    {
        this->mirror( rxOtherSet );
    } // copy constructor

    /* Mirror the setting of one parameter-set into another */
    void mirror( const Presetting& rxOtherSet )
    {
        xMatch.mirror( rxOtherSet.xMatch );
        xMisMatch.mirror( rxOtherSet.xMisMatch );
        xSeedingTechnique.mirror( rxOtherSet.xSeedingTechnique );
        xUsePairedReads.mirror( rxOtherSet.xUsePairedReads );

        xExampleCheckBox.mirror( rxOtherSet.xExampleCheckBox );
        xExamplePath.mirror( rxOtherSet.xExamplePath );
    } // method

    // int EXPORTED iGap = 4; // penalty for a DP gap opening (used in SoC width computation)
    // int EXPORTED iExtend = 2; // penalty for a DP gap extension (used in SoC width computation)
    // int EXPORTED iGap2 = 24; // penalty for a DP gap opening (used in SoC width computation)
    // int EXPORTED iExtend2 = 1; // penalty for a DP gap extension (used in SoC width computation)

    /* Updates the global parameter using the values of the current parameter set */
    void updateGlobalParameter( void )
    {
        libMA::defaults::iMatch = this->xMatch.get( );
        libMA::defaults::iMissMatch = this->xMisMatch.get( );
        libMA::defaults::sSeedSet = this->xSeedingTechnique.get( );
    } // method

    /* True if the parameter-set uses paired reads */
    bool usesPairedReads( void )
    {
        return xUsePairedReads.value;
    } // method
}; // class


/* Parameter which occur only once.
 * (Not sequencing technique / configuration related parameters)
 */
class GeneralParameter
{
  public:
    // FIXME: Set this to a standard path in the beginning (~home/ma/SAM)
    AlignerParameter<bool> bSAMOutputInReadsFolder; // SAM Output in the same folder as the reads.
    AlignerParameter<fs::path> xSAMOutputPath; // folder path

    /* Constructor */
    GeneralParameter( )
        : bSAMOutputInReadsFolder( "SAM files in same folder as reads",
                                   "If set, the SAM files are written in the folder of the reads.", true ),
          xSAMOutputPath( "Folder (path) for SAM files", "All SAM-output will be written to this folder",
                          fs::temp_directory_path( ) )

    {} // constructor

    /* Named copy Constructor */
    GeneralParameter( const GeneralParameter& rxOtherSet ) : GeneralParameter( )
    {
        mirror( rxOtherSet );
    } // copy constructor

    /* Mirror the setting of the other parameter-set into the current parameter-set */
    void mirror( const GeneralParameter& rxOtherSet )
    {
        this->xSAMOutputPath.mirror( rxOtherSet.xSAMOutputPath );
        this->bSAMOutputInReadsFolder.mirror( rxOtherSet.bSAMOutputInReadsFolder );
    } // method
}; // class


/* Management of single global parameter-set and all config (preset) parameter-sets */
class ParameterSetManager
{
  private:
    // Pointer to the currently selected parameter-set
    Presetting* pSelectedParamSet = NULL;

  public:
    GeneralParameter xGlobalParameterSet;

    // Presets for aligner configuration
    std::map<std::string, Presetting> xParametersSets;

    ParameterSetManager( )
    {
        xParametersSets.emplace( "Illumina", Presetting( ) );
        xParametersSets.emplace( "Illumina Paired", Presetting( ) );
        xParametersSets[ "Illumina Paired" ].xUsePairedReads.set( true );

        // Initially select Illumina
        this->pSelectedParamSet = &( xParametersSets[ "Illumina" ] );
    } // constructor

    Presetting& get( const std::string& sKey )
    {
        return xParametersSets[ sKey ];
    } // method

    /* Set selected parameter set by using a key.
     * TODO: Check the key for existence.
     */
    void setSelected( const std::string& sKey )
    {
        this->pSelectedParamSet = &( xParametersSets[ sKey ] );
    } // method

    /* Delivers pointer to selected parameter-set */
    Presetting* getSelected( void )
    {
        return pSelectedParamSet;
    } // method
}; // class
