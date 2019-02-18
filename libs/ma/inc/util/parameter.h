#pragma once

#include "util/default_parameters.h"
#include "util/exception.h"
#include "util/exported.h"

#if defined( __GNUC__ ) && ( __GNUC__ < 8 )
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class Presetting;

/* Generic conversion of string to different value types.
 * (Specialized for some types)
 */
template <typename VALUE_TYPE> EXPORTED VALUE_TYPE genericStringToValue( const std::string& sString );

// template <typename VALUE_TYPE> VALUE_TYPE genericStringToValue( const std::string& sString )
// {
//     std::stringstream xLineStream( sString );
//     /* Throw an exception if something goes wrong with conversion */
//     xLineStream.exceptions( std::ios::failbit );
//     VALUE_TYPE value;
//     xLineStream >> value;
//     return value;
// } // function

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
    static const char NO_SHORT_DEFINED = ' ';

    const std::string sName; // Name of parameter
    const std::string sDescription; // Description of parameter
    const char cShort;

    EXPORTED AlignerParameterBase( Presetting* pPresetting, const std::string& sName, const char cShort,
                                   const std::string& sDescription );

    void mirror( const AlignerParameterBase* pOther )
    {
        throw AnnotatedException( "trying to mirror pointers of different types" );
    } // method
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
    AlignerParameter( Presetting* pPresetting, const std::string& sName, const char cShort,
                      const std::string& sDescription,
                      const VALUE_TYPE value, // initial value
                      std::function<void( const VALUE_TYPE& )> fPredicate = predicateAlwaysOK )
        : AlignerParameterBase( pPresetting, sName, cShort, sDescription ), value( value ), fPredicate( fPredicate )
    {} // constructor

    /* Constructor */
    AlignerParameter( Presetting* pPresetting, const std::string& sName, const std::string& sDescription,
                      const VALUE_TYPE value, // initial value
                      std::function<void( const VALUE_TYPE& )> fPredicate = predicateAlwaysOK )
        : AlignerParameter( pPresetting, sName, NO_SHORT_DEFINED, sDescription, value, predicateAlwaysOK )
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

    const VALUE_TYPE get( void ) const
    {
        return this->value;
    } // method

    void mirror( const AlignerParameter<VALUE_TYPE>* pOther )
    {
        this->value = pOther->value;
    } // method
}; // class


/* Parameter class for choices (list of different textual alternatives)
 */
template <> class AlignerParameter<AlignerParameterBase::ChoicesType> : public AlignerParameterBase
{
  public:
    AlignerParameterBase::ChoicesType vChoices; // Possible text values of the parameter
    unsigned int uiSelection; // integral value of the selected choice (0 for first choice, 1 for second etc.)

    AlignerParameter( Presetting* pPresetting, const std::string& sName, const char cShort,
                      const std::string& sDescription,
                      const AlignerParameterBase::ChoicesType& rvChoices, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameterBase( pPresetting, sName, cShort, sDescription ),
          vChoices( rvChoices ),
          uiSelection( uiSelection )
    {} // constructor

    AlignerParameter( Presetting* pPresetting, const std::string& sName, const std::string& sDescription,
                      const AlignerParameterBase::ChoicesType& rvChoices, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameter( pPresetting, sName, NO_SHORT_DEFINED, sDescription, rvChoices, uiSelection )
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

    /* Get the corresponding internal string for selection */
    const std::string get( void ) const
    {
        std::cout << "Internal Setting:" << vChoices.at( uiSelection ).first << std::endl;
        return vChoices.at( uiSelection ).first;
    } // method

    void mirror( const AlignerParameter<AlignerParameterBase::ChoicesType>* pOther )
    {
        this->uiSelection = pOther->uiSelection;
    } // method
}; // class


/* Parameter class for file system paths */
template <> class AlignerParameter<fs::path> : public AlignerParameterBase
{
  public:
    fs::path xPath; // Current path

    AlignerParameter( Presetting* pPresetting, const std::string& sName, const char cShort,
                      const std::string& sDescription,
                      const fs::path& rxPath, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameterBase( pPresetting, sName, cShort, sDescription ), xPath( rxPath )
    {} // constructor

    /* Constructor */
    AlignerParameter( Presetting* pPresetting, const std::string& sName, const std::string& sDescription,
                      const fs::path& rxPath, // choices of the parameter
                      unsigned int uiSelection = 0 )
        : AlignerParameter( pPresetting, sName, NO_SHORT_DEFINED, sDescription, rxPath, uiSelection )
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

    /* Get the corresponding internal string for selection */
    const fs::path get( void ) const
    {
        return this->xPath;
    } // method

    /* Used in the case of objects copies (object mirroring) */
    void mirror( const AlignerParameter<fs::path>* pOther )
    {
        this->xPath = pOther->xPath;
    } // method
}; // class


/* Class for management of aligner parameter sets.
 * Idea: By using Reflections this code could be improved. (e.g. https://www.rttr.org/)
 */
class Presetting
{
  public:
    std::map<std::string, AlignerParameterBase*> xpAllParameters;
    std::map<char, AlignerParameterBase*> xpParametersByShort;
    AlignerParameter<int> xMatch; // score for a DP match (used in SoC width computation)
    AlignerParameter<int> xMisMatch; // score for a DP match (used in SoC width computation)
    AlignerParameter<int> xGap; // @todo
    AlignerParameter<int> xExtend; // @todo
    AlignerParameter<int> xGap2; // @todo
    AlignerParameter<int> xExtend2; // @todo
    AlignerParameter<double> xPairedBonus; // @todo
    AlignerParameter<double> xMeanPairedReadDistance; // @todo
    AlignerParameter<double> xStdPairedReadDistance; // @todo
    AlignerParameter<int> xReportN; // @todo
    AlignerParameter<int> xMaximalSeedAmbiguity; // @todo
    AlignerParameter<int> xMinSeedLength; // @todo
    AlignerParameter<int> xMinAlignmentScore; // @todo
    AlignerParameter<int> xMinimalSeedAmbiguity; // @todo
    AlignerParameter<int> xMinimalSeedSizeDrop; // @todo
    AlignerParameter<int> xMaxNumSoC; // @todo
    AlignerParameter<int> xMinNumSoC; // @todo
    AlignerParameter<int> xMaxScoreLookahead; // @todo
    AlignerParameter<int> xSwitchQlen; // @todo
    AlignerParameter<int> xMaxGapArea; // @todo
    AlignerParameter<int> xPadding; // @todo
    AlignerParameter<int> xSoCWidth; // @todo
    AlignerParameter<bool> xOptimisticGapCostEstimation; // @todo
    AlignerParameter<bool> xSkipAmbiguousSeeds; // @todo
    AlignerParameter<double> xRelMinSeedSizeAmount; // @todo
    AlignerParameter<double> xScoreDiffTolerance; // @todo
    AlignerParameter<double> xMinQueryCoverage; // @todo
    AlignerParameter<double> xSoCScoreDecreaseTolerance; // @todo
    AlignerParameter<int> xHarmScoreMin; // @todo
    AlignerParameter<double> xHarmScoreMinRel; // @todo
    AlignerParameter<int> xGenomeSizeDisable; // @todo
    AlignerParameter<bool> xDisableHeuristics; // @todo
    AlignerParameter<bool> xNoSecondary; // @todo
    AlignerParameter<bool> xNoSupplementary; // @todo
    AlignerParameter<bool> xDisableGapCostEstimationCutting; // @todo
    AlignerParameter<double> xMaxDeltaDist; // @todo
    AlignerParameter<int> xMinDeltaDist; // @todo
    AlignerParameter<double> xMaxOverlapSupplementary; // @todo
    AlignerParameter<int> xMaxSupplementaryPerPrim; // @todo
    AlignerParameter<int> xSVPenalty; // @todo
    AlignerParameter<int> xMinBandwidthGapFilling; // @todo
    AlignerParameter<int> xBandwidthDPExtension; // @todo
    AlignerParameter<double> xMaxSVRatio; // @todo
    AlignerParameter<int> xMinSVDistance; // @todo
    AlignerParameter<int> xZDrop; // @todo

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
          xMatch( this, "Match Score",
                  "Match score used in the context of Dynamic Programming\nand for SoC width computation.", 2,
                  checkPositiveValue ),
          xMisMatch( this, "Mismatch Penalty", "Penalty for a Dynamic Programming mismatch", 4, checkPositiveValue ),
          xGap( this, "", "", 4 ),
          xExtend( this, "", "", 2 ),
          xGap2( this, "", "", 24 ),
          xExtend2( this, "", "", 1 ),
          xPairedBonus( this, "", "", 1.25 ),
          xMeanPairedReadDistance( this, "", "", 400 ),
          xStdPairedReadDistance( this, "", "", 150 ),
          xReportN( this, "", "", 0 ),
          xMaximalSeedAmbiguity( this, "", "", 500 ),
          xMinSeedLength( this, "", "", 16 ),
          xMinAlignmentScore( this, "", "", 75 ),
          xMinimalSeedAmbiguity( this, "", "", 0 ),
          xMinimalSeedSizeDrop( this, "", "", 15 ),
          xMaxNumSoC( this, "", "", 30 ),
          xMinNumSoC( this, "", "", 1 ),
          xMaxScoreLookahead( this, "", "", 3 ),
          xSwitchQlen( this, "", "", 800 ),
          xMaxGapArea( this, "", "", 10000 ),
          xPadding( this, "", "", 1000 ),
          xSoCWidth( this, "", "", 0 ),
          xOptimisticGapCostEstimation( this, "", "", true ),
          xSkipAmbiguousSeeds( this, "", "", false ),
          xRelMinSeedSizeAmount( this, "", "", 0.005 ),
          xScoreDiffTolerance( this, "", "", 0.0001 ),
          xMinQueryCoverage( this, "", "", 1.1 ),
          xSoCScoreDecreaseTolerance( this, "", "", 0.1 ),
          xHarmScoreMin( this, "", "", 18 ),
          xHarmScoreMinRel( this, "", "", 0.002 ),
          xGenomeSizeDisable( this, "", "", 10000000 ),
          xDisableHeuristics( this, "", "", false ),
          xNoSecondary( this, "", "", false ),
          xNoSupplementary( this, "", "", false ),
          xDisableGapCostEstimationCutting( this, "", "", false ),
          xMaxDeltaDist( this, "", "", 0.1 ),
          xMinDeltaDist( this, "", "", 16 ),
          xMaxOverlapSupplementary( this, "", "", 0.1 ),
          xMaxSupplementaryPerPrim( this, "", "", 1 ),
          xSVPenalty( this, "", "", 100 ),
          xMinBandwidthGapFilling( this, "", "", 20 ),
          xBandwidthDPExtension( this, "", "", 512 ),
          xMaxSVRatio( this, "", "", 0.01 ),
          xMinSVDistance( this, "", "", 500 ),
          xZDrop( this, "", "", 200 ),

          xSeedingTechnique( this, "Seeding Technique", "Technique used for the initial seeding.",
                             {{"maxSpan", "Maximally Spanning"}, {"SMEMs", "SMEMs"}} ),
          xUsePairedReads( this, "Use Paired Reads", "If your reads occur as paired reads, activate this flag.",
                           false ),
          xExampleCheckBox( this, "Example Checkbox", "This is an example of a checkbox.", true ),
          xExamplePath( this, "Example Path", "This is an example of a file path.", "C:/" )
    {} // constructor

    /* Mirror the setting of one parameter-set into another */
    void mirror( const Presetting& rxOtherSet )
    {
        for( auto& rxTup : xpAllParameters )
            rxTup.second->mirror( rxOtherSet.xpAllParameters.at( rxTup.first ) );
    } // method

    /* Named copy Constructor */
    Presetting( const Presetting& rxOtherSet, const std::string& sName ) : Presetting( )
    {
        this->mirror( rxOtherSet );
    } // copy constructor

    /* True if the parameter-set uses paired reads */
    bool usesPairedReads( void )
    {
        return xUsePairedReads.value;
    } // method

    /**
     * Every parameter has to call this function from it's constructor.
     * We use this in order to generate a map of all available parameters
     */
    void registerParameter( AlignerParameterBase* pParameter )
    {
        xpAllParameters.emplace( pParameter->sName, pParameter );
        if( pParameter->cShort != pParameter->NO_SHORT_DEFINED )
        {
            auto xEmplaceRet = xpParametersByShort.emplace( pParameter->cShort, pParameter );
            if( !xEmplaceRet.second ) // there is another parameter with this shorthand...
                throw AnnotatedException( std::string( "Shorthand: " )
                                              .append( std::to_string( pParameter->cShort ) )
                                              .append( " used twice: 1) " )
                                              .append( xEmplaceRet.first->second->sName )
                                              .append( " 2)" )
                                              .append( pParameter->sName ) );
        }
    } // method

    AlignerParameterBase* byName( const std::string& rS )
    {
        try{
            return xpAllParameters.at( rS );
        } catch(std::out_of_range& rEx)
        {
            throw AlignerException(  );
        }
    } // method

    AlignerParameterBase* byShort( const char cX )
    {
        return xpParametersByShort.at( cX );
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
        : bSAMOutputInReadsFolder( nullptr, "SAM files in same folder as reads",
                                   "If set, the SAM files are written in the folder of the reads.", true ),
          xSAMOutputPath( nullptr, "Folder (path) for SAM files", "All SAM-output will be written to this folder",
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
        this->xSAMOutputPath.mirror( &rxOtherSet.xSAMOutputPath );
        this->bSAMOutputInReadsFolder.mirror( &rxOtherSet.bSAMOutputInReadsFolder );
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

    /* Delivers pointer to selected parameter-set */
    const Presetting* getSelected( void ) const
    {
        return pSelectedParamSet;
    } // method
}; // class


// parameter set manager is exported from export.cpp....
