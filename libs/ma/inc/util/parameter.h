#pragma once

#include "util/default_parameters.h"
#include "util/exception.h"
#include "util/exported.h"
#include "util/support.h"

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
    const char cShort; // Shorthand character for parameter in command line interface.
    const std::string sDescription; // Description of parameter
    const std::pair<size_t, std::string> sCategory; // Description of parameter

    /* Delete copy constructor */
    AlignerParameterBase( const AlignerParameterBase& rxOtherSet ) = delete;
    AlignerParameterBase( AlignerParameterBase&& other ) = default; // required for getting emplace with maps working
    AlignerParameterBase&
    operator=( AlignerParameterBase&& other ) = default; // required for getting emplace with maps working

    // Callback that is called on demand.
    // By default a parameter is always active.
    std::function<bool( void )> fEnabled = []( ) { return true; };

    AlignerParameterBase( const std::string& sName, const char cShort, const std::string& sDescription,
                          const std::pair<size_t, std::string>& sCategory )
        : sName( sName ), cShort( cShort ), sDescription( sDescription ), sCategory( sCategory )
    {} // constructor

    virtual void mirror( const std::shared_ptr<AlignerParameterBase> pOther )
    {
        throw AnnotatedException( "Mirroring of AlignerParameterBase objects is prohibited." );
    } // method

    virtual std::string type_name( ) const
    {
        return "unknown type";
    } // method

    virtual std::string asText( ) const
    {
        throw AnnotatedException( "AlignerParameterBase has go value" );
    } // method

    virtual void setByText( const std::string& sValueAsString )
    {
        throw AnnotatedException( "Cannot set value of AlignerParameterBase" );
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
    AlignerParameter( const std::string& sName, const char cShort, const std::string& sDescription,
                      const std::pair<size_t, std::string>& sCategory, const VALUE_TYPE value, // initial value
                      std::function<void( const VALUE_TYPE& )> fPredicate = predicateAlwaysOK )
        : AlignerParameterBase( sName, cShort, sDescription, sCategory ), value( value ), fPredicate( fPredicate )
    {} // constructor

    /* Constructor */
    AlignerParameter( const std::string& sName, const std::string& sDescription,
                      const std::pair<size_t, std::string>& sCategory,
                      const VALUE_TYPE value, // initial value
                      std::function<void( const VALUE_TYPE& )> fPredicate = predicateAlwaysOK )
        : AlignerParameter( sName, NO_SHORT_DEFINED, sDescription, sCategory, value, fPredicate )
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
    virtual void setByText( const std::string& sValueAsString )
    {
        VALUE_TYPE newValue = genericStringToValue<VALUE_TYPE>( sValueAsString );
        this->set( newValue );
    } // method

    VALUE_TYPE get( void )
    {
        return this->value;
    } // method

    /** required for exporting get to python do not use in cpp code */
    VALUE_TYPE get_py( void )
    {
        return this->get( );
    } // method

    const VALUE_TYPE get( void ) const
    {
        return this->value;
    } // method

    virtual void mirror( const std::shared_ptr<AlignerParameterBase> pOther )
    {
        auto pOtherDerived = std::dynamic_pointer_cast<const AlignerParameter<VALUE_TYPE>>( pOther );
        // std::cout << "Mirror: " << sName << " from " << this->value << " to " << pOtherDerived->value << std::endl;
        this->value = pOtherDerived->value;
    } // method

    virtual std::string type_name( ) const
    {
        return demangle( typeid( value ).name( ) );
    } // method

    virtual EXPORTED std::string asText( ) const;
}; // class

/* Parameter class for choices (list of different textual alternatives)
 */
template <> class AlignerParameter<AlignerParameterBase::ChoicesType> : public AlignerParameterBase
{
  public:
    AlignerParameterBase::ChoicesType vChoices; // Possible text values of the parameter
    unsigned int uiSelection; // integral value of the selected choice (0 for first choice, 1 for second etc.)

    AlignerParameter( const std::string& sName, const char cShort, const std::string& sDescription,
                      const std::pair<size_t, std::string>& sCategory,
                      const AlignerParameterBase::ChoicesType& rvChoices, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameterBase( sName, cShort, sDescription, sCategory ),
          vChoices( rvChoices ),
          uiSelection( uiSelection )
    {} // constructor

    AlignerParameter( const std::string& sName, const std::string& sDescription,
                      const std::pair<size_t, std::string>& sCategory,
                      const AlignerParameterBase::ChoicesType& rvChoices, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameter( sName, NO_SHORT_DEFINED, sDescription, sCategory, rvChoices, uiSelection )
    {} // constructor

    /* Throws an exception if something goes wrong */
    void set( unsigned int uiNewSelection )
    {
        this->uiSelection = uiNewSelection;
    } // method

    /* Get the corresponding internal string for selection */
    std::string get( void )
    {
        return vChoices.at( uiSelection ).first;
    } // method

    /* Get the corresponding internal string for selection */
    const std::string get( void ) const
    {
        return vChoices.at( uiSelection ).first;
    } // method

    // void mirror( const AlignerParameter<AlignerParameterBase::ChoicesType>* pOther )
    // {
    //     this->uiSelection = pOther->uiSelection;
    // } // method
    //
    virtual void mirror( const std::shared_ptr<AlignerParameterBase> pOther )
    {
        auto pOtherDerived =
            std::dynamic_pointer_cast<const AlignerParameter<AlignerParameterBase::ChoicesType>>( pOther );
        this->uiSelection = pOtherDerived->uiSelection;
    } // method

    virtual std::string type_name( ) const
    {
        std::string s;
        for( auto xPair : vChoices )
            s += xPair.first + "/";
        s.pop_back( );
        return s;
    } // method

    virtual void setByText( const std::string& sValueAsString )
    {
        this->set( stoi( sValueAsString ) );
    } // method

    virtual std::string asText( ) const
    {
        return this->get( );
    } // method
}; // class


/* Parameter class for file system paths */
template <> class AlignerParameter<fs::path> : public AlignerParameterBase
{
  public:
    fs::path xPath; // Current path

    AlignerParameter( const std::string& sName, const char cShort, const std::string& sDescription,
                      const std::pair<size_t, std::string>& sCategory,
                      const fs::path& rxPath, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameterBase( sName, cShort, sDescription, sCategory ), xPath( rxPath )
    {} // constructor

    /* Constructor */
    AlignerParameter( const std::string& sName, const std::string& sDescription,
                      const std::pair<size_t, std::string>& sCategory,
                      const fs::path& rxPath, // choices of the parameter
                      unsigned int uiSelection = 0 )
        : AlignerParameter( sName, NO_SHORT_DEFINED, sDescription, sCategory, rxPath, uiSelection )
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

    // /* Used in the case of objects copies (object mirroring) */
    // void mirror( const AlignerParameter<fs::path>* pOther )
    // {
    //     this->xPath = pOther->xPath;
    // } // method

    virtual void mirror( const std::shared_ptr<AlignerParameterBase> pOther )
    {
        auto pOtherDerived = std::dynamic_pointer_cast<const AlignerParameter<fs::path>>( pOther );
        this->xPath = pOtherDerived->xPath;
    } // method

    virtual std::string type_name( ) const
    {
        return "file_name";
    } // method

    virtual void setByText( const std::string& sValueAsString )
    {
        this->set( sValueAsString );
    } // method

    virtual std::string asText( ) const
    {
        return this->get( ).string( );
    } // method
}; // class

class ParameterSetBase; // forward declaration required for AlignerParameterPointer

template <typename VALUE_TYPE> class AlignerParameterPointer
{
  private:
    void do_register( ParameterSetBase* pPresetting );

  public:
    std::shared_ptr<AlignerParameter<VALUE_TYPE>> pContent;

    template <typename... ARGUMENT_TYPES>
    AlignerParameterPointer( ParameterSetBase* pPresetting, ARGUMENT_TYPES&&... arguments )
        : pContent( std::make_shared<AlignerParameter<VALUE_TYPE>>( arguments... ) )
    {
        do_register( pPresetting );
    } // constructor

    std::shared_ptr<AlignerParameter<VALUE_TYPE>> operator->( )
    {
        return pContent;
    } // operator

    const std::shared_ptr<AlignerParameter<VALUE_TYPE>> operator->( ) const
    {
        return pContent;
    } // operator
}; // class

/* Class for management of aligner parameter sets.
 * Idea: By using Reflections this code could be improved. (e.g. https://www.rttr.org/)
 */
class ParameterSetBase
{
  public:
    // Reflection of all parameter
    std::map<std::string, std::shared_ptr<AlignerParameterBase>> xpAllParameters;
    std::map<char, std::shared_ptr<AlignerParameterBase>> xpParametersByShort;
    std::map<std::pair<size_t, std::string>, std::vector<std::shared_ptr<AlignerParameterBase>>> xpParametersByCategory;

    /**
     * Every parameter has to call this function from it's constructor.
     * We use this in order to generate a map of all available parameters
     */
    void registerParameter( const std::shared_ptr<AlignerParameterBase> pParameter )
    {
        xpAllParameters.emplace( pParameter->sName, pParameter );
        xpParametersByCategory[ pParameter->sCategory ].push_back( pParameter );
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

    bool hasName( const std::string& rParameterName ) const
    {
        return xpAllParameters.count( rParameterName ) > 0;
    } // method

    bool hasShort( const char cX ) const
    {
        return xpParametersByShort.count( cX ) > 0;
    } // method

    std::shared_ptr<AlignerParameterBase> byName( const std::string& rParameterName )
    {
        try
        {
            return xpAllParameters.at( rParameterName );
        } // try
        catch( std::out_of_range& )
        {
            throw AnnotatedException( std::string( "Could not find parameter: " ).append( rParameterName ) );
        } // catch
    } // method

    std::shared_ptr<AlignerParameterBase> byShort( const char cX )
    {
        try
        {
            return xpParametersByShort.at( cX );
        } // try
        catch( std::out_of_range& )
        {
            throw AnnotatedException( "Could not find parameter: " + cX );
        } // catch
    } // method

    /* Mirror the setting of one parameter-set into another */
    void mirror( const ParameterSetBase& rxOtherSet )
    {
        for( auto& rxTup : xpAllParameters )
            rxTup.second->mirror( rxOtherSet.xpAllParameters.at( rxTup.first ) );
    } // method
}; // class

class Presetting : public ParameterSetBase
{
  public:
    AlignerParameterPointer<bool> xUsePairedReads; // If true, work with paired reads
    AlignerParameterPointer<int> xMatch; // score for a DP match (used in SoC width computation)
    AlignerParameterPointer<int> xMisMatch; // score for a DP match (used in SoC width computation)
    AlignerParameterPointer<int> xGap; // @todo
    AlignerParameterPointer<int> xExtend; // @todo
    AlignerParameterPointer<int> xGap2; // @todo
    AlignerParameterPointer<int> xExtend2; // @todo
    AlignerParameterPointer<double> xPairedBonus; // @todo
    AlignerParameterPointer<double> xMeanPairedReadDistance; // @todo
    AlignerParameterPointer<double> xStdPairedReadDistance; // @todo
    AlignerParameterPointer<int> xReportN; // @todo
    AlignerParameterPointer<int> xMaximalSeedAmbiguity; // @todo
    AlignerParameterPointer<int> xMinSeedLength; // @todo
    AlignerParameterPointer<int> xMinAlignmentScore; // @todo
    AlignerParameterPointer<int> xMinimalSeedAmbiguity; // @todo
    AlignerParameterPointer<int> xMinimalSeedSizeDrop; // @todo
    AlignerParameterPointer<int> xMaxNumSoC; // @todo
    AlignerParameterPointer<int> xMinNumSoC; // @todo
    AlignerParameterPointer<int> xSwitchQlen; // @todo
    AlignerParameterPointer<int> xMaxGapArea; // @todo
    AlignerParameterPointer<int> xPadding; // @todo
    AlignerParameterPointer<int> xSoCWidth; // @todo
    AlignerParameterPointer<bool> xSkipAmbiguousSeeds; // @todo
    AlignerParameterPointer<double> xRelMinSeedSizeAmount; // @todo
    AlignerParameterPointer<double> xScoreDiffTolerance; // @todo
    AlignerParameterPointer<int> xMaxScoreLookahead; // @todo
    // AlignerParameterPointer<double> xMinQueryCoverage; // @todo
    AlignerParameterPointer<double> xSoCScoreDecreaseTolerance; // @todo
    AlignerParameterPointer<int> xHarmScoreMin; // @todo
    AlignerParameterPointer<double> xHarmScoreMinRel; // @todo
    AlignerParameterPointer<int> xGenomeSizeDisable; // @todo
    AlignerParameterPointer<bool> xDisableHeuristics; // @todo
    AlignerParameterPointer<bool> xNoSecondary; // @todo
    AlignerParameterPointer<bool> xNoSupplementary; // @todo
    AlignerParameterPointer<double> xMaxDeltaDist; // @todo
    AlignerParameterPointer<int> xMinDeltaDist; // @todo
    AlignerParameterPointer<double> xMaxOverlapSupplementary; // @todo
    AlignerParameterPointer<int> xMaxSupplementaryPerPrim; // @todo
    AlignerParameterPointer<bool> xDisableGapCostEstimationCutting; // @todo
    AlignerParameterPointer<bool> xOptimisticGapCostEstimation; // @todo
    AlignerParameterPointer<int> xSVPenalty; // @todo
    AlignerParameterPointer<int> xMinBandwidthGapFilling; // @todo
    AlignerParameterPointer<int> xBandwidthDPExtension; // @todo
    // AlignerParameterPointer<double> xMaxSVRatio; // @todo
    // AlignerParameterPointer<int> xMinSVDistance; // @todo
    AlignerParameterPointer<int> xZDrop; // @todo

    AlignerParameterPointer<AlignerParameterBase::ChoicesType> xSeedingTechnique; // Seeding Technique

    static const EXPORTED std::pair<size_t, std::string> DP_PARAMETERS;
    static const EXPORTED std::pair<size_t, std::string> HEURISTIC_PARAMETERS;
    static const EXPORTED std::pair<size_t, std::string> SEEDING_PARAMETERS;
    static const EXPORTED std::pair<size_t, std::string> SOC_PARAMETERS;
    static const EXPORTED std::pair<size_t, std::string> PAIRED_PARAMETERS;
    static const EXPORTED std::pair<size_t, std::string> SAM_PARAMETERS;

    /* Delete copy constructor */
    Presetting( const Presetting& rxOtherSet ) = delete;
    Presetting( Presetting&& other ) = default; // required for getting emplace with maps working
    Presetting& operator=( Presetting&& other ) = default; // required for getting emplace with maps working

    /* Constructor */
    Presetting( )
        : // sName( sName ), //
          xUsePairedReads( this, "Use Paired Reads", 'p', "If your reads occur as paired reads, activate this flag.",
                           PAIRED_PARAMETERS, false ),
          xMatch( this, "Match Score",
                  "Match score used in the context of Dynamic Programming and for SoC width computation.",
                  DP_PARAMETERS, 2, checkPositiveValue ),
          xMisMatch( this, "Mismatch Penalty", "Penalty for a Dynamic Programming mismatch", DP_PARAMETERS, 4,
                     checkPositiveValue ),
          xGap( this, "Gap penalty", "First penalty for a DP gap opening.", DP_PARAMETERS, 4 ),
          xExtend( this, "Extend penalty", "First penalty for a DP gap extension.", DP_PARAMETERS, 2 ),
          xGap2( this, "Second gap penalty", "Second penalty for a DP gap opening.", DP_PARAMETERS, 24 ),
          xExtend2( this, "Second extend penalty", "Second penalty for a DP gap extension.", DP_PARAMETERS, 1 ),
          xPairedBonus( this, "Score factor for paired reads",
                        "This factor is multiplied to the score of successfully paired reads. Used in the context of "
                        "the computation of the mapping quality and for picking optimal alignment pairs. [val] < 1 "
                        "results in penalty; [val] > 1 results in bonus.",
                        PAIRED_PARAMETERS, 1.25 ),
          xMeanPairedReadDistance(
              this, "Mean distance of paired reads", 'd',
              "Two reads can be paired if they are within mean +- (standard deviation)*3 distance from one another on "
              "the expected strands (depends on Use Mate Pair on/off) Used in the context of the computation of the "
              "mapping quality and for picking optimal alignment pairs.",
              PAIRED_PARAMETERS, 400 ),
          xStdPairedReadDistance(
              this, "Standard deviation of paired reads", 'S',
              "<val> represents the standard deviation for the distance between paired reads. Used in the context of "
              "the computation of the mapping quality and for picking optimal alignment pairs.",
              PAIRED_PARAMETERS, 150 ),
          xReportN( this, "Max. number of Reported alignments", 'n',
                    "Do not output more than <val> alignments. 0 = no limit.", SAM_PARAMETERS, 0 ),
          xMaximalSeedAmbiguity( this, "Maximal ambiguity", "Maximal ambiguity of seeds.", SEEDING_PARAMETERS, 500 ),
          xMinSeedLength( this, "Minimal Seed length", 'l', "Minimal seed length.", SEEDING_PARAMETERS, 16 ),
          xMinAlignmentScore( this, "Minimal alignment score",
                              "Suppress the output of alignments with a score below val.", SAM_PARAMETERS, 75 ),
          xMinimalSeedAmbiguity( this, "Min ambiguity", "Stop the extension process if seeds are less ambiguous.",
                                 SEEDING_PARAMETERS, 0 ),
          xMinimalSeedSizeDrop( this, "Seeding drop-off A - Drop-off min seed size",
                                "Heuristic runtime optimization: For a given read R, let N be the number of seeds of "
                                "size >= [val]. Discard R, if N < [length(R)] * [Seeding drop-off B].",
                                SEEDING_PARAMETERS, 15 ),
          xMaxNumSoC( this, "Maximal Number of SoC's", 'N', "Only consider the <val> best scored SoC's. 0 = no limit.",
                      SOC_PARAMETERS, 30 ),
          xMinNumSoC( this, "Min Number SoC’s",
                      'M', "Always consider the first <val> SoC's no matter the Heuristic optimizations.",
                      SOC_PARAMETERS, 1 ),
          xSwitchQlen( this, "Harmonization Score dropoff - Minimal Query length",
                       "For reads of length >= [val]: Ignore all SoC's with harmonization scores lower than the "
                       "current maximal score. 0 = disabled.",
                       HEURISTIC_PARAMETERS, 800 ),
          xMaxGapArea( this, "Maximal Gap Area", "Split alignments in harmonization if gap area is larger than <val>.",
                       HEURISTIC_PARAMETERS, 10000 ),
          xPadding( this, "Padding",
                    "Padding for DP extensions. Maximal area in front and back of the alignment that gets checked "
                    "using DP extensions.",
                    DP_PARAMETERS, 1000 ),
          xSoCWidth( this, "Fixed SoC Width",
                     "Set the SoC width to a fixed value. 0 = use the formula given in the paper. This parameter is "
                     "intended for debugging purposes.",
                     SOC_PARAMETERS, 0 ),
          xSkipAmbiguousSeeds( this, "Skip ambiguous seeds",
                               "Enabled: Discard all seeds that are more ambiguous than [max ambiguity]. Disabled: "
                               "sample [max ambiguity] random seeds from too ambiguous seeds.",
                               SEEDING_PARAMETERS, false ),
          xRelMinSeedSizeAmount( this, "Seeding drop off B - Drop-off percentage",
                                 "Heuristic runtime optimization: Percentage for seed drop-off calculation. For more "
                                 "information see parameter [Seeding drop-off A]. ",
                                 SEEDING_PARAMETERS, 0.005 ),
          xScoreDiffTolerance(
              this, "Harmonization Drop-off A - Score difference",
              "Let x be the maximal encountered harmonization score. Stop harmonizing further SoC's if there are "
              "<Harmonization Drop-off B> SoC's with lower scores than x-<readlength>*<val> in a row.",
              HEURISTIC_PARAMETERS, 0.0001 ),
          xMaxScoreLookahead( this, "Harmonization Drop-off B - Lookahead", "See Harmonization Drop-off A.",
                              HEURISTIC_PARAMETERS, 3 ),
          // xMinQueryCoverage( this, "", "", 1.1 ),
          xSoCScoreDecreaseTolerance( this, "SoC Score Drop-off",
                                      "Let x be the maximal encountered SoC score. Stop harmonizing SoC's if there is "
                                      "a SoC with a score lower than <val>*x.",
                                      SEEDING_PARAMETERS, 0.1 ),
          xHarmScoreMin( this, "Minimal Harmonization Score",
                         "Discard all harmonized SoC's with scores lower than <val>.", SEEDING_PARAMETERS, 18 ),
          xHarmScoreMinRel( this, "Relative Minimal Harmonization Score",
                            "Discard all harmonized SoC's with scores lower than length(read)*<val>.",
                            SEEDING_PARAMETERS, 0.002 ),
          xGenomeSizeDisable( this, "Minimum Genome Size for Heuristics",
                              "Some heuristics can only be applied on long enough genomes. Disables: SoC score "
                              "Drop-off if the genome is shorter than <val>.",
                              SEEDING_PARAMETERS, 10000000 ),
          xDisableHeuristics( this, "Disable All Heuristics",
                              "Disables all runtime heuristics. (Intended for debugging.)", SEEDING_PARAMETERS, false ),
          xNoSecondary( this, "Omit secondary alignments", "Suppress the output of secondary alignments.",
                        SAM_PARAMETERS, false ),
          xNoSupplementary( this, "Omit supplementary alignments", "Suppress the output of supplementary alignments.",
                            SAM_PARAMETERS, false ),
          xMaxDeltaDist( this, "Artifact Filter A - Maximal Delta Distance",
                         "Filter seeds if the difference between the delta distance to it's predecessor and successor "
                         "is less then [val] percent (set to 1 to disable filter) and the delta distance to it's pre- "
                         "and successor is more than [Artifact Filter B] nt.",
                         HEURISTIC_PARAMETERS, 0.1 ),
          xMinDeltaDist( this, "Artifact Filter B - Minimal Delta Distance", "See Artifact Filter A",
                         HEURISTIC_PARAMETERS, 16 ),
          xMaxOverlapSupplementary(
              this, "Maximal supplementary overlap",
              "An non-primary alignment A is considered supplementary, if less than val percent of A overlap with the "
              "primary alignment on the query. Otherwise A is considered secondary.",
              SAM_PARAMETERS, 0.1 ),
          xMaxSupplementaryPerPrim( this, "Number Supplementary alignments",
                                    "Maximal Number of supplementary alignments per primary alignment.", SAM_PARAMETERS,
                                    1 ),
          xDisableGapCostEstimationCutting( this, "Pick Local Seed Set A - Enabled",
                                            "<val> = true enables local seed set computiaion.", HEURISTIC_PARAMETERS,
                                            false ),
          xOptimisticGapCostEstimation(
              this, "Pick Local Seed Set B - Optimistic Gap Estimation",
              "After the harmonization MA checks weather it is possible to compute a positively scored alignment from "
              "the seed set. Gaps between seeds can be estimated in two ways: Optimistic [true]: Assume that the gap "
              "can be filled using merely matches and a single insertion/deletion. Pessimistic [false]: Assume that "
              "the gap can be filled using matches and mismatches that add up to a score of 0 and a single "
              "insertion/deletion.",
              HEURISTIC_PARAMETERS, true ),
          xSVPenalty( this, "Pick Local Seed Set C - Maximal Gap Penalty",
                      "Maximal Gap cost penalty during local seed set computiaion.", HEURISTIC_PARAMETERS, 100 ),
          xMinBandwidthGapFilling( this, "Minimal bandwidth in gaps", "Minimal bandwidth for DP in gaps between seeds.",
                                   DP_PARAMETERS, 20 ),
          xBandwidthDPExtension( this, "Bandwidth for extensions", "Bandwidth for DP extensions.", DP_PARAMETERS, 512 ),
          // xMaxSVRatio( this, "", "", "CATEGORY",0.01 ),
          // xMinSVDistance( this, "", "", "CATEGORY",500 ),
          xZDrop( this, "Z Drop", "If the DP score drops faster than <val> stop the extension process.", DP_PARAMETERS,
                  200 ),

          xSeedingTechnique( this, "Seeding Technique", 's', "Technique used for the initial seeding.",
                             SEEDING_PARAMETERS,
                             AlignerParameterBase::ChoicesType{{"maxSpan", "Maximally Spanning"}, {"SMEMs", "SMEMs"}} )
    {

        xMeanPairedReadDistance->fEnabled = [this]( void ) { return this->xUsePairedReads->get( ) == true; };
        xStdPairedReadDistance->fEnabled = [this]( void ) { return this->xUsePairedReads->get( ) == true; };
        xPairedBonus->fEnabled = [this]( void ) { return this->xUsePairedReads->get( ) == true; };
    } // constructor

    /* Named copy Constructor */
    Presetting( const Presetting& rxOtherSet, const std::string& sName ) : Presetting( )
    {
        this->mirror( rxOtherSet );
    } // copy constructor

    /* True if the parameter-set uses paired reads */
    bool usesPairedReads( void )
    {
        return xUsePairedReads->value;
    } // method
}; // class

template <typename VALUE_TYPE> // @todo maybe this can be part of the constructor ?
void AlignerParameterPointer<VALUE_TYPE>::do_register( ParameterSetBase* pPresetting )
{
    if( pPresetting != nullptr )
        pPresetting->registerParameter( this->pContent );
} // method

/* Parameter which occur only once.
 * (Not sequencing technique / configuration related parameters)
 */
class GeneralParameter : public ParameterSetBase
{
  public:
    // FIXME: Set this to a standard path in the beginning (~home/ma/SAM)
    AlignerParameterPointer<bool> bSAMOutputInReadsFolder; // SAM Output in the same folder as the reads.
    AlignerParameterPointer<fs::path> xSAMOutputPath; // folder path
    AlignerParameterPointer<bool> pbUseMaxHardareConcurrency; // Exploit all cores
    AlignerParameterPointer<int> piNumberOfThreads; // selected number of threads
    AlignerParameterPointer<bool> pbPrintHelpMessage; // Print the help message to stdout

    /* Constructor */
    GeneralParameter( )
        : bSAMOutputInReadsFolder( this, "SAM files in same folder as reads",
                                   "If set, the SAM files are written in the folder of the reads.",
                                   std::pair<size_t, std::string>( ), true ),
          xSAMOutputPath( this, "Folder (path) for SAM files", 'o', "All SAM-output will be written to this folder",
                          std::pair<size_t, std::string>( ), fs::temp_directory_path( ) ),
          pbUseMaxHardareConcurrency( this, "Use all processor cores",
                                      "The number of threads used for aligning is chosen to be identical to the number "
                                      "of your processor cores.",
                                      std::pair<size_t, std::string>( ), true ),
          piNumberOfThreads( this, "Number of threads", "Number of threads used in the context of alignments.",
                             std::pair<size_t, std::string>( ), 1 ),
          pbPrintHelpMessage( this, "Print help to stdout", 'h', "Number of threads used in the context of alignments.",
                              std::pair<size_t, std::string>( ), false )

    {
        xSAMOutputPath->fEnabled = [this]( void ) { return this->bSAMOutputInReadsFolder->get( ) == false; };
        piNumberOfThreads->fEnabled = [this]( void ) { return this->pbUseMaxHardareConcurrency->get( ) == false; };
    } // constructor

    /* Named copy Constructor */
    GeneralParameter( const GeneralParameter& rxOtherSet ) : GeneralParameter( )
    {
        mirror( rxOtherSet );
    } // copy constructor
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
        xParametersSets[ "Illumina Paired" ].xUsePairedReads->set( true );

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

    /** required for exporting get to python do not use in cpp code */
    Presetting* getSelected_py( void )
    {
        return this->getSelected( );
    } // method

    /* Delivers pointer to selected parameter-set */
    const Presetting* getSelected( void ) const
    {
        return pSelectedParamSet;
    } // method

    std::shared_ptr<AlignerParameterBase> byName( const std::string& rParameterName )
    {
        if( xGlobalParameterSet.hasName( rParameterName ) )
            return xGlobalParameterSet.byName( rParameterName );
        if( this->getSelected( )->hasName( rParameterName ) )
            return this->getSelected( )->byName( rParameterName );
        throw AnnotatedException( std::string( "Could not find parameter: " ).append( rParameterName ) );
    } // method

    std::shared_ptr<AlignerParameterBase> byShort( const char cX )
    {
        if( xGlobalParameterSet.hasShort( cX ) )
            return xGlobalParameterSet.byShort( cX );
        if( this->getSelected( )->hasShort( cX ) )
            return this->getSelected( )->byShort( cX );
        throw AnnotatedException( "Could not find parameter: " + cX );
    } // method
}; // class


// parameter set manager is exported from export.cpp....
