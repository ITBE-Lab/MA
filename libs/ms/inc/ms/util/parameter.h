#pragma once

#include "util/exported.h"
#include "util/support.h"
#include <stdexcept>

#if( defined( __GNUC__ ) && ( __GNUC__ < 8 ) ) && !( defined( __clang__ ) && ( __clang_major__ >= 10 ) )
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

#include <functional>
#include <locale>
#include <map>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

/* Generic conversion of string to different value types.
 * (Specialized for some types)
 */
template <typename VALUE_TYPE> DLL_PORT( MS ) VALUE_TYPE genericStringToValue( const std::string& sString );

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
static void checkPositiveValue( const int& iValue, const std::string& rsParameterName )
{
    if( iValue < 0 )
        throw std::out_of_range( "Only positive values are allowed for the parameter " + rsParameterName + "." );
} // static method

static void checkPositiveDoubleValue( const double& dValue, const std::string& rsParameterName )
{
    if( dValue < 0 )
        throw std::out_of_range( "Only positive values are allowed for the parameter " + rsParameterName + "." );
} // static method


/* Base class of all aligner parameter classes */
class AlignerParameterBase
{
  public:
    typedef std::vector<std::pair<std::string, std::string>> ChoicesType;
    static const char NO_SHORT_DEFINED = ' ';

    const std::string sName; // Name of parameter
    const char cShort; // Shorthand character for parameter in command line interface.
    std::string sDescription; // Description of parameter
    const std::pair<size_t, std::string> sCategory; // Description of parameter
    const std::string sSetDesc; // Description of the set this parameter belongs to

    /* Delete copy constructor */
    AlignerParameterBase( const AlignerParameterBase& rxOtherSet ) = delete;
    AlignerParameterBase( AlignerParameterBase&& other ) = default; // required for getting emplace with maps working
    //AlignerParameterBase& operator=( AlignerParameterBase&& other ) = default; // required for getting emplace with maps working

    // Callback that is called on demand.
    // By default a parameter is always active.
    std::function<bool( void )> fEnabled = []( ) { return true; };

    AlignerParameterBase( const std::string& sSetDesc, const std::string& sName, const char cShort,
                          const std::string& sDescription, const std::pair<size_t, std::string>& sCategory )
        : sName( sName ), cShort( cShort ), sDescription( sDescription ), sCategory( sCategory ), sSetDesc( sSetDesc )
    {} // constructor

    virtual void mirror( const std::shared_ptr<AlignerParameterBase> pOther )
    {
        throw std::runtime_error( "Mirroring of AlignerParameterBase objects is prohibited." );
    } // method

    virtual std::string type_name( ) const
    {
        return "unknown type";
    } // method

    virtual std::string asText( ) const
    {
        throw std::runtime_error( "AlignerParameterBase has go value" );
    } // method

    virtual void setByText( const std::string& sValueAsString )
    {
        throw std::runtime_error( "Cannot set value of AlignerParameterBase" );
    } // method

    /* Clang requires a virtual destructor for classes comprising virtual methods */
    virtual ~AlignerParameterBase( )
    {} // destructor
}; // class


/* Generic aligner parameter class.
 * Suitable for primitive types like int, double etc.
 */
template <typename VALUE_TYPE> class AlignerParameter : public AlignerParameterBase
{
  public:
    VALUE_TYPE value; // Current Value of parameter
    // Universal predicate
    static void predicateAlwaysOK( const VALUE_TYPE&, const std::string& )
    {} // static method
    std::function<void( const VALUE_TYPE&, const std::string& )> fPredicate;

    /* Constructor */
    AlignerParameter( const std::string& sSetDesc, const std::string& sName, const char cShort,
                      const std::string& sDescription, const std::pair<size_t, std::string>& sCategory,
                      const VALUE_TYPE value, // initial value
                      std::function<void( const VALUE_TYPE&, const std::string& )> fPredicate = predicateAlwaysOK )
        : AlignerParameterBase( sSetDesc, sName, cShort, sDescription, sCategory ),
          value( value ),
          fPredicate( fPredicate )
    {} // constructor

    /* Constructor */
    AlignerParameter( const std::string& sSetDesc, const std::string& sName, const std::string& sDescription,
                      const std::pair<size_t, std::string>& sCategory,
                      const VALUE_TYPE value, // initial value
                      std::function<void( const VALUE_TYPE&, const std::string& )> fPredicate = predicateAlwaysOK )
        : AlignerParameter( sSetDesc, sName, NO_SHORT_DEFINED, sDescription, sCategory, value, fPredicate )
    {} // constructor

    /* Constructor */
    AlignerParameter( const std::string& sSetDesc, const std::string& sName, const std::string& sDescription,
                      const size_t& uiCategoryIndex, const std::string& sCategoryName, const VALUE_TYPE value )
        : AlignerParameter( sSetDesc, sName, NO_SHORT_DEFINED, sDescription,
                            std::make_pair( uiCategoryIndex, sCategoryName ), value )
    {} // constructor

    /* Throws an exception if something goes wrong */
    void set( VALUE_TYPE newValue )
    {
        /* Check the value. ( Error is indicated by exception )*/
        fPredicate( newValue, sName );

        /* In the case of an exception this line is skipped !*/
        this->value = newValue;
    } // method

    virtual std::string type_name( ) const
    {
        return demangle( typeid( value ).name( ) );
    } // method

    /* Throws an exception if something goes wrong */
    virtual void setByText( const std::string& sValueAsString )
    {
        try
        {
            VALUE_TYPE newValue = genericStringToValue<VALUE_TYPE>( sValueAsString );
            this->set( newValue );
        } // try
        catch( ... )
        {
            throw std::runtime_error( std::string( "Only " )
                                          .append( this->type_name( ) )
                                          .append( std::string( " values are allowed for the parameter " ) )
                                          .append( this->sName ) );
        } // catch
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

    virtual DLL_PORT( MS ) std::string asText( ) const;
}; // class

/* Parameter class for choices (list of different textual alternatives)
 */
template <> class AlignerParameter<AlignerParameterBase::ChoicesType> : public AlignerParameterBase
{
  public:
    AlignerParameterBase::ChoicesType vChoices; // Possible text values of the parameter
    unsigned int uiSelection; // integral value of the selected choice (0 for first choice, 1 for second etc.)

    AlignerParameter( const std::string& sSetDesc, const std::string& sName, const char cShort,
                      const std::string& sDescription, const std::pair<size_t, std::string>& sCategory,
                      const AlignerParameterBase::ChoicesType& rvChoices, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameterBase( sSetDesc, sName, cShort, sDescription, sCategory ),
          vChoices( rvChoices ),
          uiSelection( uiSelection )
    {} // constructor

    AlignerParameter( const std::string& sSetDesc, const std::string& sName, const std::string& sDescription,
                      const std::pair<size_t, std::string>& sCategory,
                      const AlignerParameterBase::ChoicesType& rvChoices, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameter( sSetDesc, sName, NO_SHORT_DEFINED, sDescription, sCategory, rvChoices, uiSelection )
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

    /* unoverloaded function for python */
    std::string get_py( void )
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
        // std::string s;
        // for( auto xPair : vChoices )
        //    s += xPair.first + "/";
        // s.pop_back( );
        // return s;
        return "name";
    } // method

    virtual void setByText( const std::string& sValueAsString )
    {
        this->set( stoi( sValueAsString ) );
    } // method

    virtual std::string asText( ) const
    {
        return this->get( );
    } // method

    /* Clang requires a virtual destructor for classes comprising virtual methods */
    virtual ~AlignerParameter( )
    {} // destructor
}; // class


/* Parameter class for file system paths */
template <> class AlignerParameter<fs::path> : public AlignerParameterBase
{
  public:
    fs::path xPath; // Current path

    AlignerParameter( const std::string& sSetDesc, const std::string& sName, const char cShort,
                      const std::string& sDescription, const std::pair<size_t, std::string>& sCategory,
                      const fs::path& rxPath, // choices of the parameter
                      unsigned int uiSelection = 0 ) // initially selected choice
        : AlignerParameterBase( sSetDesc, sName, cShort, sDescription, sCategory ), xPath( rxPath )
    {} // constructor

    /* Constructor */
    AlignerParameter( const std::string& sSetDesc, const std::string& sName, const std::string& sDescription,
                      const std::pair<size_t, std::string>& sCategory,
                      const fs::path& rxPath, // choices of the parameter
                      unsigned int uiSelection = 0 )
        : AlignerParameter( sSetDesc, sName, NO_SHORT_DEFINED, sDescription, sCategory, rxPath, uiSelection )
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

    /* Clang requires a virtual destructor for classes comprising virtual methods */
    virtual ~AlignerParameter( )
    {} // destructor
}; // class

class ParameterSetBase; // forward declaration required for AlignerParameterPointer

template <typename VALUE_TYPE> class AlignerParameterPointer
{
  private:
    void do_register( ParameterSetBase* pPresetting );
    std::string get_desc( ParameterSetBase* pPresetting );

  public:
    std::shared_ptr<AlignerParameter<VALUE_TYPE>> pContent;

    template <typename... ARGUMENT_TYPES>
    AlignerParameterPointer( ParameterSetBase* pPresetting, ARGUMENT_TYPES&&... arguments )
        : pContent( std::make_shared<AlignerParameter<VALUE_TYPE>>( get_desc( pPresetting ), arguments... ) )
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

    static inline std::string uniqueParameterName( std::string sSearch )
    {
        std::locale loc;
        for( std::string::size_type i = 0; i < sSearch.length( ); ++i )
            sSearch[ i ] = std::tolower( sSearch[ i ], loc );
        sSearch.erase(
            std::remove_if( sSearch.begin( ), sSearch.end( ), []( const char cX ) { return cX == ' ' || cX == '_'; } ),
            sSearch.end( ) );
        return sSearch;
    } // method

    virtual std::string getSetDescription( )
    {
        throw std::runtime_error( "This function should have been overloaded!" );
    } // method

    /**
     * Every parameter has to call this function from it's constructor.
     * We use this in order to generate a map of all available parameters
     */
    void registerParameter( const std::shared_ptr<AlignerParameterBase> pParameter )
    {
        xpAllParameters.emplace( uniqueParameterName( pParameter->sName ), pParameter );

        xpParametersByCategory[ pParameter->sCategory ].push_back( pParameter );
        if( pParameter->cShort != pParameter->NO_SHORT_DEFINED )
        {
            auto xEmplaceRet = xpParametersByShort.emplace( pParameter->cShort, pParameter );
            if( !xEmplaceRet.second ) // there is another parameter with this shorthand...
                throw std::runtime_error( std::string( "Shorthand: " )
                                              .append( std::to_string( pParameter->cShort ) )
                                              .append( " used twice: 1) " )
                                              .append( xEmplaceRet.first->second->sName )
                                              .append( " 2)" )
                                              .append( pParameter->sName ) );
        }
    } // method

    void unregisterParameter( std::shared_ptr<AlignerParameterBase> pParameter )
    {
        xpAllParameters.erase( uniqueParameterName( pParameter->sName ) );
        auto xIt1 = xpParametersByShort.begin( );
        while( xIt1 != xpParametersByShort.end( ) )
        {
            if( xIt1->second == pParameter )
                xpParametersByShort.erase( xIt1++ );
            else
                xIt1++;
        } // while
        for( auto& rPair : xpParametersByCategory )
            rPair.second.erase( std::remove_if( rPair.second.begin( ), rPair.second.end( ),
                                                [ &pParameter ]( const std::shared_ptr<AlignerParameterBase> pX ) {
                                                    return ( pX == pParameter );
                                                } ),
                                rPair.second.end( ) );
    } // method

    bool hasName( const std::string& rParameterName ) const
    {
        return xpAllParameters.count( uniqueParameterName( rParameterName ) ) > 0;
    } // method

    bool hasShort( const char cX ) const
    {
        return xpParametersByShort.count( cX ) > 0;
    } // method

    std::shared_ptr<AlignerParameterBase> byName( const std::string& rParameterName )
    {
        try
        {
            return xpAllParameters.at( uniqueParameterName( rParameterName ) );
        } // try
        catch( std::out_of_range& )
        {
            throw std::runtime_error( std::string( "Could not find parameter: " ).append( rParameterName ) );
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
            throw std::runtime_error( std::string( "Could not find parameter: " ) += cX );
        } // catch
    } // method

    /* Mirror the setting of one parameter-set into another */
    void mirror( const ParameterSetBase& rxOtherSet )
    {
        for( auto& rxTup : xpAllParameters )
            rxTup.second->mirror( rxOtherSet.xpAllParameters.at( rxTup.first ) );
    } // method

    /* Clang requires a virtual destructor for classes comprising virtual methods */
    virtual ~ParameterSetBase( )
    {} // destructor
}; // class


const std::pair<size_t, std::string> MINIMIZER_PARAMETERS = std::make_pair( 7, "Minimizer" );
const std::pair<size_t, std::string> DP_PARAMETERS = std::make_pair( 5, "Dynamic Programming" );
const std::pair<size_t, std::string> HEURISTIC_PARAMETERS = std::make_pair( 4, "Heuristics" );
const std::pair<size_t, std::string> SEEDING_PARAMETERS = std::make_pair( 1, "Seeding" );
const std::pair<size_t, std::string> SOC_PARAMETERS = std::make_pair( 2, "Strip of Consideration" );
const std::pair<size_t, std::string> PAIRED_PARAMETERS = std::make_pair( 0, "Paired Reads" );
const std::pair<size_t, std::string> SAM_PARAMETERS = std::make_pair( 3, "SAM Output" );
const std::pair<size_t, std::string> GENERAL_PARAMETER = std::make_pair( 0, "General Parameter" );
const std::pair<size_t, std::string> SV_PARAMETERS = std::make_pair( 6, "SV Parameter - ONLY RELEVANT FOR MSV" );


class Presetting : public ParameterSetBase
{
  public:
    const std::string sName;
    // DP options:
    AlignerParameterPointer<int> xPadding; // Padding
    AlignerParameterPointer<int> xBandwidthDPExtension; // Bandwidth for extensions
    AlignerParameterPointer<int> xMinBandwidthGapFilling; // Minimal bandwidth in gaps
    AlignerParameterPointer<int> xZDrop; // Z Drop
    AlignerParameterPointer<bool> xSearchInversions; // Search fro small inversions using DP
    AlignerParameterPointer<int> xZDropInversion; // Z Drop for inversion detection

    // Paired reads options:
    AlignerParameterPointer<bool> xUsePairedReads; // Use Paired Reads
    AlignerParameterPointer<double> xMeanPairedReadDistance; // Mean distance of paired reads
    AlignerParameterPointer<double> xStdPairedReadDistance; // Standard deviation of paired reads
    AlignerParameterPointer<double> xPairedBonus; // Score factor for paired reads

    // Seeding options:
    AlignerParameterPointer<AlignerParameterBase::ChoicesType> xSeedingTechnique; // Seeding Technique
    AlignerParameterPointer<int> xMinSeedLength; // Minimal Seed length
    AlignerParameterPointer<int> xMinimalSeedAmbiguity; // Minimal ambiguity
    AlignerParameterPointer<int> xMaximalSeedAmbiguity; // Maximal ambiguity
    AlignerParameterPointer<bool> xSkipAmbiguousSeeds; // Skip ambiguous seeds
    AlignerParameterPointer<int> xMinimalSeedSizeDrop; // Seeding drop-off A - Drop-off min seed size
    AlignerParameterPointer<double> xRelMinSeedSizeAmount; // Seeding drop off B - Drop-off factor

    // SoC Options:
    AlignerParameterPointer<int> xMaxNumSoC; // Maximal Number of SoCs
    AlignerParameterPointer<int> xMinNumSoC; // Min Number SoCs
    AlignerParameterPointer<int> xSoCWidth; // Fixed SoC Width
    AlignerParameterPointer<bool> xRectangularSoc; // Rectangular SoC

    // SAM Options:
    AlignerParameterPointer<int> xReportN; // Maximal number of Reported alignments
    AlignerParameterPointer<int> xMinAlignmentScore; // Minimal alignment score
    AlignerParameterPointer<bool> xNoSecondary; // Omit secondary alignments
    AlignerParameterPointer<bool> xNoSupplementary; // Omit supplementary alignments
    AlignerParameterPointer<double> xMaxOverlapSupplementary; // Maximal supplementary overlap
    AlignerParameterPointer<int> xMaxSupplementaryPerPrim; // Number Supplementary alignments
    AlignerParameterPointer<bool> xEmulateNgmlrTags; // Emulate NGMLR's tag output
    AlignerParameterPointer<bool> xOutputMCigar; // Output M symbol in cigar
    AlignerParameterPointer<bool> xCGTag; // Output long CIGARS in CG tag
    AlignerParameterPointer<bool> xSoftClip; // use soft clipping in alignments

    // SV caller options
    // @todo depre start
    AlignerParameterPointer<int> xMaxDeltaDistanceInCLuster; // maximal distance between clusters
    AlignerParameterPointer<int> xPaddingReSeeding; // re seeding padding
    // @todo depre end
    AlignerParameterPointer<int> xMinSeedSizeSV; // minimal seed size for the sv caller
    AlignerParameterPointer<int> xMaxAmbiguitySv; // max seed ambiguity for the sv caller
    AlignerParameterPointer<bool> xDoDummyJumps; // compute dummy jumps for first & last seed of read
    AlignerParameterPointer<int> xMinDistDummy; // minimal distance for seeds from read edge for dummy jumps
    AlignerParameterPointer<int> xMaxDistDummy; // minimal distance for seeds from read edge for dummy jumps
    AlignerParameterPointer<bool> xRevCompPairedReadMates; // reverse complement all paried read mates
    AlignerParameterPointer<int> xMaxSizeReseed; // minimal edge size for which we reseed
    AlignerParameterPointer<int> xMinSizeEdge; // minimal edge size for which we keep the edge (runtime improvement)
    AlignerParameterPointer<int> xMaxSizeEdge; // minimal edge size for which we keep the edge (runtime improvement)
    AlignerParameterPointer<int> xMaxFuzzinessFilter; // maximal fuzziness for sv calls
    AlignerParameterPointer<int> xMaxSuppNtShortCallFilter; // maximal number of nt for short call low support filter
    AlignerParameterPointer<int> xMaxCallSizeShortCallFilter; // maximal call size for short call low support filter
    AlignerParameterPointer<int> xMaxRegionSizeForAmbiguityFilter; // max size of the 
    AlignerParameterPointer<int> xMaxRefAmbiguityJump; // Max Ref Ambiguity Jump
    AlignerParameterPointer<int> xMMFilterMaxOcc; // xMMFilterMaxOcc
    AlignerParameterPointer<int> xMinNtInSoc; // Min NT in SoC
    AlignerParameterPointer<double> xMinNtInSocRelative; // Min NT in SoC Relative
    AlignerParameterPointer<int> xMinNtAfterReseeding; // Min NT after reseeding
    AlignerParameterPointer<double> xMinNtAfterReseedingRelative; // Min NT after reseeding Relative
    AlignerParameterPointer<int> xMinReadsInCall; // Min Reads in Call

    // Heuristic Options:
    AlignerParameterPointer<double> xSoCScoreDecreaseTolerance; // SoC Score Drop-off
    AlignerParameterPointer<int> xHarmScoreMin; // Minimal Harmonization Score
    AlignerParameterPointer<double> xHarmScoreMinRel; //
    AlignerParameterPointer<double> xScoreDiffTolerance; //
    AlignerParameterPointer<int> xMaxScoreLookahead; //
    AlignerParameterPointer<int> xSwitchQlen; //
    AlignerParameterPointer<double> xMaxDeltaDist; //
    AlignerParameterPointer<int> xMinDeltaDist; //
    AlignerParameterPointer<bool> xDisableGapCostEstimationCutting; //
    AlignerParameterPointer<bool> xOptimisticGapCostEstimation; //
    AlignerParameterPointer<int> xMaxGapArea; //
    AlignerParameterPointer<int> xGenomeSizeDisable; //
    AlignerParameterPointer<bool> xDisableHeuristics; //

    // Minimizer parameters:
    AlignerParameterPointer<short> xMinimizerK;
    AlignerParameterPointer<short> xMinimizerW;
    AlignerParameterPointer<short> xMinimizerFlag;
    AlignerParameterPointer<short> xMinimizerBucketBits;
    AlignerParameterPointer<int> xMinimizerMiniBatchSize;
    AlignerParameterPointer<uint64_t> xMinimizerBatchSize;

    /* Delete copy constructor */
    Presetting( const Presetting& rxOtherSet ) = delete;
    Presetting( Presetting&& other ) = default; // required for getting emplace with maps working
    //Presetting& operator=( Presetting&& other ) = default; // required for getting emplace with maps working

    /* Constructor */
    Presetting( const std::string sName )
        : //
          sName( sName ),
          xPadding( this, "Padding",
                    "If an alignment does not reach its read's endpoints, the missing parts can be computed via "
                    "dynamic programming. If the length of the missing parts is smaller than 'Padding', dynamic "
                    "programming is used to extend the alignment towards the endpoints of the read. Otherwise, the "
                    "unaligned parts of the read are ignored and the alignment stays unextended.",
                    DP_PARAMETERS, 1000, checkPositiveValue ),
          xBandwidthDPExtension(
              this, "Bandwidth for Extensions",
              "Bandwidth used in the context of the extension of an alignment towards the endpoints of its "
              "read. (See 'Padding')",
              DP_PARAMETERS, 512, checkPositiveValue ),
          xMinBandwidthGapFilling( this, "Minimal Bandwidth in Gaps",
                                   "Gaps between seeds are generally filled using dynamic programming. This option "
                                   "determines the minimal "
                                   "bandwidth used in the context of fillin gaps.",
                                   DP_PARAMETERS, 20, checkPositiveValue ),
          xZDrop( this, "Z Drop",
                  "If the running score during dynamic programming drops faster than <val> stop the extension process.",
                  DP_PARAMETERS, 200, checkPositiveValue ),
          xSearchInversions( this, "Detect Small Inversions",
                             "Use dynamic programming to search for small inversions that do not contain any seeds. "
                             "(Flag disabled = off)",
                             DP_PARAMETERS, false ),
          xZDropInversion(
              this, "Z Drop Inversions",
              "Check for an inversion, if the running score during dynamic programming drops faster than <val>.",
              DP_PARAMETERS, 100, checkPositiveValue ),

          // Paired Reads:
          xUsePairedReads( this, "Use Paired Reads", "For paired reads set this flag to true.", PAIRED_PARAMETERS,
                           false ),
          xMeanPairedReadDistance( this, "Mean Distance of Paired Reads", 'd',
                                   "Two reads can be paired, if they are within mean +- (standard deviation)*3 "
                                   "distance from one another on "
                                   "the expected strands (depends on Use Mate Pair on/off) Used in the context of "
                                   "the computation of the "
                                   "mapping quality and for picking optimal alignment pairs.",
                                   PAIRED_PARAMETERS, 400 ),
          xStdPairedReadDistance( this, "Standard Deviation of Paired Reads", 'S',
                                  "<val> represents the standard deviation for the distance between paired reads. "
                                  "Used in the context of "
                                  "the computation of the mapping quality and for picking optimal alignment pairs.",
                                  PAIRED_PARAMETERS, 150, checkPositiveDoubleValue ),
          xPairedBonus( this, "Score Factor for Paired Reads",
                        "This factor is multiplied to the score of successfully paired reads. Used in the context of "
                        "the computation of the mapping quality and for picking optimal alignment pairs. <val> < 1 "
                        "results in penalty; <val> > 1 results in bonus.",
                        PAIRED_PARAMETERS, 1.25, checkPositiveDoubleValue ),

          // Seeding:
          xSeedingTechnique( this, "Seeding Technique", 's',
                             "Technique used for the initial seeding. Available techniques are: maxSpan and SMEMs.",
                             SEEDING_PARAMETERS,
                             AlignerParameterBase::ChoicesType{
                                 { "maxSpan", "Maximally Spanning" }, { "SMEMs", "SMEMs" }, { "MEMs", "MEMs" } } ),
          xMinSeedLength( this, "Minimal Seed Length", 'l',
                          "All seeds with size smaller than 'minimal seed length' are discarded.", SEEDING_PARAMETERS,
                          16, checkPositiveValue ),
          xMinimalSeedAmbiguity(
              this, "Minimal Ambiguity",
              "During the extension of seeds using the FMD-index: With increasing extension width, the number of "
              "occurrences of corresponding seeds on the reference monotonically decreases. Keep extending, while "
              "the "
              "number of occurrences is higher than 'Minimal Ambiguity'.",
              SEEDING_PARAMETERS, 0, checkPositiveValue ),
          xMaximalSeedAmbiguity( this, "Maximal Ambiguity",
                                 "Discard seeds that occur more than 'Maximal ambiguity' time on the reference. Set "
                                 "this option to zero in order to disable it.",
                                 SEEDING_PARAMETERS, 100, checkPositiveValue ),
          xSkipAmbiguousSeeds( this, "Skip Ambiguous Seeds",
                               "Enabled: Discard all seeds that are more ambiguous than <Maximal Ambiguity>. Disabled: "
                               "sample <Maximal Ambiguity> random seeds from too ambiguous seeds.",
                               SEEDING_PARAMETERS, false ),
          xMinimalSeedSizeDrop( this, "Seeding Drop-off A - Minimal Seed Size",
                                "Heuristic runtime optimization: For a given read R, let N be the number of seeds of "
                                "size >= <val>. Discard R, if N < <length(R)> * <Seeding drop-off B>.",
                                SEEDING_PARAMETERS, 15, checkPositiveValue ),
          xRelMinSeedSizeAmount( this, "Seeding Drop-off B - Factor",
                                 "Heuristic runtime optimization: Factor for seed drop-off calculation. For more "
                                 "information see the parameter Seeding drop-off A. ",
                                 SEEDING_PARAMETERS, 0.005, checkPositiveDoubleValue ),

          // SoC:
          xMaxNumSoC( this, "Maximal Number of SoCs", 'N',
                      "Consider the <val> best scored SoCs merely. 0 = Consider all SoCs.", SOC_PARAMETERS, 30,
                      checkPositiveValue ),
          xMinNumSoC( this, "Minimal Number of SoCs", 'M',
                      "Always consider the first <val> SoCs no matter the Heuristic optimizations. Increasing this "
                      "parameter might improve the quality of supplementary alignments.",
                      SOC_PARAMETERS, 1, checkPositiveValue ),
          xSoCWidth( this, "Fixed SoC Width",
                     "Set the SoC width to a fixed value. 0 = use the formula given in the paper. (for debugging "
                     "purposes.)",
                     SOC_PARAMETERS, 0, checkPositiveValue ),
          xRectangularSoc( this, "Rectangular SoC",
                           "Is the SoC Rectangular shape or parallelogram shape. A rectangular shape can help with "
                           "finding inversions, but struggles more with seed-noise.",
                           SOC_PARAMETERS, true ),

          // SAM
          xReportN( this, "Maximal Number of Reported Alignments", 'n',
                    "Do not output more than <val> alignments. Set to zero for unlimited output.", SAM_PARAMETERS, 0,
                    checkPositiveValue ),
          xMinAlignmentScore( this, "Minimal Alignment Score",
                              "Suppress the output of alignments with a score below <val>.", SAM_PARAMETERS, 75 ),
          xNoSecondary( this, "Omit Secondary Alignments", "Suppress the output of secondary alignments.",
                        SAM_PARAMETERS, false ),
          xNoSupplementary( this, "Omit Supplementary Alignments", "Suppress the output of supplementary alignments.",
                            SAM_PARAMETERS, false ),
          xMaxOverlapSupplementary( this, "Maximal Supplementary Overlap",
                                    "A non-primary alignment A is considered supplementary, if less than <val> "
                                    "percent of A overlap with the "
                                    "primary alignment on the query. Otherwise A is considered secondary.",
                                    SAM_PARAMETERS, 0.1, checkPositiveDoubleValue ),
          xMaxSupplementaryPerPrim( this, "Number Supplementary Alignments",
                                    "Maximal Number of supplementary alignments per primary alignment.", SAM_PARAMETERS,
                                    1, checkPositiveValue ),
          xEmulateNgmlrTags(
              this, "Emulate NGMLR's tag output",
              "Output SAM tags as NGMLR would. Enable this flag if you want to use MA in combination with "
              "Sniffles. Enabling this flag will drastically increase the size of the SAM output file.",
              SAM_PARAMETERS, false ),
          xOutputMCigar( this, "Use M in CIGAR",
                         "Disabled: Distinguish matches and mismatches in CIGARs using '=' and 'X' operations. "
                         "Enabled: Use the 'M' operation in CIGARs.",
                         SAM_PARAMETERS, true ),
          xCGTag( this, "Output long cigars in CG tag",
                  "Some programs crash, if cigars become too long. If this flag is enabled, the CG:B:I tag is used for "
                  "the output of long cigars (cigars with more than 65536 operations).",
                  SAM_PARAMETERS, true ),
          xSoftClip( this, "Soft clip",
                     "Output the full query for each alignment, instead of omitting the sequence before and after the "
                     "alignment.",
                     SAM_PARAMETERS, false ),

          // SV
          xMaxDeltaDistanceInCLuster( this, "Maximal distance between clusters",
                                      "Maximal distance between the deltas of the closest two seeds in a cluster",
                                      SV_PARAMETERS, 200, checkPositiveValue ),
          xPaddingReSeeding( this, "re seeding padding", "padding for re seeding", SV_PARAMETERS, 100,
                             checkPositiveValue ),
          // new SV
          xMinSeedSizeSV(
              this, "Minimal Seed Size SV",
              "All seeds with size smaller than 'minimal seed length' are discarded.",
              SV_PARAMETERS, 18, checkPositiveValue ),
          xMaxAmbiguitySv( this, "Maximal Ambiguity SV",
                           "Discard seeds that occur more than 'Maximal ambiguity' time on the reference. Set "
                           "this option to zero in order to disable it."
                           "Only relevant for MSV.",
                           SV_PARAMETERS, 10000, checkPositiveValue ),
          xDoDummyJumps( this, "Do Dummy Jumps",
                         "Compute edges that connect to the sentinel vertex.", SV_PARAMETERS,
                         true ),
          xMinDistDummy(
              this, "Minimal Dummy Distance",
              "Keep edges that connect to the sentinel vertex only if their seed is at least <val> nt away from "
              "the read's ends.",
              SV_PARAMETERS, 50, checkPositiveValue ),
          xMaxDistDummy(
              this, "Maximal Dummy Distance",
              "Keep edges that connect to the sentinel vertex only if their seed is less than <val> nt away from "
              "the read's ends.",
              SV_PARAMETERS, 60, checkPositiveValue ),
          xRevCompPairedReadMates(
              this, "Paired Mate - Mate Pair",
              "Reverse complement sequences of mates while inserting into DB.", SV_PARAMETERS,
              true ),
          xMaxSizeReseed( this, "Max Size Reseed",
                          "Reseeding rectangles are limited to a width and height of <val>.",
                          SV_PARAMETERS, 50, checkPositiveValue ),
          xMinSizeEdge(
              this, "Min Size Edge",
              "Discard all matrix entries that connect breakends closer than <val> nt.",
              SV_PARAMETERS, 0, checkPositiveValue ),
          xMaxSizeEdge( this, "Max Size Edge",
                        "Discard all matrix entries that connect breakends further apart than <val> nt.  Set this "
                        "option to zero in order to disable it.",
                        SV_PARAMETERS, 0, checkPositiveValue ),
          xMaxFuzzinessFilter( this, "Max Fuzziness Filter", "- unused -", SV_PARAMETERS, 50, checkPositiveValue ),
          xMaxSuppNtShortCallFilter( this, "Max Supp Nt", "- unused -", SV_PARAMETERS, 10, checkPositiveValue ),
          xMaxCallSizeShortCallFilter( this, "Max Call Size Filter", "- unused -", SV_PARAMETERS, 20,
                                       checkPositiveValue ),
          xMaxRegionSizeForAmbiguityFilter( this, "Max Regionsize Ambiguity Filter", "", SV_PARAMETERS, 200,
                                       checkPositiveValue ),
          xMaxRefAmbiguityJump( this, "Max Ref Ambiguity Jump", "- unused -", SV_PARAMETERS, 10, checkPositiveValue ),
          xMMFilterMaxOcc(
              this, "Max Occ MM Filter",
              "Maximal number of occurrences for a k-mer during the filtering on reads.",
              SV_PARAMETERS, 200, checkPositiveValue ),
          xMinNtInSoc( this, "Min NT in SoC",
                       "Discard all SoCs with a lower accumulative size of seeds than <val> (before reseeding). Only "
                       "relevant for MSV.",
                       SV_PARAMETERS, 150, checkPositiveValue ),
          xMinNtInSocRelative( this, "Relative Min NT in SoC",
                       "Discard all SoCs with a lower accumulative size of seeds than <val>*read_size "
                       " (before reseeding). SoCs will only be discarded if this condition and Min NT in SoC both are met.",
                       SV_PARAMETERS, 1, checkPositiveValue ),
          xMinNtAfterReseeding( this, "Min NT after reseeding",
                                "Discard all SoCs with a lower accumulative size of seeds than <val> (after "
                                "reseeding).",
                                SV_PARAMETERS, 100, checkPositiveValue ),
          xMinNtAfterReseedingRelative( this, "Relative Min NT after reseeding",
                                "Discard all SoCs with a lower accumulative size of seeds than <val> (after "
                                "reseeding). SoCs will only be discarded if this condition and Min NT after reseeding both are met.",
                                SV_PARAMETERS, 1, checkPositiveValue ),
          xMinReadsInCall( this, "Min Reads in call",
                           "Keep all SV-calls that are supported by at least <val> reads.",
                           SV_PARAMETERS, 2, checkPositiveValue ),

          // Heuristic
          xSoCScoreDecreaseTolerance( this, "SoC Score Drop-off",
                                      "Let x be the maximal encountered SoC score. Stop harmonizing SoCs if there is "
                                      "a SoC with a score smaller than <val>*x.",
                                      HEURISTIC_PARAMETERS, 0.1, checkPositiveDoubleValue ),
          xHarmScoreMin( this, "Minimal Harmonization Score",
                         "Discard all harmonized SoCs with scores smaller than <val>. "
                         "Only keep detected inversions with a score >= <val> * <Match Score>.",
                         HEURISTIC_PARAMETERS, 18, checkPositiveValue ),
          xHarmScoreMinRel( this, "Relative Minimal Harmonization Score",
                            "Discard all harmonized SoCs with scores smaller than length(read)*<val>.",
                            HEURISTIC_PARAMETERS, 0.002, checkPositiveDoubleValue ),
          xScoreDiffTolerance(
              this, "Harmonization Drop-off A - Score Difference",
              "Let x be the maximal encountered harmonization score. Stop harmonizing further SoCs, if "
              "<Harmonization "
              "Drop-off B> many SoCs with scores below x - <readlength> * <val> occur consecutively.",
              HEURISTIC_PARAMETERS, 0.0001, checkPositiveDoubleValue ),
          xMaxScoreLookahead( this, "Harmonization Drop-off B - Lookahead", "See Harmonization Drop-off A.",
                              HEURISTIC_PARAMETERS, 3, checkPositiveValue ),
          xSwitchQlen( this, "Harmonization Score Drop-off - Minimal Query Length",
                       "For reads of length >= <val>: Ignore all SoCs with harmonization scores smaller than the "
                       "current maximal score. 0 = disabled.",
                       HEURISTIC_PARAMETERS, 800, checkPositiveValue ),
          xMaxDeltaDist( this, "Artifact Filter A - Maximal Delta Distance",
                         "Filter a seed, if the difference between the delta distance to its predecessor and successor "
                         "is less then <val> percent (set to 1 to disable filter) and the delta distance to its pre- "
                         "and successor is more than <Artifact Filter B> nt.",
                         HEURISTIC_PARAMETERS, 0.1, checkPositiveDoubleValue ),
          xMinDeltaDist( this, "Artifact Filter B - Minimal Delta Distance", "See Artifact Filter A",
                         HEURISTIC_PARAMETERS, 16, checkPositiveValue ),
          xDisableGapCostEstimationCutting(
              this, "Pick Local Seed Set A - Enabled",
              "Enable this flag for local seed set computation. (See Pick_Local_Seed_Set_B)", HEURISTIC_PARAMETERS,
              false ),
          xOptimisticGapCostEstimation(
              this, "Pick Local Seed Set B - Optimistic Gap Estimation",
              "After the harmonization MA checks weather it is possible to compute a positively scored alignment "
              "from "
              "the seed set. Gaps between seeds can be estimated in two ways: Optimistic [true]: Assume that the "
              "gap "
              "can be filled using merely matches and a single insertion/deletion. Pessimistic [false]: Assume "
              "that "
              "the gap can be filled using matches and mismatches that add up to a score of 0 and a single "
              "insertion/deletion.",
              HEURISTIC_PARAMETERS, true ),
          xMaxGapArea(
              this, "Maximal Gap Size",
              "If the gap between seeds is larger than <val> on query or reference, a dual extension "
              "process is used for filling the gap. Dual extension is more expensive, if the extension does not "
              "Z-drop, but more efficient otherwise.",
              HEURISTIC_PARAMETERS, 20, checkPositiveValue ),
          xGenomeSizeDisable( this, "Minimum Genome Size for Heuristics",
                              "Some heuristics can only be applied on genomes of sufficient size. The parameter "
                              "disables the SoC score "
                              "Drop-off, if the genome is shorter than <val>.",
                              HEURISTIC_PARAMETERS, 10000000, checkPositiveValue ),
          xDisableHeuristics( this, "Disable All Heuristics",
                              "Disables all runtime heuristics. (For debugging purposes)", HEURISTIC_PARAMETERS,
                              false ),
          xMinimizerK( this, "Minimizers - k", "Sequence size of minimizers", MINIMIZER_PARAMETERS, 15 ),
          xMinimizerW( this, "Minimizers - w", "Window size of minimizers", MINIMIZER_PARAMETERS, 10 ),
          xMinimizerFlag( this, "Minimizers - flag", "Flags for minimizers", MINIMIZER_PARAMETERS, 0 ),
          xMinimizerBucketBits( this, "Minimizers - bucket_bits", "@todo", MINIMIZER_PARAMETERS, 14 ),
          xMinimizerMiniBatchSize( this, "Minimizers - mini_batch_size", "@todo", MINIMIZER_PARAMETERS, 50000000 ),
          xMinimizerBatchSize( this, "Minimizers - batch_size", "@todo", MINIMIZER_PARAMETERS, 4000000000ULL )
    {
        xMeanPairedReadDistance->fEnabled = [ this ]( void ) { return this->xUsePairedReads->get( ) == true; };
        xStdPairedReadDistance->fEnabled = [ this ]( void ) { return this->xUsePairedReads->get( ) == true; };
        xPairedBonus->fEnabled = [ this ]( void ) { return this->xUsePairedReads->get( ) == true; };
        xZDropInversion->fEnabled = [ this ]( void ) { return this->xSearchInversions->get( ) == true; };
    } // constructor

    Presetting( ) : Presetting( "Unnamed" )
    {} // default constructor

    /* Named copy Constructor */
    Presetting( const Presetting& rxOtherSet, const std::string& sName ) : Presetting( rxOtherSet.sName )
    {
        this->mirror( rxOtherSet );
    } // copy constructor

    virtual std::string getSetDescription( )
    {
        return "Part of the " + sName + " presetting.";
    } // method

    /* True if the parameter-set uses paired reads */
    bool usesPairedReads( void )
    {
        return xUsePairedReads->value;
    } // method

    /* Clang requires a virtual destructor for classes comprising virtual methods */
    virtual ~Presetting( )
    {} // destructor
}; // class

template <typename VALUE_TYPE> void AlignerParameterPointer<VALUE_TYPE>::do_register( ParameterSetBase* pPresetting )
{
    if( pPresetting != nullptr )
        pPresetting->registerParameter( this->pContent );
} // method
template <typename VALUE_TYPE>
std::string AlignerParameterPointer<VALUE_TYPE>::get_desc( ParameterSetBase* pPresetting )
{
    if( pPresetting != nullptr )
        return pPresetting->getSetDescription( );
    return "";
} // method

/* Parameter which occur only once.
 * (Not sequencing technique / configuration related parameters)
 */
class GeneralParameter : public ParameterSetBase
{
  public:
    AlignerParameterPointer<AlignerParameterBase::ChoicesType>
        xSAMOutputTypeChoice; // SAM Output in the same folder as the reads.
    AlignerParameterPointer<fs::path> xSAMOutputPath; // folder path
    AlignerParameterPointer<fs::path> xSAMOutputFileName; // folder path
    AlignerParameterPointer<bool> pbUseMaxHardareConcurrency; // Exploit all cores
    AlignerParameterPointer<int> piNumberOfThreads; // selected number of threads
    AlignerParameterPointer<bool> pbPrintHelpMessage; // Print the help message to stdout
    AlignerParameterPointer<bool> pbPrintVerstion; // Print the version to stdout

    /* Constructor */
    GeneralParameter( )
        : xSAMOutputTypeChoice( this, "SAM File output", "Select output type for sam file.", GENERAL_PARAMETER,
                                AlignerParameterBase::ChoicesType{ { "Read_Folder", "In Read Folder" },
                                                                   { "Specified_Folder", "In Specified Folder" },
                                                                   { "Specified_File", "As Specified File" } } ),
          xSAMOutputPath( this, "Folder for SAM Files",
                          "Folder for SAM output in the case that the output is not directed to the reads' folder.",
                          GENERAL_PARAMETER, fs::temp_directory_path( ) ),
          xSAMOutputFileName( this, "SAM File name", 'o',
                              "Name of the SAM file that is used for the output of alignments.", GENERAL_PARAMETER,
                              "ma_out.sam" ),
          pbUseMaxHardareConcurrency( this, "Use all Processor Cores",
                                      "Number of threads used for alignments is identical to the number "
                                      "of processor cores.",
                                      GENERAL_PARAMETER, true ),
          piNumberOfThreads( this, "Number of Threads", 't',
                             "Number of threads used in the context of alignments. This options is only available, if "
                             "'use all processor cores' is off.",
                             GENERAL_PARAMETER, 1, checkPositiveValue ),
          pbPrintHelpMessage( this, "Help", 'h', "Print the complete help text.", GENERAL_PARAMETER, false ),
          pbPrintVerstion( this, "Version", 'v', "Print the version and exit.", GENERAL_PARAMETER, false )
    {
        xSAMOutputPath->fEnabled = [ this ]( void ) { return this->xSAMOutputTypeChoice->uiSelection == 1; };
        xSAMOutputFileName->fEnabled = [ this ]( void ) { return this->xSAMOutputTypeChoice->uiSelection == 2; };
        piNumberOfThreads->fEnabled = [ this ]( void ) { return this->pbUseMaxHardareConcurrency->get( ) == false; };
    } // constructor

    /* Named copy Constructor */
    GeneralParameter( const GeneralParameter& rxOtherSet ) : GeneralParameter( )
    {
        mirror( rxOtherSet );
    } // copy constructor

    virtual std::string getSetDescription( )
    {
        return "Independent of presettings.";
    } // method

    size_t getNumThreads( ) const
    {
        size_t uiConcurency = std::thread::hardware_concurrency( );
        if( !this->pbUseMaxHardareConcurrency->get( ) )
        {
            uiConcurency = this->piNumberOfThreads->get( );
        } // if
        return uiConcurency;
    } // method

    /* Clang requires a virtual destructor for classes comprising virtual methods */
    virtual ~GeneralParameter( )
    {} // destructor
}; // class


/** @brief master class for Global Parameters.
 *  @details
 *  there will only ever be one instance of this class.
 */
class GlobalParameter : public ParameterSetBase
{
  public:
    // jump parameters
    AlignerParameterPointer<double> xJumpS; // s parameter for fuzziness of jumps
    AlignerParameterPointer<double> xJumpSNeg; // s parameter for fuzziness of jumps in negative direciton
    AlignerParameterPointer<double> xJumpM; // m parameter for fuzziness of jumps
    AlignerParameterPointer<double> xJumpH; // h parameter for fuzziness of jumps
    AlignerParameterPointer<int> xSeedDirFuzziness; // m parameter for fuzziness of jumps
    // parameters from default_parameters @todo these should all be checked and removed if possible...
    AlignerParameterPointer<int> iMatch;
    AlignerParameterPointer<int> iMissMatch;
    AlignerParameterPointer<int> iGap;
    AlignerParameterPointer<int> iExtend;
    AlignerParameterPointer<int> iGap2;
    AlignerParameterPointer<int> iExtend2;
    AlignerParameterPointer<int> xCoverageBinSize; // Coveragebin size
    AlignerParameterPointer<int> uiSVPenalty;

    /* Constructor */
    GlobalParameter( )
        : xJumpS( this, "fuzziness-s", "Maximal fuzziness for entries.", SV_PARAMETERS, 400 ),
          xJumpSNeg( this, "fuzziness-s-neg", "Maximal fuzziness for entries.", SV_PARAMETERS,
                     400 ),
          xJumpM( this, "fuzziness-m", "Fuzziness slope.", SV_PARAMETERS, 0.25 ),
          xJumpH( this, "fuzziness-h", "Fuzziness zero-point.", SV_PARAMETERS, 25 ),
          xSeedDirFuzziness( this, "Seed Dir Fuzziness", "Absolute fuzziness in seed direction.", SV_PARAMETERS, 3,
                             checkPositiveValue ),
          // DP:
          iMatch( this, "Match Score",
                  "Match score. (Used in the context of Dynamic Programming and for SoC width computation.)",
                  DP_PARAMETERS, 2, checkPositiveValue ),
          iMissMatch( this, "Mismatch Penalty", "Penalty for mismatch. ", DP_PARAMETERS, 4, checkPositiveValue ),
          iGap( this, "Gap penalty", "First penalty for gap opening. (Two piece affine gap costs)", DP_PARAMETERS, 4,
                checkPositiveValue ),
          iExtend( this, "Extend Penalty", "First penalty for gap extension.  (Two piece affine gap costs)",
                   DP_PARAMETERS, 2, checkPositiveValue ),
          iGap2( this, "Second Gap Penalty", "Second penalty for gap opening. (Two piece affine gap costs)",
                 DP_PARAMETERS, 24, checkPositiveValue ),
          iExtend2( this, "Second Extend Penalty", "Second penalty for gap extension. (Two piece affine gap costs)",
                    DP_PARAMETERS, 1, checkPositiveValue ),
          xCoverageBinSize( this, "Coveragebin size",
                           "Binsize for storing the seed coverage on the genome",
                           SV_PARAMETERS, 100, checkPositiveValue ),
          uiSVPenalty( this, "Pick Local Seed Set C - Maximal Gap Penalty",
                       "Maximal gap cost penalty during local seed set computation.", HEURISTIC_PARAMETERS, 100,
                       checkPositiveValue )
    {} // constructor

    /* Named copy Constructor */
    GlobalParameter( const GlobalParameter& rxOtherSet ) = delete;

    virtual std::string getSetDescription( )
    {
        return "Independent of presettings.";
    } // method

    /* Clang requires a virtual destructor for classes comprising virtual methods */
    virtual ~GlobalParameter( )
    {} // destructor
}; // class


/** @brief Single master object for Global Parameters. */
const std::shared_ptr<GlobalParameter> pGlobalParams = std::make_shared<GlobalParameter>( );

/* Management of single global parameter-set and all config (preset) parameter-sets */
class ParameterSetManager
{
  private:
    // Pointer to the currently selected parameter-set
    std::shared_ptr<Presetting> pSelectedParamSet = NULL;

  public:
    std::shared_ptr<GeneralParameter> pGeneralParameterSet;

    // Presets for aligner configuration
    std::map<std::string, std::shared_ptr<Presetting>> xParametersSets;

    void setCommonSvParameters(std::shared_ptr<Presetting> pPre)
    {
        pPre->xMaximalSeedAmbiguity->set( 1 );
        pPre->xRectangularSoc->set( false );
    }

    ParameterSetManager( ) : pGeneralParameterSet( std::make_shared<GeneralParameter>( ) )
    {
        xParametersSets.emplace( "default", std::make_shared<Presetting>( "Default" ) );

        xParametersSets.emplace( "illumina", std::make_shared<Presetting>( "Illumina" ) );
        xParametersSets[ "illumina" ]->xSeedingTechnique->set( 1 ); // use SMEMs
        xParametersSets[ "illumina" ]->xMaximalSeedAmbiguity->set( 500 );
        xParametersSets[ "illumina" ]->xMinNumSoC->set( 10 );
        xParametersSets[ "illumina" ]->xMaxNumSoC->set( 20 );

        xParametersSets.emplace( "illuminapaired", std::make_shared<Presetting>( "Illumina Paired" ) );
        xParametersSets[ "illuminapaired" ]->xSeedingTechnique->set( 1 ); // use SMEMs
        xParametersSets[ "illuminapaired" ]->xUsePairedReads->set( true );
        xParametersSets[ "illuminapaired" ]->xMaximalSeedAmbiguity->set( 500 );
        xParametersSets[ "illuminapaired" ]->xMinNumSoC->set( 10 );
        xParametersSets[ "illuminapaired" ]->xMaxNumSoC->set( 20 );

        xParametersSets.emplace( "pacbio", std::make_shared<Presetting>( "PacBio" ) );
        xParametersSets[ "pacbio" ]->xMaxSupplementaryPerPrim->set( 100 );
        xParametersSets[ "pacbio" ]->xMinNumSoC->set( 5 );


        xParametersSets.emplace( "nanopore", std::make_shared<Presetting>( "Nanopore" ) );
        xParametersSets[ "nanopore" ]->xSeedingTechnique->set( 1 );
        xParametersSets[ "nanopore" ]->xMaxSupplementaryPerPrim->set( 100 );
        xParametersSets[ "nanopore" ]->xMinNumSoC->set( 5 );

        xParametersSets.emplace( "sv-illumina", std::make_shared<Presetting>( "SV-Illumina" ) );
        setCommonSvParameters( xParametersSets[ "sv-illumina" ] );
        xParametersSets[ "sv-illumina" ]->xDoDummyJumps->set( false );
        xParametersSets[ "sv-illumina" ]->xMinReadsInCall->set( 10 );
        //xParametersSets[ "sv-illumina" ]->xMaxSizeEdge->set( 200 );
        xParametersSets[ "sv-illumina" ]->xHarmScoreMin->set( 10 );
        xParametersSets[ "sv-illumina" ]->xHarmScoreMinRel->set( 0.1 );
        xParametersSets[ "sv-illumina" ]->xMinNtAfterReseedingRelative->set( 0.4 );
        xParametersSets[ "sv-illumina" ]->xMinNtInSoc->set( 25 );
        xParametersSets[ "sv-illumina" ]->xMinNtInSocRelative->set( 0.1 );


        xParametersSets.emplace( "sv-pacbio", std::make_shared<Presetting>( "SV-PacBio" ) );
        setCommonSvParameters( xParametersSets[ "sv-pacbio" ] );
        xParametersSets[ "sv-pacbio" ]->xSoCWidth->set( 3000 );
        xParametersSets[ "sv-pacbio" ]->xMaxSizeReseed->set( 1000 );
        //xParametersSets[ "sv-pacbio" ]->xMinSizeEdge->set( 200 );
        xParametersSets[ "sv-pacbio" ]->xMinNtAfterReseeding->set( 600 );
        xParametersSets[ "sv-pacbio" ]->xHarmScoreMinRel->set( 0 );
        xParametersSets[ "sv-pacbio" ]->xHarmScoreMin->set( 25 );
        xParametersSets[ "sv-pacbio" ]->xMinNtInSoc->set( 25 );

        // Initially select Illumina
        this->pSelectedParamSet = xParametersSets[ "default" ];
    } // constructor

    ParameterSetManager( const ParameterSetManager& ) = delete; // delete copy constructor

    std::shared_ptr<Presetting> get( const std::string& sKey )
    {
        if( xParametersSets.count( ParameterSetBase::uniqueParameterName( sKey ) ) == 0 )
            throw std::runtime_error( "The presetting '" + sKey + "' can not be found." );
        return xParametersSets[ ParameterSetBase::uniqueParameterName( sKey ) ];
    } // method

    /* Set selected parameter set by using a key.
     */
    void setSelected( const std::string& sKey )
    {
        if( xParametersSets.count( ParameterSetBase::uniqueParameterName( sKey ) ) == 0 )
            throw std::runtime_error( "The presetting '" + sKey + "' can not be found." );
        this->pSelectedParamSet = xParametersSets[ ParameterSetBase::uniqueParameterName( sKey ) ];
    } // method

    /* Delivers pointer to selected parameter-set */
    std::shared_ptr<Presetting> getSelected( void )
    {
        return pSelectedParamSet;
    } // method

    /** required for exporting get to python do not use in cpp code */
    std::shared_ptr<Presetting> getSelected_py( void )
    {
        return this->getSelected( );
    } // method

    /* Delivers pointer to selected parameter-set */
    const std::shared_ptr<Presetting> getSelected( void ) const
    {
        return pSelectedParamSet;
    } // method

    std::shared_ptr<AlignerParameterBase> byName( const std::string& rParameterName )
    {
        if( pGlobalParams->hasName( rParameterName ) )
            return pGlobalParams->byName( rParameterName );
        if( pGeneralParameterSet->hasName( rParameterName ) )
            return pGeneralParameterSet->byName( rParameterName );
        if( this->getSelected( )->hasName( rParameterName ) )
            return this->getSelected( )->byName( rParameterName );
        throw std::runtime_error( std::string( "Could not find parameter: " ).append( rParameterName ) );
    } // method

    std::shared_ptr<AlignerParameterBase> byShort( const char cX )
    {
        if( pGlobalParams->hasShort( cX ) )
            return pGlobalParams->byShort( cX );
        if( pGeneralParameterSet->hasShort( cX ) )
            return pGeneralParameterSet->byShort( cX );
        if( this->getSelected( )->hasShort( cX ) )
            return this->getSelected( )->byShort( cX );
        throw std::runtime_error( std::string( "Could not find parameter: " ) += cX );
    } // method

    size_t getNumThreads( ) const
    {
        return pGeneralParameterSet->getNumThreads( );
    } // method

    void addSetting( const std::string& rsName )
    {
        xParametersSets.emplace( rsName, std::make_shared<Presetting>( rsName ) );
    } // method
}; // class


// parameter set manager is exported from export.cpp....
