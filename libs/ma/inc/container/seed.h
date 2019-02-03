/**
 * @file seed.h
 * @brief Implements Seed.
 * @author Markus Schmidt
 */
#ifndef SEED_H
#define SEED_H

#include "container/container.h"
#include "container/interval.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <list>
/// @endcond

#define DELTA_CACHE ( 1 )
#define CONTIG_ID_CACHE ( 0 )

namespace libMA
{
///@brief any index on the query or reference nucleotide sequence is given in this datatype
typedef uint64_t nucSeqIndex;

/**
 * @brief A seed.
 * @details
 * A extracted seed, that comprises two intervals, one on the query one on the reference.
 * Both intervals are equal in size.
 * @note the overloaded functions of Interval refer to the Interval on the query.
 * @ingroup container
 */
class Seed : public Container, public Interval<nucSeqIndex>
{
  public:
    ///@brief the beginning of the match on the reference
    nucSeqIndex uiPosOnReference;
    unsigned int uiAmbiguity;
#if DELTA_CACHE == ( 1 )
    nucSeqIndex uiDelta = 0;
#endif
#if CONTIG_ID_CACHE == ( 1 )
    size_t uiContigId = 0;
#endif

    /**
     * @brief Creates a new Seed.
     */
    Seed( const nucSeqIndex uiPosOnQuery, const nucSeqIndex uiLength,
          const nucSeqIndex uiPosOnReference )
        : Interval( uiPosOnQuery, uiLength ), uiPosOnReference( uiPosOnReference ), uiAmbiguity( 0 )
    {} // constructor

    /**
     * @brief Creates a new Seed.
     */
    Seed( const nucSeqIndex uiPosOnQuery, const nucSeqIndex uiLength,
          const nucSeqIndex uiPosOnReference, const unsigned int uiAmbiguity )
        : Interval( uiPosOnQuery, uiLength ),
          uiPosOnReference( uiPosOnReference ),
          uiAmbiguity( uiAmbiguity )
    {} // constructor

    /**
     * @brief Copys from a Seed.
     */
    Seed( const Seed &rOther )
        : Interval( rOther ),
          uiPosOnReference( rOther.uiPosOnReference ),
          uiAmbiguity( rOther.uiAmbiguity )
#if DELTA_CACHE == ( 1 )
          ,
          uiDelta( rOther.uiDelta )
#endif
#if CONTIG_ID_CACHE == ( 1 )
          ,
          uiContigId( rOther.uiContigId )
#endif
    {} // copy constructor

    /**
     * @brief Default Constructor.
     */
    Seed( ) : Interval( )
    {} // default constructor

    /**
     * @brief Returns the beginning of the seed on the reference.
     */
    inline nucSeqIndex start_ref( ) const
    {
        return uiPosOnReference;
    } // function

    /**
     * @brief Returns the end of the seed on the reference.
     */
    inline nucSeqIndex end_ref( ) const
    {
        return uiPosOnReference + size( );
    } // function

    /**
     * @brief Returns the value of the seed.
     * @details
     * A seeds value corresponds to its size.
     */
    inline nucSeqIndex getValue( ) const
    {
        return size( );
    } // function

    /**
     * @brief Copys from another Seed.
     */
    inline Seed &operator=( const Seed &rxOther )
    {
        Interval::operator=( rxOther );
        uiPosOnReference = rxOther.uiPosOnReference;
        uiAmbiguity = rxOther.uiAmbiguity;
#if DELTA_CACHE == ( 1 )
        uiDelta = rxOther.uiDelta;
#endif
#if CONTIG_ID_CACHE == ( 1 )
        uiContigId = rxOther.uiContigId;
#endif
        return *this;
    } // operator

    /*
     * @brief compares two Seeds.
     * @returns true if start and size are equal, false otherwise.
     */
    inline bool operator==( const Seed &rxOther )
    {
        return Interval::operator==( rxOther ) && uiPosOnReference == rxOther.uiPosOnReference;
    } // operator

    // overload
    inline bool canCast( const std::shared_ptr<Container> &c ) const
    {
        return std::dynamic_pointer_cast<Seed>( c ) != nullptr;
    } // function

    // overload
    inline std::string getTypeName( ) const
    {
        return "Seed";
    } // function

    // overload
    inline std::shared_ptr<Container> getType( ) const
    {
        return std::shared_ptr<Container>( new Seed( ) );
    } // function

}; // class

class Alignment;
/**
 * @brief Used to store some statistics to each alignment
 * @details
 * Intended for figuring out optimal thresholds.
 */
class AlignmentStatistics
{
  public:
    unsigned int index_of_strip;
    unsigned int num_seeds_in_strip;
    unsigned int anchor_size;
    unsigned int anchor_ambiguity;
    std::weak_ptr<Alignment> pOther;
    bool bFirst;
    bool bSetMappingQualityToZero = false;
    std::string sName;
    nucSeqIndex uiInitialQueryBegin;
    nucSeqIndex uiInitialRefBegin;
    nucSeqIndex uiInitialQueryEnd;
    nucSeqIndex uiInitialRefEnd;

    AlignmentStatistics( )
        : index_of_strip( 0 ),
          num_seeds_in_strip( 0 ),
          anchor_size( 0 ),
          anchor_ambiguity( 0 ),
          pOther( ),
          bFirst( false ),
          sName( "unknown" ),
          uiInitialQueryBegin( 0 ),
          uiInitialRefBegin( 0 ),
          uiInitialQueryEnd( 0 ),
          uiInitialRefEnd( 0 )
    {}

    // void operator=( const AlignmentStatistics &rOther )
    // {
    //     index_of_strip = rOther.index_of_strip;
    //     num_seeds_in_strip = rOther.num_seeds_in_strip;
    //     anchor_size = rOther.anchor_size;
    //     anchor_ambiguity = rOther.anchor_ambiguity;
    //     pOther = rOther.pOther;
    //     bFirst = rOther.bFirst;
    //     sName = rOther.sName;
    //     uiInitialQueryBegin = rOther.uiInitialQueryBegin;
    //     uiInitialRefBegin = rOther.uiInitialRefBegin;
    //     uiInitialQueryEnd = rOther.uiInitialQueryEnd;
    //     uiInitialRefEnd = rOther.uiInitialRefEnd;
    // } // function
}; // class

class SVInfo
{
  public:
    /**
     * @details
     * Contains indices of seeds, where a structural variant (deletion or insertion)
     * follows the seed.
     */
    std::vector<size_t> vSeedIndicesOfSVIndels;
}; // class

/**
 * @brief A set with Seed elements.
 * @details
 * Also holds the summed up score of the seeds within the list.
 * @ingroup Container
 */
class Seeds : public Container
{
  private:
    typedef std::vector<Seed> TP_VEC;
    ///@brief the content of the seed set
    TP_VEC vContent;

  public:
    typedef typename TP_VEC::value_type value_type;
    typedef typename TP_VEC::size_type size_type;
    typedef typename TP_VEC::difference_type difference_type;
    typedef typename TP_VEC::iterator iterator;

    nucSeqIndex mem_score = 0;
    // some statistics
    AlignmentStatistics xStats;
    /// @brief is the seed set consistent (set to true after harmonization)
    bool bConsistent = false;

    Seeds( std::initializer_list<value_type> init ) : vContent( init )
    {} // initializer list constructor

    template <class InputIt> Seeds( InputIt xBegin, InputIt xEnd ) : vContent( xBegin, xEnd )
    {} // iterator constructor

    Seeds( ) : vContent( )
    {} // default constructor

    Seeds( const std::shared_ptr<Seeds> pOther ) : vContent( pOther->vContent )
    {} // constructor

    Seeds( size_t numElements ) : vContent( numElements )
    {} // constructor

    // overload
    inline bool canCast( std::shared_ptr<Container> c ) const
    {
        return std::dynamic_pointer_cast<Seeds>( c ) != nullptr;
    } // method

    // overload
    inline std::string getTypeName( ) const
    {
        return "Seeds";
    } // method

    // overload
    inline std::shared_ptr<Container> getType( ) const
    {
        return std::shared_ptr<Container>( new Seeds( ) );
    } // method

    /// @brief returns the sum off all scores within the list
    inline nucSeqIndex getScore( ) const
    {
        nucSeqIndex iRet = 0;
        for( const Seed &rS : *this )
            iRet += rS.getValue( );
        return iRet;
    } // function

    /// @brief append another seed set
    inline void append( const std::shared_ptr<Seeds> pOther )
    {
        for( Seed &rS : pOther->vContent )
            push_back( rS );
    } // method

    /// @briefreturn wether this seed set is larger according to getScore()
    inline bool larger( const std::shared_ptr<Container> pOther ) const
    {
        const std::shared_ptr<Seeds> pSeeds = std::dynamic_pointer_cast<Seeds>( pOther );
        if( pSeeds == nullptr )
            return true;
        return getScore( ) > pSeeds->getScore( );
    } // method

    // setter
    inline value_type &operator[]( size_type uiI )
    {
        return vContent[ uiI ];
    } // operator

    // getter
    inline const value_type &operator[]( size_type uiI ) const
    {
        return vContent[ uiI ];
    } // operator

    inline void push_back( const value_type &value )
    {
        vContent.push_back( value );
    } // method

    inline void pop_back( void )
    {
        vContent.pop_back( );
    } // method

    template <class... Args> inline void emplace_back( Args &&... args )
    {
        vContent.emplace_back( args... );
    } // method

    inline size_type size( void ) const
    {
        return vContent.size( );
    } // method

    inline bool empty( void ) const
    {
        return vContent.empty( );
    } // method

    inline value_type &front( void )
    {
        return vContent.front( );
    } // method

    inline value_type &back( void )
    {
        return vContent.back( );
    } // method

    inline const value_type &front( void ) const
    {
        return vContent.front( );
    } // method

    inline const value_type &back( void ) const
    {
        return vContent.back( );
    } // method

    inline TP_VEC::iterator begin( void ) noexcept
    {
        return vContent.begin( );
    } // method

    inline TP_VEC::iterator end( void ) noexcept
    {
        return vContent.end( );
    } // method

    inline TP_VEC::const_iterator begin( void ) const noexcept
    {
        return vContent.begin( );
    } // method

    inline TP_VEC::const_iterator end( void ) const noexcept
    {
        return vContent.end( );
    } // method

    inline void erase( TP_VEC::iterator pos )
    {
        vContent.erase( pos );
    } // method

    inline void erase( TP_VEC::iterator first, TP_VEC::iterator last )
    {
        vContent.erase( first, last );
    } // method

    inline TP_VEC::iterator insert( TP_VEC::const_iterator pos, const value_type &value )
    {
        return vContent.insert( pos, value );
    } // method

    template <class InputIt>
    inline TP_VEC::iterator insert( TP_VEC::const_iterator pos, InputIt first, InputIt last )
    {
        return vContent.insert( pos, first, last );
    } // method

    void reserve( size_type uiN )
    {
        vContent.reserve( uiN );
    } // method

    void clear( ) noexcept
    {
        vContent.clear( );
    } // method
}; // class
} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief exports the Seed and Seedlist classes to python.
 * @ingroup export
 */
#ifdef WITH_BOOST
void exportSeed( );
#else
void exportSeed( py::module& rxPyModuleId );
#endif
#endif

#endif