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
#define AMBIGUITY ( 0 )

namespace libMA
{
///@brief any index on the query or reference nucleotide sequence is given in this datatype
typedef uint64_t nucSeqIndex;

class Pack;
class NucSeq;
/**
 * @brief A seed.
 * @details
 * A extracted seed, that comprises two intervals, one on the query one on the reference.
 * Both intervals are equal in size.
 * @note the overloaded functions of Interval refer to the Interval on the query.
 * @ingroup container
 */
class Seed : public Interval<nucSeqIndex>
{
  public:
    ///@brief the beginning of the match on the reference
    nucSeqIndex uiPosOnReference;
#if AMBIGUITY == ( 1 )
    unsigned int uiAmbiguity = 0;
#endif
#if DELTA_CACHE == ( 1 )
    nucSeqIndex uiDelta = 0;
#endif
#if CONTIG_ID_CACHE == ( 1 )
    size_t uiContigId = 0;
#endif

    /**
     * @brief Creates a new Seed.
     */
    Seed( const nucSeqIndex uiPosOnQuery, const nucSeqIndex uiLength, const nucSeqIndex uiPosOnReference )
        : Interval( uiPosOnQuery, uiLength ), uiPosOnReference( uiPosOnReference )
    {} // constructor

#if AMBIGUITY == ( 1 )
    /**
     * @brief Creates a new Seed.
     */
    Seed( const nucSeqIndex uiPosOnQuery, const nucSeqIndex uiLength, const nucSeqIndex uiPosOnReference,
          const unsigned int uiAmbiguity )
        : Interval( uiPosOnQuery, uiLength ), uiPosOnReference( uiPosOnReference ), uiAmbiguity( uiAmbiguity )
    {} // constructor
#endif

    /**
     * @brief Copys from a Seed.
     */
    Seed( const Seed& rOther )
        : Interval( rOther ),
          uiPosOnReference( rOther.uiPosOnReference )
#if AMBIGUITY == ( 1 )
          ,
          uiAmbiguity( rOther.uiAmbiguity )
#endif
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
    {
#if DEBUG_LEVEL == 0
        static_assert( sizeof( Seed ) == 32, "" );
#endif
    } // default constructor

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
    inline Seed& operator=( const Seed& rxOther )
    {
        Interval::operator=( rxOther );
        uiPosOnReference = rxOther.uiPosOnReference;
#if AMBIGUITY == ( 1 )
        uiAmbiguity = rxOther.uiAmbiguity;
#endif
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
    inline bool operator==( const Seed& rxOther )
    {
        return Interval::operator==( rxOther ) && uiPosOnReference == rxOther.uiPosOnReference;
    } // operator

}; // class

#if 0
    template<int s> struct CheckSizeOfSeeds;
    CheckSizeOfSeeds<sizeof(Seed)> xCheckSizeOfSeeds;
#endif

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
    bool bFirst; // for paired alignments: is the a alignment for the first mate read of the second one
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
    {} // constructor
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

    /**
     * @brief sorts by reference position (asc) and for equal reference pos by size (desc)
     */
    inline void sortByRefPos( )
    {
        std::sort( vContent.begin( ),
                   vContent.end( ),
                   []( const Seed& rA, const Seed& rB ) {
                       if( rA.start_ref( ) != rB.start_ref( ) )
                           return rA.start_ref( ) < rB.start_ref( );
                       return rA.size( ) > rB.size( );
                   } // lambda
        );
    } // method

    /**
     * @brief sorts by query position (asc)
     * @details checks wether the vector is already sorted
     */
    inline void sortByQPos( )
    {
        for( size_t uiI = 1; uiI < vContent.size( ); uiI++ )
            if( vContent[ uiI - 1 ].start( ) > vContent[ uiI ].start( ) ) // check if not sorted
            {
                // if not sorted then sort...
                std::sort( vContent.begin( ),
                           vContent.end( ),
                           []( const Seed& rA, const Seed& rB ) { return rA.start( ) < rB.start( ); } // lambda
                ); // std::sort call
                return;
            } // if
    } // method

    /// @brief returns the sum off all scores within the list
    inline nucSeqIndex getScore( ) const
    {
        nucSeqIndex iRet = 0;
        for( const Seed& rS : *this )
            iRet += rS.getValue( );
        return iRet;
    } // function

    /// @brief append another seed set
    inline void append( const std::shared_ptr<Seeds> pOther )
    {
        for( Seed& rS : pOther->vContent )
            push_back( rS );
    } // method

    // setter
    inline value_type& operator[]( size_type uiI )
    {
        return vContent[ uiI ];
    } // operator

    // getter
    inline const value_type& operator[]( size_type uiI ) const
    {
        return vContent[ uiI ];
    } // operator

    inline void push_back( const value_type& value )
    {
        vContent.push_back( value );
    } // method

    inline void pop_back( void )
    {
        vContent.pop_back( );
    } // method

    template <class... Args> inline void emplace_back( Args&&... args )
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

    inline value_type& front( void )
    {
        return vContent.front( );
    } // method

    inline value_type& back( void )
    {
        return vContent.back( );
    } // method

    inline const value_type& front( void ) const
    {
        return vContent.front( );
    } // method

    inline const value_type& back( void ) const
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

    inline TP_VEC::iterator insert( TP_VEC::const_iterator pos, const value_type& value )
    {
        return vContent.insert( pos, value );
    } // method

    template <class InputIt> inline TP_VEC::iterator insert( TP_VEC::const_iterator pos, InputIt first, InputIt last )
    {
        return vContent.insert( pos, first, last );
    } // method

    void reserve( size_type uiN )
    {
        vContent.reserve( uiN );
    } // method

    void resize( size_type uiN )
    {
        vContent.resize( uiN );
    } // method

    void clear( ) noexcept
    {
        vContent.clear( );
    } // method

    /// @brief returns unique seeds in this; shared seeds; unique seeds in pOther
    std::tuple<std::shared_ptr<Seeds>, std::shared_ptr<Seeds>, std::shared_ptr<Seeds>>
    splitSeedSets( std::shared_ptr<Seeds> pOther )
    {
        auto xRet =
            std::make_tuple( std::make_shared<Seeds>( ), std::make_shared<Seeds>( ), std::make_shared<Seeds>( ) );

        std::sort( vContent.begin( ), vContent.end( ), []( const Seed& rA, const Seed& rB ) {
            if( rA.start( ) != rB.start( ) )
                return rA.start( ) < rB.start( );
            if( rA.size( ) != rB.size( ) )
                return rA.size( ) < rB.size( );
            return rA.start_ref( ) < rB.start_ref( );
        } );

        std::sort( pOther->vContent.begin( ), pOther->vContent.end( ), []( const Seed& rA, const Seed& rB ) {
            if( rA.start( ) != rB.start( ) )
                return rA.start( ) < rB.start( );
            if( rA.size( ) != rB.size( ) )
                return rA.size( ) < rB.size( );
            return rA.start_ref( ) < rB.start_ref( );
        } );

        size_t uiI = 0;
        size_t uiJ = 0;
        while( uiI < vContent.size( ) && uiJ < pOther->size( ) )
        {
            if( vContent[ uiI ].start( ) == pOther->vContent[ uiJ ].start( ) &&
                vContent[ uiI ].size( ) == pOther->vContent[ uiJ ].size( ) &&
                vContent[ uiI ].start_ref( ) == pOther->vContent[ uiJ ].start_ref( ) )
            {
                std::get<1>( xRet )->push_back( vContent[ uiI ] );
                uiI++;
                uiJ++;
            } // if
            else if( vContent[ uiI ].start( ) < pOther->vContent[ uiJ ].start( ) ||
                     ( vContent[ uiI ].start( ) == pOther->vContent[ uiJ ].start( ) &&
                       ( vContent[ uiI ].size( ) < pOther->vContent[ uiJ ].size( ) ||
                         ( vContent[ uiI ].size( ) == pOther->vContent[ uiJ ].size( ) &&
                           vContent[ uiI ].start_ref( ) < pOther->vContent[ uiJ ].start_ref( ) ) ) ) )
            {
                // std::raise(SIGINT);
                std::get<0>( xRet )->push_back( vContent[ uiI ] );
                uiI++;
            } // if
            else
            {
                std::get<2>( xRet )->push_back( pOther->vContent[ uiJ ] );
                uiJ++;
            } // if
        } // while

        while( uiI < vContent.size( ) )
            std::get<0>( xRet )->push_back( vContent[ uiI++ ] );
        while( uiJ < pOther->size( ) )
            std::get<2>( xRet )->push_back( pOther->vContent[ uiJ++ ] );

        return xRet;
    } // method


    /// @brief returns unique seeds in this; shared seeds; unique seeds in pOther
    std::tuple<size_t, size_t, size_t> compareSeedSets( std::shared_ptr<Seeds> pOther )
    {
        auto xTemp = splitSeedSets( pOther );

        return std::make_tuple( std::get<0>( xTemp )->size( ), std::get<1>( xTemp )->size( ),
                                std::get<2>( xTemp )->size( ) );
    } // method

    void EXPORTED confirmSeedPositions( std::shared_ptr<NucSeq> pQuery, std::shared_ptr<Pack> pRef,
                                        bool bIsMaxExtended );
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