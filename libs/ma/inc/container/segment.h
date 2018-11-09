/**
 * @file segment.h
 * @brief Implements a the IntervalTree used for segmentation and various other related classes.
 * @author Markus Schmidt
 */
#ifndef INTERVALTREE_H
#define INTERVALTREE_H

#include "container/fMIndex.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <algorithm>
#include <numeric>
#include <thread>
/// @endcond

#define confMETA_MEASURE_DURATION ( 1 )

#define MEASURE_DURATIONS ( 1 )

namespace libMA
{
/**
 * @brief A Suffix Array Segment.
 * @details
 * A Suffix Array Segment is made up of two Intervals.
 * @li @c a SAInterval.
 * @li @c a Interval representing the position of the sequence on the query.
 * @ingroup container
 */
class Segment : public Container, public Interval<nucSeqIndex>
{
  public:
    SAInterval xSaInterval;
    /**
     * @brief Creates a new Segment.
     * @details Creates a new Segment on the base of a SAInterval and the
     * respective indices on the quey.
     */
    Segment( nucSeqIndex uiStart, nucSeqIndex uiSize, SAInterval xSaInterval )
        : Interval( uiStart, uiSize ), xSaInterval( xSaInterval )
    {} // constructor

    Segment( ) : Interval( ), xSaInterval( )
    {} // constructor

    Segment( const Segment& other ) : Interval( other ), xSaInterval( other.xSaInterval )
    {} // copy constructor


    // overload
    bool canCast( std::shared_ptr<Container> c ) const
    {
        return std::dynamic_pointer_cast<Segment>( c ) != nullptr;
    } // function

    // overload
    std::string getTypeName( ) const
    {
        return "Segment";
    } // function

    // overload
    std::shared_ptr<Container> getType( ) const
    {
        return std::shared_ptr<Container>( new Segment( ) );
    } // function

    /**
     * @brief The bwt interval within.
     * @returns the bwt interval within.
     */
    inline const SAInterval& saInterval( ) const
    {
        return xSaInterval;
    } // function

    /**
     * @brief Copys from another Segment.
     * @note override
     */
    inline Segment& operator=( const Segment& rxOther )
    {
        Interval::operator=( rxOther );
        xSaInterval = rxOther.xSaInterval;
        return *this;
    } // operator

}; // class ( Segment )


/**
 * @brief The segment tree.
 * @details
 * The segment "tree" is actually a doubly linked list.
 * The tree only exists logically,
 * meaning that the segments within the list represent the first layer of the tree initally.
 * Then after each iteration, the segments within the list represent the next layer down of the
 * tree.
 * @ingroup container
 */
class SegmentVector : public Container
{
  private:
    typedef std::vector<Segment> TP_VEC;
    ///@brief the content of the seed set
    TP_VEC vContent;

  public:
    typedef typename TP_VEC::value_type value_type;
    typedef typename TP_VEC::size_type size_type;
    typedef typename TP_VEC::difference_type difference_type;
    typedef typename TP_VEC::iterator iterator;
#if MEASURE_DURATIONS == ( 1 )
    // is default constructed
    double fExtraction = 0, fSorting = 0, fLinesweep = 0;
#endif

    SegmentVector( std::initializer_list<value_type> init ) : vContent( init )
    {} // initializer list constructor

    template <class InputIt> SegmentVector( InputIt xBegin, InputIt xEnd ) : vContent( xBegin, xEnd )
    {} // iterator constructor

    SegmentVector( ) : vContent( )
    {} // default constructor

    SegmentVector( const std::shared_ptr<SegmentVector> pOther ) : vContent( pOther->vContent )
    {} // constructor

    SegmentVector( size_t numElements ) : vContent( numElements )
    {} // constructor

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

    void clear( ) noexcept
    {
        vContent.clear( );
    } // method

    inline size_t numSeedsLarger( size_t uiMinSize ) const
    {
        return std::accumulate( this->begin( ), this->end( ),
                                (size_t)0, // initial value for the accumulation
                                [&uiMinSize]( size_t uiSum, const Segment& rSegment ) {
                                    // the seed is scored by how many times it is larger than the
                                    // uiMinSize.
                                    // @note rounded downwards by the cast.
                                    return uiSum + ( (size_t)rSegment.size( ) / uiMinSize );
                                } // lambda
        );
    } // method

    // overload
    bool canCast( std::shared_ptr<Container> c ) const
    {
        return std::dynamic_pointer_cast<SegmentVector>( c ) != nullptr;
    } // function

    // overload
    std::string getTypeName( ) const
    {
        return "SegmentVector";
    } // function

    // overload
    std::shared_ptr<Container> getType( ) const
    {
        return std::shared_ptr<SegmentVector>( new SegmentVector( ) );
    } // function

    /**
     * @brief Extracts all seeds from the tree.
     * @details
     * Calls fDo for all recorded hits.
     * fDo shall return false to terminate the iteration.
     * @Note pushBackBwtInterval records an interval of hits
     */
    template <class FUNCTOR>
    void forEachSeed( FMIndex& rxFMIndex, // std::shared_ptr<FMIndex> pxFMIndex,
                      nucSeqIndex uiQLen, size_t uiMAxAmbiguity, size_t uiMinLen, bool bSkip, FUNCTOR&& fDo
                      // std::function<bool(const Seed& s)> fDo
    )
    {
        // iterate over all the intervals that have been recorded using pushBackBwtInterval()
        for( const Segment& rSegment : *this )
        {
            if( rSegment.size( ) < uiMinLen )
                continue;
            // if the interval contains more than uiMAxAmbiguity hits it's of no importance and will
            // produce nothing but noise

            // if bSkip is not set uiJump by is used to not return more than uiMAxAmbiguity

            t_bwtIndex uiJumpBy = 1;
            if( rSegment.saInterval( ).size( ) > (t_bwtIndex)uiMAxAmbiguity && uiMAxAmbiguity != 0 )
            {
                if( bSkip )
                    continue;
                uiJumpBy = rSegment.saInterval( ).size( ) / uiMAxAmbiguity;
            } // if

            // iterate over the interval in the BWT
            for( auto ulCurrPos = rSegment.saInterval( ).start( ); ulCurrPos < rSegment.saInterval( ).end( );
                 ulCurrPos += uiJumpBy )
            {
                // calculate the referenceIndex using pxUsedFmIndex->bwt_sa() and call fDo for every
                // match individually
                nucSeqIndex ulIndexOnRefSeq = rxFMIndex.bwt_sa( ulCurrPos );
                nucSeqIndex uiPosOnQuery = rSegment.start( );
                bool bOnForw = ulIndexOnRefSeq < rxFMIndex.getRefSeqLength( ) / 2;
                if( !bOnForw )
                {
                    ulIndexOnRefSeq = rxFMIndex.getRefSeqLength( ) - ( ulIndexOnRefSeq + rSegment.size( ) + 1 );
                    assert( uiPosOnQuery < uiQLen );
                    uiPosOnQuery = uiQLen - (uiPosOnQuery + rSegment.size( ));
                } // if
                assert( ulIndexOnRefSeq < rxFMIndex.getRefSeqLength( ) / 2 );
                // call the given function
                if( !fDo( Seed( uiPosOnQuery, rSegment.size( ) + 1, ulIndexOnRefSeq,
                                (unsigned int)rSegment.saInterval( ).size( ), bOnForw ) ) )
                    return;
            } // for
        } // for
    } // function


    /**
     * @brief Extracts all seeds from the tree.
     * @details
     * Calls fDo for all recorded hits.
     * fDo shall return false to terminate the iteration.
     * @Note pushBackBwtInterval records an interval of hits
     */
    template <class FUNCTOR>
    void emplaceAllEachSeeds( FMIndex& rxFMIndex, nucSeqIndex uiQLen, size_t uiMAxAmbiguity, size_t uiMinLen,
                              Seeds& rvSeedVector,
                              FUNCTOR&& fDo // this function is called after each seed is emplaced
    )
    {
        forEachSeed( rxFMIndex, uiQLen, uiMAxAmbiguity, uiMinLen, true, [&]( Seed&& rS ) {
            rvSeedVector.push_back( rS );
            return fDo( );
        } );
    } // function

    /**
     * @brief Extracts all seeds from the segment list.
     */
    std::shared_ptr<Seeds> extractSeeds( std::shared_ptr<FMIndex> pxFMIndex, unsigned int uiMAxAmbiguity,
                                         unsigned int uiMinLen, nucSeqIndex uiQLen, bool bSkip = true )
    {
        std::shared_ptr<Seeds> pRet = std::shared_ptr<Seeds>( new Seeds( ) );
        forEachSeed( *pxFMIndex, uiQLen, uiMAxAmbiguity, uiMinLen, bSkip,
                     [&pRet]( const Seed& s ) {
                         pRet->push_back( s );
                         return true;
                     } // lambda
        ); // for each
        return pRet;
    } // function


    /**
     * @brief returns the number of seeds
     */
    inline uint64_t numSeeds( unsigned int max_size ) const
    {
        uint64_t uiTotal = 0;
        for( const Segment& rSegment : *this )
            if( max_size == 0 || rSegment.xSaInterval.size( ) <= max_size )
                uiTotal += rSegment.xSaInterval.size( );
        return uiTotal;
    } // function
}; // class
} // namespace libMA

#ifdef WITH_PYTHON
/**
 * @brief Exposes the SegmentVector to boost python.
 * @ingroup export
 */
#ifdef WITH_BOOST
void exportSegment( );
#else
void exportSegment( py::module& rxPyModuleId );
#endif
#endif


#endif