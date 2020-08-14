/**
 * @file container.h
 * @brief Implements Container.
 * @author Markus Schmidt
 */

#pragma once

#include "ms/util/export.h"
#include "util/debug.h"
#include "util/support.h"
#include <stdexcept>

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <memory>

/// @endcond

#if DEBUG_LEVEL >= 1
#define TOMBSTONE_VAL_ALIVE 1010101010101
#define TOMBSTONE_VAL_DEAD 2020202020202
#endif

namespace libMS
{

/**
 * @defgroup container Containers
 * @brief All classes containing data elements.
 * @details
 * All classes that store data should inherit from Container.
 * These classes are then added to this Group.
 */

/**
 * @brief Abstract class intended to hold data objects used by @ref Module "modules".
 * @details
 * All classes containing data should inherit from this class.
 * @ingroup container
 */
class Container
{
  public:
#if DEBUG_LEVEL >= 1
    size_t uiTombStone = TOMBSTONE_VAL_ALIVE;

    virtual ~Container( )
    {
        uiTombStone = TOMBSTONE_VAL_DEAD;
    } // deconstructor
#else
    /* Clang requires a virtual destructor for classes comprising virtual methods */
    virtual ~Container( )
    {} // deconstructor
#endif

    virtual void dummy( )
    {} // make class abstract
}; // class

/**
 * @brief A vector of containers, that itself is a Container
 * @details
 * This can be used to define a vector as input and output of Modules.
 * @ingroup container
 */
template <class TP_CONTENT> class ContainerVector : public Container
{
  public:
    typedef std::vector<TP_CONTENT> TP_VEC;
    TP_VEC vContent;

    typedef typename TP_VEC::value_type value_type;
    typedef typename TP_VEC::size_type size_type;
    typedef typename TP_VEC::difference_type difference_type;
    typedef typename TP_VEC::iterator iterator;


    ContainerVector( std::initializer_list<TP_CONTENT> init ) : vContent( init )
    {} // initializer list constructor

    template <class InputIt> ContainerVector( InputIt xBegin, InputIt xEnd ) : vContent( xBegin, xEnd )
    {} // iterator constructor

    ContainerVector( ) : vContent( )
    {} // container vector

    ContainerVector( TP_CONTENT contentType ) : vContent( )
    {} // container vector

    ContainerVector( const std::shared_ptr<ContainerVector<TP_CONTENT>> pOther ) : vContent( pOther->vContent )
    {} // container vector

    ContainerVector( std::shared_ptr<std::vector<TP_CONTENT>> pContent ) : vContent( *pContent )
    {} // container vector

    ContainerVector( size_t numElements ) : vContent( numElements )
    {} // container vector

    // delete copy constructor
    ContainerVector( const ContainerVector& rOther ) = delete;

    std::vector<TP_CONTENT> get( )
    {
        return vContent;
    } // function

    // setter
    inline TP_CONTENT& operator[]( size_t uiI )
    {
        return vContent[ uiI ];
    } // operator

    // getter
    inline const TP_CONTENT& operator[]( size_t uiI ) const
    {
        return vContent[ uiI ];
    } // operator

    inline void push_back( const TP_CONTENT& value )
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

    inline size_t size( void ) const
    {
        return vContent.size( );
    } // method

    inline bool empty( void ) const
    {
        return vContent.empty( );
    } // method

    inline TP_CONTENT& front( void )
    {
        return vContent.front( );
    } // method

    inline TP_CONTENT& back( void )
    {
        return vContent.back( );
    } // method

    inline const TP_CONTENT& front( void ) const
    {
        return vContent.front( );
    } // method

    inline const TP_CONTENT& back( void ) const
    {
        return vContent.back( );
    } // method

    inline typename TP_VEC::iterator begin( void ) noexcept
    {
        return vContent.begin( );
    } // method

    inline typename TP_VEC::iterator end( void ) noexcept
    {
        return vContent.end( );
    } // method

    inline typename TP_VEC::const_iterator begin( void ) const noexcept
    {
        return vContent.begin( );
    } // method

    inline typename TP_VEC::const_iterator end( void ) const noexcept
    {
        return vContent.end( );
    } // method

    inline void erase( typename TP_VEC::iterator pos )
    {
        vContent.erase( pos );
    } // method

    inline void erase( typename TP_VEC::iterator first, typename TP_VEC::iterator last )
    {
        vContent.erase( first, last );
    } // method

    inline typename TP_VEC::iterator insert( typename TP_VEC::const_iterator pos, const TP_CONTENT& value )
    {
        return vContent.insert( pos, value );
    } // method

    template <class InputIt>
    inline typename TP_VEC::iterator insert( typename TP_VEC::const_iterator pos, InputIt first, InputIt last )
    {
        return vContent.insert( pos, first, last );
    } // method

    inline void reserve( size_type n )
    {
        vContent.reserve( n );
    } // method

}; // class

template <class TP_CONTENT>
inline bool operator==( const ContainerVector<TP_CONTENT>& rLeft, const ContainerVector<TP_CONTENT>& rRight )
{
    return rLeft.vContent == rRight.vContent;
} // function

template <class TP_CONTENT>
inline bool operator!=( const ContainerVector<TP_CONTENT>& rLeft, const ContainerVector<TP_CONTENT>& rRight )
{
    return rLeft.vContent != rRight.vContent;
} // function

template <class TP_CONTENT>
inline bool operator<( const ContainerVector<TP_CONTENT>& rLeft, const ContainerVector<TP_CONTENT>& rRight )
{
    return rLeft.vContent < rRight.vContent;
} // function

template <class TP_CONTENT>
inline bool operator<=( const ContainerVector<TP_CONTENT>& rLeft, const ContainerVector<TP_CONTENT>& rRight )
{
    return rLeft.vContent <= rRight.vContent;
} // function

template <class TP_CONTENT>
inline bool operator>( const ContainerVector<TP_CONTENT>& rLeft, const ContainerVector<TP_CONTENT>& rRight )
{
    return rLeft.vContent > rRight.vContent;
} // function

template <class TP_CONTENT>
inline bool operator>=( const ContainerVector<TP_CONTENT>& rLeft, const ContainerVector<TP_CONTENT>& rRight )
{
    return rLeft.vContent >= rRight.vContent;
} // function

/**
 * @brief input for all python modules.
 * @details
 * Vector of containers.
 */
class PyContainerVector : public ContainerVector<std::shared_ptr<Container>>
{
    // use the constructors of ContainerVector<std::shared_ptr<Container>>
    using ContainerVector<std::shared_ptr<Container>>::ContainerVector;
}; // class

} // namespace libMS


#ifdef WITH_PYTHON
/**
 * @brief Function to export Container to boost python.
 * @ingroup export
 */

void exportContainer( libMS::SubmoduleOrganizer& xOrganizer );
#endif
