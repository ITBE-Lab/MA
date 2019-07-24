#pragma once

#include "util/support.h"
#include <vector>

namespace libMA
{
/**
 * @brief squeezed vector
 * @details
 * necessary for the line sweep clustering of the
 */
template <class TP_CONTENT> class SqueezedVector
{
  public:
    typedef std::vector<TP_CONTENT> TP_VEC;
    size_t uiNumElements;
    size_t uiSqueezeFactor;
    size_t uiCenterStripUp;
    size_t uiCenterStripDown;
    TP_VEC vContent; // the actual container that holds the content

    typedef typename TP_VEC::value_type value_type;
    typedef typename TP_VEC::size_type size_type;
    typedef typename TP_VEC::difference_type difference_type;
    typedef typename TP_VEC::iterator iterator;


    inline size_t physical_size_upper_squeezed_section( ) const
    {
        if( uiNumElements > uiCenterStripUp )
            return ( 1 + uiNumElements - uiCenterStripUp ) / uiSqueezeFactor;
        return 0;
    } // method

    inline size_t physical_size_lower_squeezed_section( ) const
    {
        if( uiNumElements > uiCenterStripDown )
            return ( 1 + uiNumElements - uiCenterStripDown ) / uiSqueezeFactor;
        return 0;
    } // method

    inline size_t physical_size( ) const
    {
        return physical_size_upper_squeezed_section( ) + physical_size_lower_squeezed_section( ) + uiCenterStripUp +
               uiCenterStripDown + 1;
    } // method

    inline size_t to_physical_coord( int64_t iX, int64_t iY ) const
    {
        // coordinate transformation
        iY -= iX;

        // check if we are above or below logical zero
        if( iY > (int64_t)uiCenterStripUp )
            iY = ( iY - uiCenterStripUp ) / uiSqueezeFactor + uiCenterStripUp;
        else if( iY < -((int64_t)uiCenterStripUp) )
            iY = ( iY + uiCenterStripDown ) / uiSqueezeFactor - uiCenterStripDown;
        // else
        //  no further transformation needed since we are in the non-squeezed part of the vector

        // move the logically negative part upwards so that the start of the coordinate system is zero
        iY += uiCenterStripDown + physical_size_lower_squeezed_section( );

        return (size_t)iY;
    } // method

    SqueezedVector( size_t uiNumElements, size_t uiSqueezeFactor, size_t uiCenterStripUp, size_t uiCenterStripDown )
        : uiNumElements( uiNumElements ),
          uiSqueezeFactor( uiSqueezeFactor ),
          uiCenterStripUp( uiCenterStripUp ),
          uiCenterStripDown( uiCenterStripDown ),
          vContent( this->physical_size( ) )
    {} // container vector

    // delete copy constructor
    SqueezedVector( const SqueezedVector& rOther ) = delete;

    TP_VEC& get( )
    {
        return vContent;
    } // function

    /**
     * @brief setter/getter; usage: vSqueezedVector[{iX, iY}];
     */
    inline TP_CONTENT& operator[]( std::tuple<int64_t, int64_t> tXY )
    {
        return vContent[ to_physical_coord( std::get<0>( tXY ), std::get<1>( tXY ) ) ];
    } // operator

    /**
     * @brief setter/getter; usage: vSqueezedVector[{iX, iY, iP}];
     * @details
     * iP is a physical offset
     */
    inline TP_CONTENT& operator[]( std::tuple<int64_t, int64_t, int64_t> tXYP )
    {
        return vContent[ to_physical_coord( std::get<0>( tXYP ), std::get<1>( tXYP ) ) + std::get<2>( tXYP ) ];
    } // operator

#if 0
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

    inline size_t size( void ) const
    {
        return vContent.size( );
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
#endif
}; // class

}; // namespace libMA