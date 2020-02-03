/**
 * @file geom.h
 * @brief Implements a generic Interval class.
 * @author Markus Schmidt
 */
#pragma once

#include "exception.h"
#include "support.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <algorithm>
#include <cstring>
#include <functional>
#include <iostream>
/// @endcond DOXYGEN_SHOW_SYSTEM_INCLUDES

namespace geomUtil
{
/**
 * @brief A generic multipurpose Interval.
 */
template <typename T> class Interval
{
  public:
    /// @brief Start position of interval.
    T iStart = 0;
    /// @brief Size of interval.
    T iSize = 0;
    /**
     * @brief Creates a new Interval.
     */
    Interval( T iStart, T iSize ) : iStart( iStart ), iSize( iSize )
    {} // constructor

    /**
     * @brief Copys from another Interval.
     */
    Interval( const Interval& c ) : iStart( c.iStart ), iSize( c.iSize )
    {} // copy constructor

    /**
     * @brief Default empty constructor.
     */
    Interval( )
    {} // default constructor

    inline static Interval start_end( T start, T end )
    {
        Interval xRet( start, 0 );
        xRet.end( end );
        return xRet;
    } // function

    /**
     * @returns the end of the interval.
     * @brief Interval end.
     * @note end = start + size
     */
    inline const T end( ) const
    {
        return iStart + iSize;
    } // function

    ///@brief Wrapper for boost-python.
    inline const T end_boost1( ) const
    {
        return iStart + iSize;
    } // function

    /**
     * @brief Allows chaning the end of the interval.
     * @note end = start + size
     */
    inline void end( const T iVal )
    {
        iSize = iVal - iStart;
    } // function

    ///@brief Wrapper for boost-python.
    inline void end_boost2( const T iVal )
    {
        end( iVal );
    } // function

    /**
     * @returns the start of the interval.
     * @brief Interval start.
     */
    inline const T start( ) const
    {
        return iStart;
    } // function

    ///@brief Wrapper for boost-python.
    inline const T start_boost1( ) const
    {
        return iStart;
    } // function

    /**
     * @brief allows chaning the beginning of the interval.
     */
    inline void start( const T iVal )
    {
        T iEnd = end( );
        iStart = iVal;
        end( iEnd );
    } // function

    ///@brief Wrapper for boost-python.
    inline void start_boost2( const T iVal )
    {
        start( iVal );
    } // function

    /**
     * @brief allows chaning the size of the interval.
     */
    inline void size( const T iVal )
    {
        iSize = iVal;
    } // function

    ///@brief wrapper for boost-python
    inline void size_boost2( const T iVal )
    {
        iSize = iVal;
    } // function

    /**
     * @returns the center of the interval.
     * @brief The center of the interval.
     */
    inline const T center( ) const
    {
        return start( ) + size( ) / 2;
    } // function

    /**
     * @returns the size of the interval.
     * @brief The size of the interval.
     */
    inline const T size( ) const
    {
        return iSize;
    } // function

    ///@brief wrapper for boost-python
    inline const T size_boost1( ) const
    {
        return iSize;
    } // function

    /**
     * @brief allows chaning the start and size of the interval.
     */
    inline void set( T iStart, T iSize )
    {
        start( iStart );
        size( iSize );
    } // function

    /*
     * @returns the interval start and end for i = 0 and 1 respectively.
     * @brief Interval start for 0 and end for 1.
     * @details
     * Throws a nullpointer exception for any other i.
     */
    inline const T operator[]( const std::size_t i ) const
    {
        if( i == 0 )
            return start( );
        if( i == 1 )
            return end( );
        throw AnnotatedException( "can only access index 0 and 1 of interval" );
    } // operator

    /*
     * @brief copys from another Interval.
     */
    inline Interval& operator=( const Interval& rxOther )
    {
        iStart = rxOther.iStart;
        iSize = rxOther.iSize;
        return *this;
    } // operator

    /*
     * @brief compares two Intervals.
     * @returns true if start and size are equal, false otherwise.
     */
    inline bool operator==( const Interval& rxOther )
    {
        return iStart == rxOther.iStart && iSize == rxOther.iSize;
    } // operator
}; // class

/**
 * @brief well known binary (WKB) of a polygon
 * @details
 * Object that holds a WKB of a polygon.
 * Can be encoded or decoded from a rectangle.
 * Main purpose of this is to hold the WKB data while an insert is performed into a database.
 */
template <size_t SIZE> class WKBPolygon
{
  public:
    std::array<uint8_t, SIZE> aData;

    WKBPolygon( )
    {} // constructor

    WKBPolygon( const WKBPolygon& rOther ) : aData( rOther.aData )
    {} // copy constructor

    inline void set( size_t uiPos, uint8_t uiData )
    {
        aData[ uiPos ] = uiData;
    } // method

    inline void setDouble( size_t uiPos, double fData )
    {
        assert( SIZE >= uiPos + sizeof( double ) );
        std::memcpy( &aData[ uiPos ], &fData, sizeof( double ) );
    } // method

    inline size_t get( size_t uiI ) const
    {
        return aData[ uiI ];
    } // method

    inline double getDouble( size_t uiI ) const
    {
        assert( SIZE >= uiI + sizeof( double ) );

        double fRes;
        std::memcpy( &fRes, &aData[ uiI ], sizeof( double ) );
        return fRes;
    } // method

    inline void* getData( ) const
    {
        return (void*)&aData[ 0 ];
    } // method

    inline void setData( void* auiData )
    {
        std::copy_n( reinterpret_cast<uint8_t*>( auiData ), SIZE, aData.begin( ) );
    } // method
}; // class

/**
 * @brief a rectangle represented by two intervals.
 */
template <typename T> class Rectangle
{
  public:
    /// @brief number of bytes in a WKB header of a polygon with 1 ring (add one extra byte for each extra ring)
    const static size_t uiSizeWKBHeader = 1 + // big/little endian a 1 byte
                                          4 + // geom type a 4 bytes
                                          4 + // number of rings
                                          4; // number of points in ring

    /// @brief number of bytes in a WKB representation of a rectangle
    /// @details header + 5 points (first and last point are equal) a two coordinates of type double.
    const static size_t uiSizeWKB = uiSizeWKBHeader + 5 * 2 * sizeof( double );

  private:
    /**
     * @brief index of the X coordinate of point with index uiIdx in a WKB representation of a rectangle
     * @details
     * the coordinate of each point is several bytes long, this returns the index of the first byte
     */
    static size_t posOfPointX( size_t uiIdx )
    {
        assert( sizeof( double ) == 8 ); // assert that float has the correct size
        return uiSizeWKBHeader + uiIdx * 2 * sizeof( double ); // each point is expressed by two doubles (x and y);
    } // method

    /**
     * @brief index of the Y coordinate of point with index uiIdx in a WKB representation of a rectangle
     * @details
     * the coordinate of each point is several bytes long, this returns the index of the first byte
     */
    static size_t posOfPointY( size_t uiIdx )
    {
        return posOfPointX( uiIdx ) + sizeof( double );
    } // method

  public:
    Interval<T> xXAxis, xYAxis;
    /**
     * @brief Creates a new Rectangle.
     */
    Rectangle( T iStartX, T iStartY, T iSizeX, T iSizeY ) : xXAxis( iStartX, iSizeX ), xYAxis( iStartY, iSizeY )
    {} // constructor

    Rectangle( ) : Rectangle( 0, 0, 0, 0 )
    {} // default constructor

    /**
     * @brief encode rectangle as WKB polygon
     * @details
     * WKB representation is a follows:
     * - data is given as big (0) / little (1) endian          a 1 byte
     * - geometry type (3 == polygon)                          a 4 byte
     * - number of rings in the polygon (1 for a rectangle)    a 4 byte
     * - number of points in the polygon (5 for a rectangle)   a 80 byte (=8*5*2)
     * The first point must be repeated at the end to show that the polygon is "closed".
     * Each point has two coordinates (x and y) a 8 byte.
     */
    inline WKBPolygon<uiSizeWKB> getWKB( )
    {
        WKBPolygon<uiSizeWKB> xData;
        xData.set( 0, is_big_endian( ) ? 0x00 : 0x01 ); // big/little endian a 1 byte

        xData.set( 1, is_big_endian( ) ? 0x00 : 0x03 ); // geom type a 4 bytes
        xData.set( 2, 0 ); // geom type a 4 bytes
        xData.set( 3, 0 ); // geom type a 4 bytes
        xData.set( 4, is_big_endian( ) ? 0x03 : 0x00 ); // geom type a 4 bytes

        xData.set( 5, is_big_endian( ) ? 0x00 : 0x01 ); // number of rings a 4 bytes
        xData.set( 6, 0 ); // number of rings a 4 bytes
        xData.set( 7, 0 ); // number of rings a 4 bytes
        xData.set( 8, is_big_endian( ) ? 0x01 : 0x00 ); // number of rings a 4 bytes

        xData.set( 9, is_big_endian( ) ? 0x00 : 0x05 ); // number of points a 4 bytes
        xData.set( 10, 0 ); // number of points a 4 bytes
        xData.set( 11, 0 ); // number of points a 4 bytes
        xData.set( 12, is_big_endian( ) ? 0x05 : 0x00 ); // number of points a 4 bytes

        // counterclockwise:
        // bottom left
        xData.setDouble( posOfPointX( 0 ), (double)xXAxis.start( ) );
        xData.setDouble( posOfPointY( 0 ), (double)xYAxis.start( ) );

        // bottom right
        xData.setDouble( posOfPointX( 1 ), (double)xXAxis.end( ) );
        xData.setDouble( posOfPointY( 1 ), (double)xYAxis.start( ) );

        // top right
        xData.setDouble( posOfPointX( 2 ), (double)xXAxis.end( ) );
        xData.setDouble( posOfPointY( 2 ), (double)xYAxis.end( ) );

        // top left
        xData.setDouble( posOfPointX( 3 ), (double)xXAxis.start( ) );
        xData.setDouble( posOfPointY( 3 ), (double)xYAxis.end( ) );

        // bottom left (again)
        xData.setDouble( posOfPointX( 4 ), (double)xXAxis.start( ) );
        xData.setDouble( posOfPointY( 4 ), (double)xYAxis.start( ) );

        return xData;
    } // method

    /**
     * @brief decode rectangle from WKB polygon
     * @details
     * See getWKB.
     * Throws exception if WKB polygon is not a rectangle.
     */
    inline void fromWKB( WKBPolygon<uiSizeWKB>& xData )
    {
        // @todo might cause trouble if data from systems with different endian is inserted into the DB
        // in this case code needs to be written to change the endian of xData.
        if( xData.get( 0 ) != ( is_big_endian( ) ? 0x00 : 0x01 ) ) // check endian
            throw std::runtime_error( "WKB endian of DB does not match endian of system" );
        if( is_big_endian( ) && xData.get( 4 ) != 0x03 )
            throw std::runtime_error( "WKB is no polygon" );
        else if( !is_big_endian( ) && xData.get( 1 ) != 0x03 )
            throw std::runtime_error( "WKB is no polygon" );

        if( is_big_endian( ) && xData.get( 8 ) != 0x01 )
            throw std::runtime_error( "WKB polygon has more than one ring" );
        else if( !is_big_endian( ) && xData.get( 5 ) != 0x01 )
            throw std::runtime_error( "WKB polygon has more than one ring" );

        if( is_big_endian( ) && xData.get( 12 ) != 0x05 )
            throw std::runtime_error( "WKB polygon does not have 4 points in ring" );
        else if( !is_big_endian( ) && xData.get( 9 ) != 0x05 )
            throw std::runtime_error( "WKB polygon does not have 4 points in ring" );

        if( xData.getDouble( posOfPointX( 0 ) ) != xData.getDouble( posOfPointX( 4 ) ) )
            throw std::runtime_error(
                "WKB polygon is no (closed) rectangle (i.e. first points x does not match last points x)" );

        if( xData.getDouble( posOfPointY( 0 ) ) != xData.getDouble( posOfPointY( 4 ) ) )
            throw std::runtime_error(
                "WKB polygon is no (closed) rectangle (i.e. first points y does not match last points y)" );

        T uiXStart = (T)xData.getDouble( posOfPointX( 0 ) );
        if( (T)xData.getDouble( posOfPointX( 3 ) ) != uiXStart )
            throw std::runtime_error( "WKB polygon is no rectangle (i.e. angles are not rectangular uiXStart)" );

        T uiXEnd = (T)xData.getDouble( posOfPointX( 1 ) );
        if( (T)xData.getDouble( posOfPointX( 2 ) ) != uiXEnd )
            throw std::runtime_error( "WKB polygon is no rectangle (i.e. angles are not rectangular uiXEnd)" );
        if( uiXStart > uiXEnd )
            throw std::runtime_error( "WKB rectangle is in wrong order (i.e. uiXStart >= uiXEnd)" );

        T uiYStart = (T)xData.getDouble( posOfPointY( 0 ) );
        if( (T)xData.getDouble( posOfPointY( 1 ) ) != uiYStart )
            throw std::runtime_error( "WKB polygon is no rectangle (i.e. angles are not rectangular uiYStart)" );

        T uiYEnd = (T)xData.getDouble( posOfPointY( 2 ) );
        if( (T)xData.getDouble( posOfPointY( 3 ) ) != uiYEnd )
            throw std::runtime_error( "WKB polygon is no rectangle (i.e. angles are not rectangular uiYEnd)" );
        if( uiYStart > uiYEnd )
            throw std::runtime_error( "WKB rectangle is in wrong order (i.e. uiYStart >= uiYEnd)" );

        xXAxis.iStart = uiXStart;
        xYAxis.iStart = uiYStart;
        xXAxis.iSize = uiXEnd - uiXStart;
        xYAxis.iSize = uiYEnd - uiYStart;
    } // method

    /*
     * @brief compares two Rectangles.
     * @returns true if origin and size are equal, false otherwise.
     */
    inline bool operator==( const Rectangle<T>& rxOther )
    {
        return xXAxis == rxOther.xXAxis && xYAxis == rxOther.xYAxis;
    } // operator
}; // class

} // namespace geomUtil