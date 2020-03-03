/* Authors: Arne Kutzner and Markus Schmidt
 * Created: Jan. 2020
 * MIT License
 * @file wkb_spatial.h
 * @brief Support of spatial types via WKB representation for SQL databases.
 */

#pragma once

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include "geom.h" // import of rectangle class
#include <algorithm>
#include <functional>
#include <iostream>
/// @endcond DOXYGEN_SHOW_SYSTEM_INCLUDES

class WKBPoint
{
  public:
    const static size_t uiSizeWKBHeader = 1 + // big/little endian a 1 byte
                                          4; // bytes for representation number of points in ring
    const static size_t uiSizeWKB = uiSizeWKBHeader + 2 * sizeof( double );

    std::array<uint8_t, uiSizeWKB> aData; // byte buffer

  protected:
    static size_t posOfPointX( size_t uiIdx )
    {
        assert( sizeof( double ) == 8 ); // assert that float has the correct size
        return uiSizeWKBHeader + uiIdx * 2 * sizeof( double ); // each point is expressed by two doubles (x and y);
    } // method

    static size_t posOfPointY( size_t uiIdx )
    {
        return posOfPointX( uiIdx ) + sizeof( double );
    } // method

    /** @brief Set WKBPoint header */
    void initHeader( )
    {
        set( 0, is_big_endian( ) ? 0x00 : 0x01 ); // big/little endian a 1 byte

        set( 1, is_big_endian( ) ? 0x00 : 0x01 ); // geom type a 4 bytes
        set( 2, 0 ); // geom type a 4 bytes
        set( 3, 0 ); // geom type a 4 bytes
        set( 4, is_big_endian( ) ? 0x01 : 0x00 ); // geom type a 4 bytes

    } // method

  public:
    /** @brief Constructor initializes header */
    WKBPoint( )
    {
        initHeader( );
    } // constructor

    WKBPoint( double fX, double fY ) : WKBPoint()
    {
        setPoint(0, fX, fY);
    } // constructor

    inline void set( size_t uiPos, uint8_t uiData )
    {
        aData[ uiPos ] = uiData;
    } // method

    inline void setDouble( size_t uiPos, double fData )
    {
        assert( uiSizeWKB >= uiPos + sizeof( double ) );
        memcpy( &aData[ uiPos ], &fData, sizeof( double ) );
    } // method

    inline size_t get( size_t uiI ) const
    {
        return aData[ uiI ];
    } // method

    inline double getDouble( size_t uiI ) const
    {
        assert( uiSizeWKB >= uiI + sizeof( double ) );

        double fRes;
        memcpy( &fRes, &aData[ uiI ], sizeof( double ) );
        return fRes;
    } // method

    inline void* getData( ) const
    {
        return (void*)&aData[ 0 ];
    } // method

    /** @brief Used by queries to set the WKBData */
    inline void setData( void* auiData )
    {
        std::copy_n( reinterpret_cast<uint8_t*>( auiData ), uiSizeWKB, aData.begin( ) );
    } // method

    template <typename IntegralType> inline void setPoint( size_t uiPointIdx, IntegralType valX, IntegralType valY )
    {
        setDouble( posOfPointX( uiPointIdx ), static_cast<double>( valX ) );
        setDouble( posOfPointY( uiPointIdx ), static_cast<double>( valY ) );
    } // method
}; // class

template <size_t NUM_POINTS> class WKBPolygon
{
  public:
    const static size_t uiSizeWKBHeader = 1 + // big/little endian a 1 byte
                                          4 + // 4 bytes for representation of geometry type
                                          4 + // bytes for representation if number of rings
                                          4; // bytes for representation number of points in ring
    const static size_t uiSizeWKB = uiSizeWKBHeader + NUM_POINTS * 2 * sizeof( double );

    std::array<uint8_t, uiSizeWKB> aData; // byte buffer

  protected:
    static size_t posOfPointX( size_t uiIdx )
    {
        assert( sizeof( double ) == 8 ); // assert that float has the correct size
        return uiSizeWKBHeader + uiIdx * 2 * sizeof( double ); // each point is expressed by two doubles (x and y);
    } // method

    static size_t posOfPointY( size_t uiIdx )
    {
        return posOfPointX( uiIdx ) + sizeof( double );
    } // method

    /** @brief Set WKBPolygon header */
    void initHeader( )
    {
        set( 0, is_big_endian( ) ? 0x00 : 0x01 ); // big/little endian a 1 byte

        set( 1, is_big_endian( ) ? 0x00 : 0x03 ); // geom type a 4 bytes
        set( 2, 0 ); // geom type a 4 bytes
        set( 3, 0 ); // geom type a 4 bytes
        set( 4, is_big_endian( ) ? 0x03 : 0x00 ); // geom type a 4 bytes

        set( 5, is_big_endian( ) ? 0x00 : 0x01 ); // number of rings a 4 bytes
        set( 6, 0 ); // number of rings a 4 bytes
        set( 7, 0 ); // number of rings a 4 bytes
        set( 8, is_big_endian( ) ? 0x01 : 0x00 ); // number of rings a 4 bytes

        set( 9, is_big_endian( ) ? 0x00 : NUM_POINTS ); // number of points a 4 bytes
        set( 10, 0 ); // number of points a 4 bytes
        set( 11, 0 ); // number of points a 4 bytes
        set( 12, is_big_endian( ) ? NUM_POINTS : 0x00 ); // number of points a 4 bytes
    } // method

  public:
    /** @brief Constructor initializes header */
    WKBPolygon( )
    {
        initHeader( );
    } // constructor

    inline void set( size_t uiPos, uint8_t uiData )
    {
        aData[ uiPos ] = uiData;
    } // method

    inline void setDouble( size_t uiPos, double fData )
    {
        assert( uiSizeWKB >= uiPos + sizeof( double ) );
        memcpy( &aData[ uiPos ], &fData, sizeof( double ) );
    } // method

    inline size_t get( size_t uiI ) const
    {
        return aData[ uiI ];
    } // method

    inline double getDouble( size_t uiI ) const
    {
        assert( uiSizeWKB >= uiI + sizeof( double ) );

        double fRes;
        memcpy( &fRes, &aData[ uiI ], sizeof( double ) );
        return fRes;
    } // method

    inline void* getData( ) const
    {
        return (void*)&aData[ 0 ];
    } // method

    /** @brief Used by queries to set the WKBData */
    inline void setData( void* auiData )
    {
        std::copy_n( reinterpret_cast<uint8_t*>( auiData ), uiSizeWKB, aData.begin( ) );
    } // method

    template <typename IntegralType> inline void setPoint( size_t uiPointIdx, IntegralType valX, IntegralType valY )
    {
        setDouble( posOfPointX( uiPointIdx ), static_cast<double>( valX ) );
        setDouble( posOfPointY( uiPointIdx ), static_cast<double>( valY ) );
    } // method
}; // class


template <typename T> class WKBRectangle : public WKBPolygon<5>
{
    /** @brief Uses the current rectangle data for computing the corresponding WKB data */
    void computeWKBDataFromRectangle( const geom::Rectangle<T>& rxRec )
    {
        // counterclockwise:
        // bottom left
        setPoint( 0, rxRec.xXAxis.start( ), rxRec.xYAxis.start( ) );
        // xData.setDouble(posOfPointX(0), static_cast<double>(this->xXAxis.start()));
        // xData.setDouble(posOfPointY(0), static_cast<double>(this->xYAxis.start()));

        // bottom right
        setPoint( 1, rxRec.xXAxis.end( ), rxRec.xYAxis.start( ) );
        // xData.setDouble(posOfPointX(1), static_cast<double>(this->xXAxis.end()));
        // xData.setDouble(posOfPointY(1), static_cast<double>(this->xYAxis.start()));

        // top right
        setPoint( 2, rxRec.xXAxis.end( ), rxRec.xYAxis.end( ) );
        // xData.setDouble(posOfPointX(2), static_cast<double>(this->xXAxis.end()));
        // xData.setDouble(posOfPointY(2), static_cast<double>(this->xYAxis.end()));

        // top left
        setPoint( 3, rxRec.xXAxis.start( ), rxRec.xYAxis.end( ) );
        // xData.setDouble(posOfPointX(3), static_cast<double>(this->xXAxis.start()));
        // xData.setDouble(posOfPointY(3), static_cast<double>(this->xYAxis.end()));

        // bottom left (again)
        setPoint( 4, rxRec.xXAxis.start( ), rxRec.xYAxis.start( ) );
        // xData.setDouble(posOfPointX(4), static_cast<double>(this->xXAxis.start()));
        // xData.setDouble(posOfPointY(4), static_cast<double>(this->xYAxis.start()));
    } // method

  public:
    WKBRectangle( const geom::Rectangle<T>& rxRect ) : WKBPolygon<5>( )
    {
        computeWKBDataFromRectangle( rxRect );
    } // constructor

    /** @brief (0, 0, 0, 0) rectangle. */
    WKBRectangle( ) : WKBRectangle<T>( geom::Rectangle<T>( 0, 0, 0, 0 ) )
    {} // default constructor

    geom::Rectangle<T> getRect( ) const
    {
        if( this->get( 0 ) != ( is_big_endian( ) ? 0x00 : 0x01 ) ) // check endian
            throw std::runtime_error( "WKB endian of DB does not match endian of system" );
        if( is_big_endian( ) && this->get( 4 ) != 0x03 )
            throw std::runtime_error( "WKB is no polygon" );
        else if( !is_big_endian( ) && this->get( 1 ) != 0x03 )
            throw std::runtime_error( "WKB is no polygon" );

        if( is_big_endian( ) && this->get( 8 ) != 0x01 )
            throw std::runtime_error( "WKB polygon has more than one ring" );
        else if( !is_big_endian( ) && this->get( 5 ) != 0x01 )
            throw std::runtime_error( "WKB polygon has more than one ring" );

        if( is_big_endian( ) && this->get( 12 ) != 0x05 )
            throw std::runtime_error( "WKB polygon does not have 4 points in ring" );
        else if( !is_big_endian( ) && this->get( 9 ) != 0x05 )
            throw std::runtime_error( "WKB polygon does not have 4 points in ring" );

        if( this->getDouble( posOfPointX( 0 ) ) != this->getDouble( posOfPointX( 4 ) ) )
            throw std::runtime_error(
                "WKB polygon is no (closed) rectangle (i.e. first points x does not match last points x)" );

        if( this->getDouble( posOfPointY( 0 ) ) != this->getDouble( posOfPointY( 4 ) ) )
            throw std::runtime_error(
                "WKB polygon is no (closed) rectangle (i.e. first points y does not match last points y)" );

        T uiXStart = (T)this->getDouble( posOfPointX( 0 ) );
        if( (T)this->getDouble( posOfPointX( 3 ) ) != uiXStart )
            throw std::runtime_error( "WKB polygon is no rectangle (i.e. angles are not rectangular uiXStart)" );

        T uiXEnd = (T)this->getDouble( posOfPointX( 1 ) );
        if( (T)this->getDouble( posOfPointX( 2 ) ) != uiXEnd )
            throw std::runtime_error( "WKB polygon is no rectangle (i.e. angles are not rectangular uiXEnd)" );
        if( uiXStart > uiXEnd )
            throw std::runtime_error( "WKB rectangle is in wrong order (i.e. uiXStart >= uiXEnd)" );

        T uiYStart = (T)this->getDouble( posOfPointY( 0 ) );
        if( (T)this->getDouble( posOfPointY( 1 ) ) != uiYStart )
            throw std::runtime_error( "WKB polygon is no rectangle (i.e. angles are not rectangular uiYStart)" );

        T uiYEnd = (T)this->getDouble( posOfPointY( 2 ) );
        if( (T)this->getDouble( posOfPointY( 3 ) ) != uiYEnd )
            throw std::runtime_error( "WKB polygon is no rectangle (i.e. angles are not rectangular uiYEnd)" );
        if( uiYStart > uiYEnd )
            throw std::runtime_error( "WKB rectangle is in wrong order (i.e. uiYStart >= uiYEnd)" );

        return geom::Rectangle<T>( uiXStart, uiYStart, uiXEnd - uiXStart, uiYEnd - uiYStart );
    } // method
}; // class

typedef WKBRectangle<uint64_t> WKBUint64Rectangle;

/// 
inline std::ostream& operator<<( std::ostream& xOS, const WKBUint64Rectangle& xRectWKB )
{
    const auto xRect = xRectWKB.getRect( );
	xOS << std::dec << "Rectangle: x: " << xRect.xXAxis << " y: " << xRect.xYAxis << std::endl;
    xOS << "WKBPolygon: ";
    for( auto uiI : xRectWKB.aData )
        xOS << std::hex << (int)uiI << " ";
    xOS << std::endl;
    return xOS;
} // operator