/**
 * @file nucSeq.h
 * @brief Implements NucSeq.
 * @author Arne Kutzner
 */
#pragma once

#include "container/seed.h"
#include "util/sqlite3.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <algorithm>
#include <array>
#include <cstring>
#include <memory>
#include <numeric>
/// @endcond

namespace libMA
{
// class GeneticSequence; //DEPRECATED
class NucSeq;

/** 32bit rounding to the next exponent as define
 */
#ifndef kroundup32
#define kroundup32( x )                                                                                                \
    ( --( x ),                                                                                                         \
      ( x ) |= ( x ) >> 1,                                                                                             \
      ( x ) |= ( x ) >> 2,                                                                                             \
      ( x ) |= ( x ) >> 4,                                                                                             \
      ( x ) |= ( x ) >> 8,                                                                                             \
      ( x ) |= ( x ) >> 16,                                                                                            \
      ++( x ) )
#endif

/** Generic reverse function, as it occurs in std::algorithms
 */
template <class T> void reverse( T word[], size_t length )
{
    char temp;
    for( size_t i = 0; i < length / 2; i++ )
    {
        temp = word[ i ];
        word[ i ] = word[ length - i - 1 ];
        word[ length - i - 1 ] = temp;
    } // for
} // reverse

#define WITH_QUALITY ( 1 )

/**
 * @brief Contains a genetic sequence made out of nucleotides (A, C, G, T).
 * @details
 * Class for genetic sequence that consist of nucleotides. (A, C, G, T)
 * @ingroup container
 */
class NucSeq : public Container
{
  public:
    /** The encapsulated sequence
     */
    uint8_t* pxSequenceRef = nullptr;
#if WITH_QUALITY
    uint8_t* pxQualityRef = nullptr;
#endif

    /** Current size of the content of the encapsulated sequence
     */
    size_t uiSize = 0;

    /** Current size of the buffer.
     */
    size_t uxCapacity = 0;


    /** Resets all protected attributes to its initial values.
     */
    inline void vReleaseMemory( )
    {
        /** Allocated memory will be released!
         */
        if( pxSequenceRef != NULL )
            free( pxSequenceRef );
#if WITH_QUALITY
        if( pxQualityRef != NULL )
            free( pxQualityRef );
#endif
    } // protected method

    void vResetProtectedAttributes( )
    {
        pxSequenceRef = NULL;
#if WITH_QUALITY
        pxQualityRef = NULL;
#endif
        uiSize = 0;
        uxCapacity = 0;
    } // protected method

#if WITH_QUALITY
    bool bHasQuality( ) const
    {
        return pxQualityRef != nullptr;
    } // method
#endif

    /** Tries to allocate the requested amount of memory and throws an exception if this process
     * fails. uxRequestedSize is expressed in "number of requested elements.
     */
    void vReserveMemory( size_t uxRequestedSize )
    {
        /* TO DO: This should be a bit more sophisticated ...
         */
        kroundup32( uxRequestedSize );

        /* We try to reserve the requested memory.
         * See:
         * http://stackoverflow.com/questions/1986538/how-to-handle-realloc-when-it-fails-due-to-memory
         */
        auto pxReallocRef = (uint8_t*)realloc( pxSequenceRef, uxRequestedSize * sizeof( uint8_t ) );
        uint8_t* pxReallocRef2 = NULL;
#if WITH_QUALITY
        if( bHasQuality( ) )
            pxReallocRef2 = (uint8_t*)realloc( pxQualityRef, uxRequestedSize * sizeof( uint8_t ) );
#endif


        if( pxReallocRef == NULL
#if WITH_QUALITY
            || ( bHasQuality( ) && pxReallocRef2 == NULL )
#endif
        )
        {
            throw AnnotatedException(
                ( std::string( "Memory Reallocation Failed for requested size " ) + std::to_string( uxRequestedSize ) )
                    .c_str( ) );
        } // if

        pxSequenceRef = pxReallocRef;
#if WITH_QUALITY
        if( bHasQuality( ) )
            pxQualityRef = pxReallocRef2;
#endif
        uxCapacity = uxRequestedSize;
    } // method

  public:
    /** The table used to translate from base pairs to numeric codes for nucleotides
     */
    static const EXPORTED unsigned char xNucleotideTranslationTable[ 256 ];

    std::string sName = "unknown";
    int64_t iId = -1;
    DEBUG( size_t uiFromLine = 0; ) // DEBUG

#if WITH_QUALITY
    void addQuality( )
    {
        assert( uiSize == 0 );
        size_t uxRequestedSize = 100;
        /* TO DO: This should be a bit more sophisticated ...
         */
        kroundup32( uxRequestedSize );

        /* We try to reserve the requested memory.
         * See:
         * http://stackoverflow.com/questions/1986538/how-to-handle-realloc-when-it-fails-due-to-memory
         */
        auto pxReallocRef = (uint8_t*)realloc( pxSequenceRef, uxRequestedSize * sizeof( uint8_t ) );
        auto pxReallocRef2 = (uint8_t*)realloc( pxQualityRef, uxRequestedSize * sizeof( uint8_t ) );


        if( pxReallocRef == NULL || pxReallocRef2 == NULL )
        {
            throw AnnotatedException(
                ( std::string( "Memory Reallocation Failed for requested size " ) + std::to_string( uxRequestedSize ) )
                    .c_str( ) );
        } // if

        pxSequenceRef = pxReallocRef;
        pxQualityRef = pxReallocRef2;
        uxCapacity = uxRequestedSize;
    } // methdo
#endif

    /** Default constructor
     */
    NucSeq( )
    {
        vResetProtectedAttributes( );
    } // default constructor

    /** Constructor that get the initial content of the sequence in text form.
     * FIX ME: This can be done a bit more efficient via the GeneticSequence class.
     */
    NucSeq( const std::string& rsInitialText )
    {
        vResetProtectedAttributes( );
        vAppend( rsInitialText.c_str( ) );
    } // constructor

#if WITH_QUALITY
    NucSeq( const std::string& rsInitialText, const uint8_t* pInitialQuality )
    {
        vResetProtectedAttributes( );
        addQuality( );
        vAppend( rsInitialText.c_str( ), pInitialQuality );
    } // constructor
#endif

    virtual ~NucSeq( )
    {
        /* Release all allocated memory.
         */
        vReleaseMemory( );
    } // destructor

    /** Move constructor on the foundation of text sequences.
     * Reuses the space of the text-sequence! TO DO: move & to &&
     */
    NucSeq( const NucSeq& rOther )
    {
        vResetProtectedAttributes( );
        vAppend( rOther.pxSequenceRef, rOther.uiSize );
    } // constructor


    /** is implicitly deleted by geneticSequence but boost python needs to know */
    // AK: VSC++ complained: NucSeq( const NucSeq& ) = delete;


    /** This moves the ownership of the protected attributes to another object.
     * The receiver of pxSequenceRef is responsible for its deletion.
     */
    void vTransferOwnership( NucSeq& rReceivingSequence )
    {
        /* We transport the three protected attributes to the receiver ...
         */
        rReceivingSequence.pxSequenceRef = this->pxSequenceRef;
#if WITH_QUALITY
        rReceivingSequence.pxQualityRef = this->pxQualityRef;
#endif
        rReceivingSequence.uiSize = this->uiSize;
        rReceivingSequence.uxCapacity = this->uxCapacity;

        /* ... and delete the knowledge here
         */
        vResetProtectedAttributes( );
    } // protected method

    /** Clears the inner sequence, but does not deallocate the memory.
     */
    inline void vClear( )
    {
        uiSize = 0;
    } // method

    /** Returns whether the sequence is empty or not.
     */
    inline bool bEmpty( )
    {
        return uiSize == 0;
    } // method

    inline bool empty( ) const
    {
        return uiSize == 0;
    } // method

    /** Fast getter and setter for element access.
     * If assertions activated we do a range check.
     */
    inline uint8_t operator[]( size_t uiSubscript ) const
    {
        assert( uiSubscript < uiSize );
        return pxSequenceRef[ uiSubscript ];
    } // method (get)
    inline uint8_t& operator[]( size_t uiSubscript )
    {
        assert( uiSubscript < uiSize );
        return pxSequenceRef[ uiSubscript ];
    } // method (set)

#if WITH_QUALITY
    inline uint8_t& quality( size_t uiSubscript )
    {
        assert( uiSubscript < uiSize );
        return pxQualityRef[ uiSubscript ];
    } // method

    inline uint8_t getQuality( size_t uiSubscript )
    {
        assert( uiSubscript < uiSize );
        return pxQualityRef[ uiSubscript ];
    } // method
#endif

    /** Resizes the internal buffer of the sequence to the requested value.
     */
    inline void resize( size_t uiRequestedSize ) // throws exception
    { /* Check, whether we have enough capacity, if not reserve memory
       */
        if( uxCapacity < uiRequestedSize )
        {
            vReserveMemory( uiRequestedSize );
        } // if

        uiSize = uiRequestedSize;
    } // method

    /** Because we want the reference to the sequence private we offer a getter method.
     * WARNING! Here you can get a null-pointer.
     */
    inline const uint8_t* const pGetSequenceRef( ) const
    {
        return this->pxSequenceRef;
    } // method

    /** Because we want to keep the size private we offer a getter method.
     */
    inline const size_t uxGetSequenceSize( ) const
    {
        return this->uiSize;
    } // method

    inline const size_t length( ) const
    {
        return this->uiSize;
    } // method

    /** Reverse the elements of the plain sequence.
     */
    inline void vReverse( )
    {
        reverse( pxSequenceRef, uiSize );

#if WITH_QUALITY
        if( bHasQuality( ) )
            reverse( pxQualityRef, uiSize );
#endif
    } // method

    inline void vReverse( size_t uiFrom, size_t uiTo )
    {
        reverse( pxSequenceRef + uiFrom, uiTo - uiFrom );
#if WITH_QUALITY
        if( bHasQuality( ) )
            reverse( pxQualityRef + uiFrom, uiTo - uiFrom );
#endif
    } // method

    inline void vReverseAll( )
    {
        vReverse( );
    } // method

#if WITH_QUALITY
    /** WARNING: the inner string might not null-terminated after this operation.
     */
    inline NucSeq& vAppend( const uint8_t* pSequence, const uint8_t* pQuality, size_t uxNumberOfElements )
    {
        assert( bHasQuality( ) );
        size_t uxRequestedSize = uxNumberOfElements + this->uiSize;

        if( uxCapacity < uxRequestedSize )
        {
            vReserveMemory( uxRequestedSize );
        } // if

        /** WARNING: If we work later with non 8-bit data we have to be careful here
         */
        memcpy( this->pxSequenceRef + uiSize, pSequence, uxNumberOfElements * sizeof( uint8_t ) );
        memcpy( this->pxQualityRef + uiSize, pQuality, uxNumberOfElements * sizeof( uint8_t ) );

        uiSize = uxRequestedSize;

        return *this;
    } // method
#endif

    /** WARNING: the inner string might not null-terminated after this operation.
     */
    inline NucSeq& vAppend( const uint8_t* pSequence, size_t uxNumberOfElements )
    {
        size_t uxRequestedSize = uxNumberOfElements + this->uiSize;

        if( uxCapacity < uxRequestedSize )
        {
            vReserveMemory( uxRequestedSize );
        } // if

        /** WARNING: If we work later with non 8-bit data we have to be careful here
         */
        memcpy( this->pxSequenceRef + uiSize, pSequence, uxNumberOfElements * sizeof( uint8_t ) );

        uiSize = uxRequestedSize;

        return *this;
    } // method

    /** WARNING: the inner string might not null-terminated after this operation.
     */
    inline NucSeq& vAppend( NucSeq& rOther )
    {
#if WITH_QUALITY
        assert( this->bHasQuality( ) == rOther.bHasQuality( ) );
        if( rOther.bHasQuality( ) )
            return this->vAppend( rOther.pxSequenceRef, rOther.pxQualityRef, rOther.length( ) );
        else
#endif
            return this->vAppend( rOther.pxSequenceRef, rOther.length( ) );
    } // method

    /** Push back of a single symbol.
     */
    inline void push_back( const uint8_t xElement )
    {
#if WITH_QUALITY
        assert( !bHasQuality( ) );
#endif
        if( this->uiSize >= this->uxCapacity )
        {
            vReserveMemory( this->uiSize + 1 );
        } // if

        pxSequenceRef[ uiSize ] = xElement;
        uiSize++;
    } // method

#if WITH_QUALITY
    /** Push back of a single symbol.
     */
    inline void push_back( const uint8_t xElement, const uint8_t xQuality )
    {
        assert( bHasQuality( ) );
        if( this->uiSize >= this->uxCapacity )
        {
            vReserveMemory( this->uiSize + 1 );
        } // if

        pxSequenceRef[ uiSize ] = xElement;
        pxQualityRef[ uiSize ] = xQuality;
        uiSize++;
    } // method
#endif

    /** Compares two sequences for equality
     */
    inline bool equal( const NucSeq& rOtherSequence )
    {
        if( this->uiSize == rOtherSequence.uiSize )
        {
            return memcmp( this->pxSequenceRef, rOtherSequence.pxSequenceRef, sizeof( uint8_t ) * uiSize ) == 0
#if WITH_QUALITY
                   && this->bHasQuality( ) == rOtherSequence.bHasQuality( ) &&
                   ( !this->bHasQuality( ) ||
                     memcmp( this->pxQualityRef, rOtherSequence.pxQualityRef, sizeof( uint8_t ) * uiSize ) == 0 )
#endif
                ;
        } // if
        return false;
    } // method

    // overload
    bool canCast( std::shared_ptr<Container> c ) const
    {
        return std::dynamic_pointer_cast<NucSeq>( c ) != nullptr;
    } // function

    // overload
    std::string getTypeName( ) const
    {
        return "NucSeq";
    } // function

    // overload
    std::shared_ptr<Container> getType( ) const
    {
        return std::shared_ptr<Container>( new NucSeq( ) );
    } // function

    /** Delivers the complement of a single nucleotide.
     */
    static inline char nucleotideComplement( char iNucleotide )
    {
        /* Complements of nucleotides
         *                               0  1  2  3
         */
        static const char chars[ 4 ] = {3, 2, 1, 0};

        return ( iNucleotide < 4 ) ? chars[ (int)iNucleotide ] : 5;
    } // static method

    /** Iterates over all base pairs in the sequence and creates the complement.
     * (A -> T, T -> A, C -> G, G -> C)
     */
    void vSwitchAllBasePairsToComplement( )
    {
        for( size_t uxIterator = 0; uxIterator < uiSize; uxIterator++ )
        {
            pxSequenceRef[ uxIterator ] = nucleotideComplement( pxSequenceRef[ uxIterator ] );
        } // for
    } // function

    /** transforms the character representation into a representation on the foundation of digits.
     */
    inline void vTranslateToNumericFormUsingTable( const unsigned char* alphabetTranslationTable, size_t uxStartIndex )
    {
        for( size_t uxIterator = uxStartIndex; uxIterator < uiSize; uxIterator++ )
        {
            pxSequenceRef[ uxIterator ] = alphabetTranslationTable[ pxSequenceRef[ uxIterator ] ];
        } // for
    } // method

    /** Gives the textual representation for some numeric representation.
     * Important: Keep this inline, so that it is not compiled into a function of its own.
     */
    static inline char translateACGTCodeToCharacter( uint8_t uiNucleotideCode )
    {
        static const char chars[ 4 ] = {'A', 'C', 'G', 'T'};
        if( uiNucleotideCode < 4 )
        {
            return chars[ uiNucleotideCode ];
        } // if
        else
        {
            return 'N';
        } // else
    } // static method

    /** transforms the numeric representation into a character representation.
     */
    inline void vTranslateToCharacterForm( size_t uxStartIndex )
    {
        for( size_t uxIterator = uxStartIndex; uxIterator < uiSize; uxIterator++ )
        {
            pxSequenceRef[ uxIterator ] = (uint8_t)translateACGTCodeToCharacter( pxSequenceRef[ uxIterator ] );
        } // for
    } // method

    /** transforms the character representation into a representation on the foundation of digits.
     */
    inline void vTranslateToNumericForm( size_t uxStartIndex )
    {
        vTranslateToNumericFormUsingTable( xNucleotideTranslationTable, uxStartIndex );
    } // method

    /** transforms the character representation into a representation on the foundation of digits.
     */
    inline void vTranslateToNumericForm( )
    {
        vTranslateToNumericForm( 0 );
    } // method

    /** transforms the numeric representation into a character representation.
     */
    inline void vTranslateToCharacterForm( )
    {
        vTranslateToCharacterForm( 0 );
    } // method

    /** The symbol on some position in textual form.
     * We count starting from 0.
     */
    inline char charAt( nucSeqIndex uxPosition )
    {
        if( uxPosition >= uiSize )
        {
            throw AnnotatedException( "Index out of range (charAt)" );
        } // if

        return translateACGTCodeToCharacter( pxSequenceRef[ uxPosition ] );
    } // method

    /** The symbol on some position in textual form.
     * We count starting from 0.
     */
    inline char compCharAt( size_t uxPosition )
    {
        if( uxPosition >= uiSize )
        {
            throw AnnotatedException( "Index out of range (compCharAt)" );
        } // if

        return translateACGTCodeToCharacter( nucleotideComplement( pxSequenceRef[ uxPosition ] ) );
    } // method

    /** Appends a string containing nucleotides as text and automatically translates the symbols.
     */
    void vAppend( const char* pcString )
    {
        size_t uxSizeBeforeAppendOperation = this->uiSize;
#if WITH_QUALITY
        assert( !bHasQuality( ) );
#endif

        /* WARNING! char and uint8_t must have the same size or we get a serious problem here!
         */
        vAppend( (const uint8_t*)pcString, strlen( pcString ) );

        vTranslateToNumericFormUsingTable( xNucleotideTranslationTable, uxSizeBeforeAppendOperation );
    } // method

#if WITH_QUALITY
    /** Appends a string containing nucleotides as text and automatically translates the symbols.
     */
    void vAppend( const char* pcString, const uint8_t* pQuality )
    {
        assert( bHasQuality( ) );
        size_t uxSizeBeforeAppendOperation = this->uiSize;

        /* WARNING! char and uint8_t must have the same size or we get a serious problem here!
         */
        vAppend( (const uint8_t*)pcString, pQuality, strlen( pcString ) );

        vTranslateToNumericFormUsingTable( xNucleotideTranslationTable, uxSizeBeforeAppendOperation );
    } // method
#endif

    /** wrapper for boost
     */
    void vAppend_boost( const char* pcString )
    {
        vAppend( pcString );
    } // method

    std::string toString( )
    {
        std::string ret = "";
        for( unsigned int i = 0; i < length( ); i++ )
            ret += charAt( i );
        return ret;
    } // function

    std::string fromToComplement( nucSeqIndex uiStart, nucSeqIndex uiEnd )
    {
        std::string ret = "";
        // for (unsigned int i = uiStart; i < uiEnd && i < length(); i++)
        for( nucSeqIndex i = uiEnd; i > uiStart; i-- )
            ret += compCharAt( i - 1 );
        return ret;
    } // function

    std::string fromTo( nucSeqIndex uiStart, nucSeqIndex uiEnd )
    {
        std::string ret = "";
        for( nucSeqIndex i = uiStart; i < uiEnd && i < length( ); i++ )
            ret += charAt( i );
        return ret;
    } // function

    std::string fromToQual( nucSeqIndex uiStart, nucSeqIndex uiEnd )
    {
#if WITH_QUALITY
        if( bHasQuality( ) )
        {
            std::string ret = "";
            for( nucSeqIndex i = uiStart; i < uiEnd && i < length( ); i++ )
                ret += (char)quality( i );
            return ret;
        } // if
#endif
        return "*";
    } // function

    std::string toQualString( )
    {
#if WITH_QUALITY
        if( bHasQuality( ) )
        {
            std::string ret = "";
            for( nucSeqIndex i = 0; i < length( ); i++ )
                ret += (char)quality( i );
            return ret;
        } // if
#endif
        return "*";
    } // function

    /* TO DO: Make the 5 a class constant!
     */
    inline uint8_t uxAlphabetSize( ) const
    {
        /* The Alphabet size for sequences of nucleotides is 5.
         */
        return 5; // eContentType == SEQUENCE_IS_NUCLEOTIDE ? 5 : 20;
    } // method

    std::string fastaq( )
    {
        std::string sRet = ">" + sName + "\n";
        for( unsigned int i = 0; i < length( ); i++ )
            sRet += charAt( i );
        sRet += "\n";
#if WITH_QUALITY
        if( bHasQuality( ) )
        {
            sRet += "+\n";
            for( unsigned int i = 0; i < length( ); i++ )
                sRet += (char)quality( i );
            sRet += "\n";
        } // if
#endif
        return sRet;
    } // method

    std::string fastaq_l( unsigned int uiLineLength )
    {
        std::string sRet = ">" + sName;
        for( unsigned int i = 0; i < length( ); i++ )
        {
            if( i % uiLineLength == 0 )
                sRet += "\n";
            sRet += charAt( i );
        } // for
        sRet += "\n";
#if WITH_QUALITY
        if( bHasQuality( ) )
        {
            sRet += "+";
            for( unsigned int i = 0; i < length( ); i++ )
            {
                if( i % uiLineLength == 0 )
                    sRet += "\n";
                sRet += (char)quality( i );
            } // for
            sRet += "\n";
        } // if
#endif
        return sRet;
    } // method

    /**
     * checks for untranslated characters in the sequence..
     */
    inline void check( )
    {
        for( unsigned int i = 0; i < length( ); i++ )
        {
            if( pxSequenceRef[ i ] > 4 )
            {
                // if was not allow print error and throw exception
                std::cerr << "Having invalid character in string: '" << pxSequenceRef[ i ] << "' at position: " << i
                          << " full fastaq: " << fastaq( ) << std::endl;
                throw AnnotatedException( "Found invalid character in nucSeq." );
            } // if
        } // for
    } // method

    /**
     * DEPRECATED
     * can not deal with 'N's
     */
    inline std::vector<uint8_t> as4Bit( nucSeqIndex uiFrom, nucSeqIndex uiTo, bool bReversed ) const
    {
        assert( uiTo <= length( ) );
        assert( uiFrom <= uiTo );
        DEBUG_3( for( size_t i = uiFrom; i < uiTo; i++ ) {
            assert( pxSequenceRef[ i ] < 4 );
            std::cout << (int)pxSequenceRef[ i ] << " ";
        } // for
                     std::cout
                     << std::endl; ) // DEBUG
        static const uint8_t aTranslate[ 4 ] = {1, 2, 4, 8};
        std::vector<uint8_t> vRet( uiTo - uiFrom - 1 );

        for( size_t i = 0; i < vRet.size( ); i++ )
            vRet[ bReversed ? vRet.size( ) - ( i + 1 ) : i ] = aTranslate[ pxSequenceRef[ i + uiFrom ] ];

        DEBUG_3( for( size_t i = 0; i < vRet.size( ); i++ ) std::cout << (int)vRet[ i ] << " ";
                 std::cout << std::endl; ) // DEBUG
        return vRet;
    } // method
}; // class NucSeq

class NucSeqSql : public SQL_BLOB
{
  public:
    std::shared_ptr<NucSeq> pNucSeq;

    NucSeqSql( std::shared_ptr<NucSeq> pNucSeq ) : pNucSeq( pNucSeq )
    {} // constructor

    NucSeqSql( )
    {} // default constructor

    const unsigned char* toBlob( ) const
    {
        if( pNucSeq == nullptr )
            return nullptr;
        return (unsigned char*)pNucSeq->pxSequenceRef;
    } // method

    const size_t blobSize( ) const
    {
        if( pNucSeq == nullptr )
            return 0;
        return pNucSeq->uiSize;
    } // method

    void fromBlob( const unsigned char* ucBlob, const size_t uiSize )
    {
        if( uiSize == 0 )
            return;
        pNucSeq = std::make_shared<NucSeq>( );
        pNucSeq->vClear( );
        pNucSeq->vAppend( (const uint8_t*)ucBlob, uiSize );
    } // method

}; // class

} // namespace libMA

#ifdef WITH_PYTHON
void exportNucSeq( py::module& rxPyModuleId );
#endif
