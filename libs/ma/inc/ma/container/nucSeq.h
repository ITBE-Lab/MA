/**
 * @file nucSeq.h
 * @brief Implements NucSeq.
 * @author Arne Kutzner
 */
#pragma once

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <algorithm>
#include <array>
#include <cstring>
#include <memory>
#include <numeric>
/// @endcond

#include "ma/container/seed.h"
#include <db_base.h> // import database support

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
class NucSeq : public libMS::Container
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
            throw std::runtime_error(
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
    static const DLL_PORT( MA ) unsigned char xNucleotideTranslationTable[ 256 ];

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
            throw std::runtime_error(
                ( std::string( "Memory Reallocation Failed for requested size " ) + std::to_string( uxRequestedSize ) )
                    .c_str( ) );
        } // if

        pxSequenceRef = pxReallocRef;
        pxQualityRef = pxReallocRef2;
        uxCapacity = uxRequestedSize;
    } // method
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

    /** @brief Move constructor on the foundation of text sequences.
     * @details
     * Reuses the space of the text-sequence! TO DO: move & to &&
     */
    NucSeq( const NucSeq& rOther )
    {
        vResetProtectedAttributes( );
        vAppend( rOther.pxSequenceRef, rOther.uiSize );
    } // constructor

    /**  @brief creates a copy of the subsequence [uiFrom, uiTo)
     */
    NucSeq( const NucSeq& rOther, nucSeqIndex uiFrom, nucSeqIndex uiTo )
    {
        vResetProtectedAttributes( );
        vAppend( rOther.pxSequenceRef + uiFrom, uiTo - uiFrom );
    } // constructor

    /** @brief creates a copy of the subsequence [uiFrom, uiTo)
     */
    NucSeq( std::shared_ptr<NucSeq> pOther, nucSeqIndex uiFrom, nucSeqIndex uiTo ) : NucSeq( *pOther, uiFrom, uiTo )
    {} // constructor


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

    inline bool isEqual( NucSeq& xOther )
    {
        if( length( ) != xOther.length( ) )
            return false;
        for( size_t uiI = 0; uiI < length( ); uiI++ )
            if( ( *this )[ uiI ] != xOther[ uiI ] )
                return false;
#if WITH_QUALITY
        for( size_t uiI = 0; uiI < length( ); uiI++ )
            if( this->getQuality( uiI ) != xOther.getQuality( uiI ) )
                return false;
#endif
        return true;
    }

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
    bool canCast( std::shared_ptr<libMS::Container> c ) const
    {
        return std::dynamic_pointer_cast<NucSeq>( c ) != nullptr;
    } // function

    // overload
    std::string getTypeName( ) const
    {
        return "NucSeq";
    } // function

    // overload
    std::shared_ptr<libMS::Container> getType( ) const
    {
        return std::shared_ptr<libMS::Container>( new NucSeq( ) );
    } // function

    /** Delivers the complement of a single nucleotide.
     */
    static inline char nucleotideComplement( char iNucleotide )
    {
        /* Complements of nucleotides
         *                               0  1  2  3
         */
        static const char chars[ 4 ] = { 3, 2, 1, 0 };

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
        static const char chars[ 4 ] = { 'A', 'C', 'G', 'T' };
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
    inline char charAt( nucSeqIndex uxPosition ) const
    {
        if( uxPosition >= uiSize )
        {
            throw std::runtime_error( "Index out of range (charAt)" );
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
            throw std::runtime_error( "Index out of range (compCharAt)" );
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

    std::string toString( ) const
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

    std::string toStringComplement( )
    {
        return fromToComplement( 0, length( ) );
    } // function

    std::string fromTo( nucSeqIndex uiStart, nucSeqIndex uiEnd ) const
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
                throw std::runtime_error( "Found invalid character in nucSeq." );
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
        static const uint8_t aTranslate[ 4 ] = { 1, 2, 4, 8 };
        std::vector<uint8_t> vRet( uiTo - uiFrom - 1 );

        for( size_t i = 0; i < vRet.size( ); i++ )
            vRet[ bReversed ? vRet.size( ) - ( i + 1 ) : i ] = aTranslate[ pxSequenceRef[ i + uiFrom ] ];

        DEBUG_3( for( size_t i = 0; i < vRet.size( ); i++ ) std::cout << (int)vRet[ i ] << " ";
                 std::cout << std::endl; ) // DEBUG
        return vRet;
    } // method
}; // class NucSeq

/** @brief Compressed representation of nucleotide sequences without quality information.
 *  @details
 *  Header: uint64 value that tells the size of the uncompressed sequence.
 *  Compression scheme:
 * - bit 7 = 1               : 2 bit representation
 * -             bit 6 = 0   : bits 0 - 5 encode three A,C,T,G and next byte is not special (first symbol in bits 5,4)
 * -             bit 6 = 1   : like above, but additionally, the following byte encodes 4 A,C,G,T (fourth symbol in
 * -                           bits 7,6)
 * - bit 7 = 0   bit 6 = 1   : N symbols, the number of N symbols is encoded bits 0 - 5. PLEASE NOTE:
 *                             Because there is at least one N we store the number of following N. (occurrences(N) - 1 )
 * - bit 7 = 0   bit 6 = 0   : NT symbol is in bits 0 - 3 (see coding table for NT sequences)
 *
 * | Symbol       | Description                   | Bases represented  |    | Complement |Code|
 * |--------------|-------------------------------|--------------------|----|------------|----|
 * | A            | Adenine                       | A                  | 1  | T          | 0  |
 * | C            | Cytosine                      |     C              |    | G          | 1  |
 * | G            | Guanine                       |         G          |    | C          | 2  |
 * | T            | Thymine                       |             T      |    | A          | 3  |
 * | U            | Uracil                        |             U      |    | A          |    |
 * | W            | Weak                          | A           T      | 2  | W          |    |
 * | S            | Strong                        |     C   G          |    | S          |    |
 * | M            | aMino                         | A   C              |    | K          |    |
 * | K            | Keto                          |         G   T      |    | M          |    |
 * | R            | puRine                        | A       G          |    | Y          |    |
 * | Y            | pYrimidine                    |     C       T      |    | R          |    |
 * | B            | not A (B comes after A)       |     C   G   T      | 3  | V          |    |
 * | D            | not C (D comes after C)       | A       G   T      |    | H          |    |
 * | H            | not G (H comes after G)       | A   C       T      |    | D          |    |
 * | V            | not T (V comes after T and U) | A   C   G          |    | B          |    |
 * | N            | any Nucleotide (not a gap)    | A   C   G   T      | 4  | N          | 4  |
 * | Z            | Zero                          |                    | 0  | Z          |    |
 */
class CompressedNucSeq
{
    /** @brief Checks if next NUM_REQUESTED symbols starting at pUncomp are all A,C,G,T. uiSize is the number of
     *   remaining symbols in the uncompressed sequence.
     *  @details Template encourages loop unrolling.
     */
    template <size_t NUM_REQUESTED> inline size_t countAcgt( uint8_t* pUncomp, size_t uiSize )
    {
        size_t uiAcgtCount = 0;
        for( ; uiAcgtCount < uiSize; uiAcgtCount++ )
        {
            if( ( pUncomp[ uiAcgtCount ] > 3 ) // is not A,C,G or T
                || ( uiAcgtCount >= NUM_REQUESTED ) ) // exceeded maximum
                break;
            // DEL: else
            // DEL:     std::cout << uiAcgtCount << std::endl;
        }
        assert( uiAcgtCount <= NUM_REQUESTED );
        return uiAcgtCount;
    } // method

    /** @brief Counts and returns the number of N-symbols starting from of pUncompItr. The counting is limited by
     *  MAX_COUNT.
     */
    template <size_t MAX_COUNT> inline size_t countAnyNucleotide( uint8_t* pNucSecBuf, size_t uiSize )
    {
        size_t uiNCount = 0;
        for( ; uiNCount < uiSize; uiNCount++ )
            if( ( pNucSecBuf[ uiNCount ] != 4 ) // is not N
                || ( uiNCount >= MAX_COUNT ) ) // exceeded maximum
                break;
        return uiNCount;
    } // method

#if 0
    /** @brief Returns uiNuc shifted by SHIT bits in uiDest. */
    template <size_t SHIFT> inline uint8_t acgtShift( const uint8_t uiDest, const uint8_t uiNuc )
    {
        assert( uiNuc < 4 ); // must be ACGT
        return uiDest | ( uiNuc << SHIFT );
    } // method

    /** @brief Compresses 3 consecutive A,C,G,T into the bits 0-6 of a byte, shifted by SHIFT bits. */
    template <size_t SHIFT> inline uint8_t acgt3Compress( uint8_t* pNucSecBuf )
    {
        return acgtShift<0 + SHIFT>(
            acgtShift<2 + SHIFT>( acgtShift<4 + SHIFT>( 0x00, *pNucSecBuf ), *( pNucSecBuf + 1 ) ),
            *( pNucSecBuf + 2 ) );
    } // method
#else
    /** @brief Returns uiNuc shifted by SHIT bits in uiDest. */
    inline uint8_t acgtShift( const uint8_t uiDest, const uint8_t uiNuc, const uint8_t SHIFT )
    {
        assert( uiNuc < 4 ); // must be ACGT
        return uiDest | ( uiNuc << SHIFT );
    } // method

    /** @brief Compresses 3 consecutive A,C,G,T into the bits 0-6 of a byte, shifted by SHIFT bits. */
    template <size_t SHIFT> inline uint8_t acgt3Compress( uint8_t* pNucSecBuf )
    {
        return acgtShift( acgtShift( acgtShift( 0x00, *pNucSecBuf, 4 + SHIFT ), *( pNucSecBuf + 1 ), 2 + SHIFT ),
                          *( pNucSecBuf + 2 ), 0 + SHIFT );
    } // method
#endif

    /** @brief Compresses 4 consecutive A,C,G,T into a single byte. */
    inline uint8_t acgt4Compress( uint8_t* pNucSecBuf )
    {
        return acgtShift( acgt3Compress<2>( pNucSecBuf ), *( pNucSecBuf + 3 ), 0 );
    } // method

    /** @brief Write header to buffer. The header is the size of the uncompressed sequence as uint64 in big-endian. */
    inline void writeHeader( uint8_t* pComprBuf, uint64_t uiDecompSeqSize )
    {
        memcpy( pComprBuf, &uiDecompSeqSize, sizeof uiDecompSeqSize );
        // little-endian to big-endian on Intel and friends
        std::reverse( pComprBuf, pComprBuf + sizeof uiDecompSeqSize );
    } // method

    /** @brief Reads the size of the uncompressed sequence from the header. */
    inline uint64_t readHeader( uint8_t* pComprBuf )
    {
        uint64_t uiSize;
        // little-endian to big-endian on Intel and friends
        std::reverse( pComprBuf, pComprBuf + sizeof( uint64_t ) );
        memcpy( &uiSize, pComprBuf, sizeof uiSize );
        return uiSize;
    } // method

    /** @brief Decompress the bits 0-5 of *rpComp into three symbols in rpUncomp.
     *  Move the pointers forward accordingly.
     *  IMPROVEMENT: Via a table based approach, the decompression times could be improved.
     */
    inline void acgt3decompress( uint8_t*& rpComp, uint8_t*& rpUncomp )
    {
        auto uiFstByte = *( rpComp++ );
        *( rpUncomp + 2 ) = uiFstByte & ( 0x03 );
        uiFstByte = uiFstByte >> 2;
        *( rpUncomp + 1 ) = uiFstByte & ( 0x03 );
        uiFstByte = uiFstByte >> 2;
        *( rpUncomp + 0 ) = uiFstByte & ( 0x03 );
        rpUncomp += 3;
    } // method

    /** @brief Decompress the bits 0-7 of *rpComp into four symbols in rpUncomp.
     *  Move the pointers forward accordingly.
     *  IMPROVEMENT: Via a table based approach, the decompression times could be improved.
     */
    inline void acgt7decompress( uint8_t*& rpComp, uint8_t*& rpUncomp )
    {
        acgt3decompress( rpComp, rpUncomp );
        auto uiSecByte = *( rpComp++ );
        *( rpUncomp + 3 ) = uiSecByte & ( 0x03 );
        uiSecByte = uiSecByte >> 2;
        *( rpUncomp + 2 ) = uiSecByte & ( 0x03 );
        uiSecByte = uiSecByte >> 2;
        *( rpUncomp + 1 ) = uiSecByte & ( 0x03 );
        uiSecByte = uiSecByte >> 2;
        *( rpUncomp + 0 ) = uiSecByte & ( 0x03 );
        rpUncomp += 4;
    } // method

  public:
    size_t uiSizeCompSeq = 0; // size of final compressed sequence after compressing
    ByteBuffer xCompSeqBuf; // Internal buffer for storage of the compressed sequence.

    CompressedNucSeq( )
    {}

    CompressedNucSeq( const CompressedNucSeq& ) = delete; // no copies of compressed NucSeqs

    // CompressedNucSeq(const CompressedNucSeq&& other )
    // {
    //     std::cout << "MOVE CONSTRUCTOR" << std::endl;
    // } // move constructor

    CompressedNucSeq& operator=( CompressedNucSeq&& other )
    {
        // this may be called once or twice
        // if called twice, 'other' is the just-moved-from V subobject
        std::cout << "Used Move assignment" << std::endl;
        return *this;
    }

    /** @brief Returns the size of the compressed sequence. */
    inline size_t size( ) const
    {
        return uiSizeCompSeq;
    } // method

    /** @brief Returns a pointer to the internal buffer. */
    inline uint8_t* get( ) const
    {
        return (uint8_t*)xCompSeqBuf.get( );
    } // method

    /** @brief Compresses the nucleotide sequence given as argument and stores the result internally in the current
     *  object.
     */
    void compress( const NucSeq& rxNucSeq )
    {
        size_t uiReqBufSize = rxNucSeq.uiSize // size of nucleotide seq
                              + sizeof( uint64_t ); // storage for seq size memorizing

        xCompSeqBuf.resize( uiReqBufSize ); // compressed seq won't become longer than rxNucSeq itself
        const auto pCompStart = (uint8_t*)xCompSeqBuf.get( );
        auto pCompItr = pCompStart;
        writeHeader( pCompItr, rxNucSeq.uiSize );
        pCompItr += sizeof( uint64_t ); // forward the pointer after the header
        const uint8_t* pUncompEnd = rxNucSeq.pxSequenceRef + rxNucSeq.uiSize;
        // Iterate over the nucleotide sequence to be compressed
        for( uint8_t* pUncompItr = rxNucSeq.pxSequenceRef; pUncompItr < pUncompEnd; )
        {
            assert( pUncompItr - rxNucSeq.pxSequenceRef >= 0 );
            size_t uiSize = static_cast<size_t>( pUncompEnd - pUncompItr ); // remaining length of nuc seq
            size_t uiAcgtCount = countAcgt<7>( pUncompItr, uiSize );
            if( uiAcgtCount == 7 )
            {
                // Next seven symbols are all A,C,G,T
                // DEBUG: std::cout << "7";
                *( pCompItr++ ) = 0xC0 | acgt3Compress<0>( pUncompItr ); // write next 3 symbols
                *( pCompItr++ ) = acgt4Compress( pUncompItr + 3 ); // write next 4 symbols
                pUncompItr += 7; // move 7 symbols forward
            } // if (7 consecutive A,C,G<T)
            else if( uiAcgtCount >= 3 )
            {
                // AT least next three symbols are all A,C,G,T
                // DEBUG: std::cout << "3";
                *( pCompItr++ ) = 0x80 | acgt3Compress<0>( pUncompItr ); // write next 3 symbols
                pUncompItr += 3; // move 3 symbols forward
            } // else if (3 consecutive A,C,G<T)
            else if( size_t uiNum = countAnyNucleotide<64>( pUncompItr, uiSize ) )
            {
                // Consecutive sequence of uiNum N symbols (maximal 64 symbols)
                // DEBUG:std::cout << "N(" << uiNum << ")";
                assert( uiNum <= 64 && ( uiNum > 0 ) );
                size_t numAnyNucleotideAhead = uiNum - 1;
                assert( numAnyNucleotideAhead >= 0 && numAnyNucleotideAhead < 64 );
                *( pCompItr++ ) = 0x40 | static_cast<uint8_t>( numAnyNucleotideAhead );
                pUncompItr += uiNum;
            } // else if (consecutive N)
            else
            {
                // Encode a single symbol
                // DEBUG: std::cout << "*";
                assert( *pUncompItr < 32 );
                *( pCompItr++ ) = *( pUncompItr++ );
            } // else (single symbol)
        } // for

        uiSizeCompSeq = pCompItr - pCompStart;
        // DEBUG: std::cout << "Compression rate: " << (double)uiSizeCompSeq / (double)rxNucSeq.uiSize << std::endl;
    } // method

    /** @brief Decompresses the internally stored compressed nucleotide sequence into the object given as argument.
     *  The sequence that is hold by the argument is overwritten by the outcome of the decompression.
     */
    void decompress( NucSeq& rxNucSeq, uint8_t* pCompItr = NULL, size_t uiExtBufSize = 0 )
    {
        // If there is specific external buffer for decomposition, use the internal one that is used for compression.
        if( pCompItr == NULL )
            pCompItr = (uint8_t*)xCompSeqBuf.get( );
        uint64_t uiNucSeqSize = readHeader( pCompItr );
#ifndef NDEBUG
        auto pCompEnd = pCompItr + ( ( uiExtBufSize > 0 ) ? uiExtBufSize : this->uiSizeCompSeq );
#endif
        pCompItr += sizeof( uint64_t );

        rxNucSeq.resize( uiNucSeqSize );
        uint8_t* pUncompItr = rxNucSeq.pxSequenceRef;
        const uint8_t* pUncompEnd = rxNucSeq.pxSequenceRef + rxNucSeq.uiSize;
        // DEBUG: std::cout << uiNucSeqSize << std::endl;
        // Iterate so long decompression has not finished
        while( pUncompItr < pUncompEnd )
        {
            assert( pUncompEnd > pUncompItr );
#ifndef NDEBUG
            size_t uiSize = static_cast<size_t>( pUncompEnd - pUncompItr ); // remaining length of nuc seq
#endif
            if( *pCompItr & 0x80 ) // first check for 2 bit compression
            {
                // 2 bit compression
                if( *pCompItr & 0x40 )
                {
                    // 7 consecutive A,C,G,T
                    assert( uiSize >= 7 );
                    acgt7decompress( pCompItr, pUncompItr ); // pointers are moved via reference
                } // if
                else
                {
                    // 3 consecutive A,C,G,T
                    assert( uiSize >= 3 );
                    acgt3decompress( pCompItr, pUncompItr ); // pointers are moved via reference
                }
            } // if (2 bit encoded)
            else if( *pCompItr & 0x40 ) // check for consecutive N
            {
                // consecutive N symbols
                size_t uiNumSyms =
                    *pCompItr & 0x3F; // least 6 bit encode the number of N's
                                      // DEBUG:std::cout << "N decompress with size: " << uiNumSyms << std::endl;
                uiNumSyms++; // the stored value is reduced by one for efficiency reasons.
                assert( uiNumSyms <= uiSize );
                for( size_t uiCount = 0; uiCount < uiNumSyms; uiCount++ )
                    *( pUncompItr++ ) = 4; // 4 encodes N
                pCompItr++;
            } // else if (consecutive N)
            else // top two bits are zero, so we have a single symbol in bits 0-6
            {
                // copy single symbol
                *( pUncompItr++ ) = *( pCompItr++ );
            } // else (single symbol)
        } // while
        assert( pCompItr == pCompEnd );
    } // method

    std::shared_ptr<NucSeq> pUncomNucSeq = nullptr;

    /** @brief Decompresses from the buffer given as argument into the shared NucSeq available as attribute */
    void decompress( uint8_t* pCompExtBuf, size_t uiExtBufSize )
    {
        this->pUncomNucSeq = std::make_shared<NucSeq>( ); // create shared pointer for decompression.
        this->decompress( *pUncomNucSeq.get( ), pCompExtBuf, uiExtBufSize );
    } // method
}; // class CompressedNucSeq

inline std::shared_ptr<CompressedNucSeq> makeSharedCompNucSeq( const NucSeq& rxNucSeq )
{
    std::shared_ptr<CompressedNucSeq> pCompdNucSeq = std::make_shared<CompressedNucSeq>( );
    pCompdNucSeq->compress( rxNucSeq );
    return pCompdNucSeq;
} // function

/** @brief like makeSharedCompNucSeq( const NucSeq& rxNucSeq ) but can deal with nullptrs */
inline std::shared_ptr<CompressedNucSeq> makeSharedCompNucSeq( const std::shared_ptr<NucSeq> pNucSeq )
{
    if( pNucSeq == nullptr )
        return nullptr;
    return makeSharedCompNucSeq( *pNucSeq );
} // function

} // namespace libMA

inline std::string buf_to_hex( char* pBuf, size_t uiSize )
{
    static const char hex_digits[] = "0123456789ABCDEF";

    std::string output;
    output.reserve( uiSize * 2 );
    for( auto pItr = pBuf; pItr < pBuf + uiSize; pItr++ )
    {
        uint8_t c = (uint8_t)*pItr;
        // std::cout << (int)(c >> 4) << " " << (int)(c & 15) << std::endl;
        output.push_back( hex_digits[ c >> 4 ] );
        output.push_back( hex_digits[ c & 15 ] );
        output.push_back( ' ' );
    }
    return output;
} // function


#ifdef WITH_MYSQL
/* Integration of shared pointers to CompressedNucSeq objects as data-type in the MySQL interface.
 */
using CompNucSeqSharedPtr = std::shared_ptr<libMA::CompressedNucSeq>;

// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string MySQLConDB::TypeTranslator::getSQLTypeName<CompNucSeqSharedPtr>( )
{
    return "LONGBLOB";
} // specialized method

// Part 2: Input arguments: Set the start of the blob (void *), size of the blob and type of the blob.
template <> inline void MySQLConDB::StmtArg::set( const CompNucSeqSharedPtr& rxCompSeq )
{
    if( rxCompSeq != nullptr )
    {
        // On MSVC the size of unsigned long is 32 bit merely ...
        if( rxCompSeq->size( ) > static_cast<size_t>( std::numeric_limits<unsigned long>::max( ) ) )
            std::runtime_error( "MySQLConDB::StmtArg::set: Overflow for rxCompSeq.size() " );
        unsigned long uiCompSeqSize = static_cast<unsigned long>( rxCompSeq->size( ) );
        this->uiLength = uiCompSeqSize;
        pMySQLBind->buffer_length = uiCompSeqSize;
        pMySQLBind->buffer = (void*)( rxCompSeq->get( ) );
        pMySQLBind->buffer_type = MYSQL_TYPE_LONG_BLOB; // this type must be equal to the type in Part 3.
    } // if
    else
    {
        pMySQLBind->buffer = NULL;
        pMySQLBind->buffer_type = MYSQL_TYPE_NULL; // this type must be equal to the type in Part 3.
    } // else
} // specialized method

// Part 3: Code for supporting query output:
//         1. Via the third argument of the call of init, set the MySQL datatype for your cell type.
//         2. Using storeVarSizeCel, fetch the blob from the byte-buffer of the cell.
template <>
struct /* MySQLConDB:: */ RowCell<CompNucSeqSharedPtr> : public /* MySQLConDB::*/ RowCellBase<CompNucSeqSharedPtr>
{
    inline void init( MYSQL_BIND* pMySQLBind, CompNucSeqSharedPtr* pCellValue, size_t uiColNum )
    {
        *pCellValue = std::make_shared<libMA::CompressedNucSeq>( );
        RowCellBase<CompNucSeqSharedPtr>::init( pMySQLBind, pCellValue, MYSQL_TYPE_LONG_BLOB, uiColNum );
    } // method

    // Decompress the nucleotide sequence directly from the buffer.
    inline void storeVarSizeCell( )
    {
        // DEBUG: std::cout << buf_to_hex( this->pVarLenBuf.get( ), this->uiLength ) << std::endl;
        if( !( this->is_null ) )
        {
            if( *pCellValue == nullptr )
                *pCellValue = std::make_shared<libMA::CompressedNucSeq>( );
            ( *pCellValue )->decompress( reinterpret_cast<uint8_t*>( this->pVarLenBuf.get( ) ), this->uiLength );
        }
        else
            *pCellValue = nullptr;
    } // method
}; // specialized class
#endif

#ifdef POSTGRESQL
/* Integration of shared pointers to CompressedNucSeq objects as data-type in the MySQL interface.
 */
using CompNucSeqSharedPtr = std::shared_ptr<libMA::CompressedNucSeq>;

// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string PostgreSQLDBCon::TypeTranslator::getSQLTypeName<CompNucSeqSharedPtr>( )
{
    return "bytea";
} // specialized method

#ifdef _MSC_VER
template <> inline std::string PostgreSQLDBCon::TypeTranslator::getSQLTypeName<CompNucSeqSharedPtr&>( )
{
    return PostgreSQLDBCon::TypeTranslator::getSQLTypeName<CompNucSeqSharedPtr>( );
} // specialized method
#endif

// Part 2: Input arguments: Set the start of the blob (void *), size of the blob and type of the blob.
template <> inline void PostgreSQLDBCon::StmtArg::set( const CompNucSeqSharedPtr& rxCompSeq )
{
    if( rxCompSeq != nullptr )
    {
        unsigned long uiCompSeqSize = static_cast<unsigned long>( rxCompSeq->size( ) );
        if( uiCompSeqSize > (size_t)std::numeric_limits<int>::max( ) )
            throw PostgreSQLError( "PG: Length of NucSeq exceeds maximum of type integer." );

        rpParamValue = (char*)( rxCompSeq->get( ) );
        riParamLength = static_cast<int>( uiCompSeqSize );
        riParamFormat = PG_BINARY_ARG;
    } // if
    else
    {
        rpParamValue = (char*)NULL;
        riParamLength = 0;
        riParamFormat = PG_BINARY_ARG;
    } // else
} // specialized method

// Part 3: Code for supporting query output:
//         1. Via the third argument of the call of init, set the MySQL datatype for your cell type.
//         2. Using storeVarSizeCel, fetch the blob from the byte-buffer of the cell.
template <> struct PGRowCell<CompNucSeqSharedPtr> : public PGRowCellBase<CompNucSeqSharedPtr>
{
    inline void init( CompNucSeqSharedPtr* pCellValue, size_t uiColNum )
    {
        *pCellValue = std::make_shared<libMA::CompressedNucSeq>( );
        PGRowCellBase<CompNucSeqSharedPtr>::init( pCellValue, uiColNum );
    } // method

    // Decompress the nucleotide sequence directly from the buffer.
    inline void store( const PGresult* pPGRes )
    {
        // DEBUG: std::cout << buf_to_hex( this->pVarLenBuf.get( ), this->uiLength ) << std::endl;
        if( !( this->isNull ) )
        {
            if( *pCellValue == nullptr )
                *pCellValue = std::make_shared<libMA::CompressedNucSeq>( );
            ( *pCellValue )
                ->decompress( reinterpret_cast<uint8_t*>( this->getValPtr( pPGRes ) ), this->getValLength( pPGRes ) );
        } // if
    } // method
}; // specialized class
#endif


/** @brief Implements the binary representation of compressed sequences */
template <typename DBCon>
inline std::string csvPrint( DBCon& rxDBCon, const std::shared_ptr<libMA::CompressedNucSeq>& pCompNucSep )
{
    return rxDBCon.blobAsQuotedSafeString( pCompNucSep->xCompSeqBuf.get( ),
                                           pCompNucSep->uiSizeCompSeq ); // mysql_real_escape_string_quote()
}; // method

namespace libMA
{
/** @brief Encapsulates a shared pointer to a nucleotide sequence */
class NucSeqSql // not required any longer : public SQL_BLOB
{
  public:
    std::shared_ptr<NucSeq> pNucSeq = nullptr;

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
        pNucSeq = std::make_shared<NucSeq>( );
        pNucSeq->vClear( );
        if( uiSize == 0 )
            return;
        pNucSeq->vAppend( (const uint8_t*)ucBlob, uiSize );
    } // method
}; // class

} // namespace libMA

#ifdef WITH_MYSQL
/* DATABASE INTEGRATION MySQL */

/* Integration of NucSeqSql as data-type in the MySQL interface.
 */
// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string MySQLConDB::TypeTranslator::getSQLTypeName<libMA::NucSeqSql>( )
{
    return "LONGBLOB";
} // specialized method

// Part 2: Input arguments: Set the start of the blob (void *), size of the blob and type of the blob.
template <> inline void MySQLConDB::StmtArg::set( const libMA::NucSeqSql& rxBlob )
{
    this->uiLength = static_cast<unsigned long>( rxBlob.blobSize( ) );
    pMySQLBind->buffer_length = static_cast<unsigned long>( rxBlob.blobSize( ) );
    pMySQLBind->buffer_type = MYSQL_TYPE_LONG_BLOB; // this type must be equal to the type in Part 3.
    pMySQLBind->buffer = (void*)( rxBlob.toBlob( ) );
} // specialized method

// Part 3: Code for supporting query output:
//         1. Via the third argument of the call of init, set the MySQL datatype for your cell type.
//         2. Using storeVarSizeCel, fetch the blob from the byte-buffer of the cell.
template <> struct /* MySQLConDB:: */ RowCell<libMA::NucSeqSql> : public /* MySQLConDB::*/ RowCellBase<libMA::NucSeqSql>
{
    inline void init( MYSQL_BIND* pMySQLBind, libMA::NucSeqSql* pCellValue, size_t uiColNum )
    {
        RowCellBase<libMA::NucSeqSql>::init( pMySQLBind, pCellValue, MYSQL_TYPE_LONG_BLOB, uiColNum );
    } // method

    // Fetch the blob from the buffer.
    inline void storeVarSizeCell( )
    {
        if( !( this->is_null ) )
            pCellValue->fromBlob( reinterpret_cast<unsigned char*>( this->pVarLenBuf.get( ) ), this->uiLength );
    } // method
}; // specialized class
#endif

#ifdef POSTGRESQL
/* DATABASE INTEGRATION PostgreSQL
 */

/* Integration of NucSeqSql as data-type in the MySQL interface.
 */
// Part1 : Specify the corresponding MySQL-type for your blob.
template <> inline std::string PostgreSQLDBCon::TypeTranslator::getSQLTypeName<libMA::NucSeqSql>( )
{
    return "bytea";
} // specialized method

#ifdef _MSC_VER
template <> inline std::string PostgreSQLDBCon::TypeTranslator::getSQLTypeName<libMA::NucSeqSql&>( )
{
    return PostgreSQLDBCon::TypeTranslator::getSQLTypeName<libMA::NucSeqSql>( );
} // specialized method
#endif

// Part 2: Input arguments: Set the start of the blob (void *), size of the blob and type of the blob.
template <> inline void PostgreSQLDBCon::StmtArg::set( const libMA::NucSeqSql& rxBlob )
{
    if( rxBlob.blobSize( ) > (size_t)std::numeric_limits<int>::max( ) )
        throw PostgreSQLError( "PG: Length of NucSeq exceeds maximum of type integer." );
    rpParamValue = (char*)( rxBlob.toBlob( ) );
    riParamLength = static_cast<int>( rxBlob.blobSize( ) );
    riParamFormat = PG_BINARY_ARG;
} // specialized method

// Part 3: Code for supporting query output:
template <> struct PGRowCell<libMA::NucSeqSql> : public PGRowCellBase<libMA::NucSeqSql>
{
    inline void init( libMA::NucSeqSql* pCellValue, size_t uiColNum )
    {
        PGRowCellBase<libMA::NucSeqSql>::init( pCellValue, uiColNum );
    } // method

    inline void store( const PGresult* pPGRes )
    {
        if( !( this->isNull ) )
            pCellValue->fromBlob( reinterpret_cast<unsigned char*>( this->getValPtr( pPGRes ) ),
                                  this->getValLength( pPGRes ) );
    } // method
}; // specialized class
#endif


#ifdef WITH_PYTHON
void exportNucSeq( libMS::SubmoduleOrganizer& xOrganizer );
#endif
