#pragma once

#include <memory>
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <numeric>
#include <cmath>
#include <cstring>
#include <functional>
#include <string>
#include "exception.h"
#include "support.h"
#include "interval.h"

#include <boost/log/trivial.hpp>
#include <boost/python.hpp>
#include "container.h"


/* 32bit rounding to the next exponent as define
 */
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/* Generic reverse function, as it occurs in std::algorithms
 */
template <class T>
void reverse(T word[], size_t length)
{
	char temp;
	for ( size_t i = 0; i < length / 2; i++ )
	{
		temp = word[i];
		word[i] = word[length - i - 1 ];
		word[length - i - 1] = temp;
	} // for
} // reverse

/* Class for the management of sequences (genetic or text)
 * Special string class, for sequence handling. 
 */
template <class ELEMENT_TYPE>
class PlainSequence
{
private :
	friend GeneticSequence;

	/* Normally we avoid the copy of PlainSequence Objects, except in the context of the GeneticSequence class.
	 */
	PlainSequence( const PlainSequence &rSequence )
	{
		/* We reset all attributes of our fresh sequence object.
		 */
		vResetProtectedAttributes();

		/* Now we copy all elements from the source sequence to the current sequence.
		 */
		vAppend( rSequence.pGetSequenceRef(), rSequence.uxGetSequenceSize() );
	} // copy constructor

	/* Sophisticated copy constructor that takes additionally and interval and predicate as argument.
	 * The predicate gets a relative index as input and has to decided whether the element of this index makes it into the result.
	 * TODO: range check for the interval.
	 */
	PlainSequence( const PlainSequence &rSequence, 
				   const DeprecatedIntervalDescriptor &rInterval, 
				   const std::function<bool (size_t)>& predciateFunction )
	{
		/* We reset all attributes of our fresh sequence object.
		 */
		vResetProtectedAttributes();

		/* We iterate over the complete interval and apply the predicate to each element of the interval.
		 */
		for (size_t uxIterator = rInterval.uxStart; uxIterator < rInterval.uxEnd; uxIterator++ )
		{
			if ( predciateFunction( uxIterator - rInterval.uxStart ) == true )
			{
				vAppend( rSequence.pxSequenceRef + uxIterator, 1 );
			} // if
		} // for
	} // copy constructor

	/* Copy constructor that copies some subsection indicated by the interval.
	 * The interval will not be ranged-checked!
	 */
	PlainSequence( const PlainSequence &rSequence, // the sequence that is taken as copy source
				   const ZeroBasedIntervalDescriptor &rInterval
				 )
	{
		/* We reset all attributes of our fresh sequence object.
		 */
		vResetProtectedAttributes();

		//// BOOST_LOG_TRIVIAL( warning ) << "rInterval : " << rInterval.uxSize << " " << rInterval.uxOffsetInHostSection << " " << rSequence.uiSize;

		/* We iterate over the complete interval and apply the predicate to each element of the interval.
		 */
		for ( size_t uxIterator = 0; uxIterator < rInterval.uxSize; uxIterator++ )
		{
			vAppend( rSequence.pxSequenceRef + ( rInterval.uxOffsetInHostSection + uxIterator ), 1 );
		} // for
	} // copy constructor

protected :
	/* The encapsulated sequence
	 */
	ELEMENT_TYPE *pxSequenceRef;

	/* Current size of the content of the encapsulated sequence
	 */
	size_t uiSize;

	/* Current size of the buffer.
	 */
	size_t uxCapacity;

	/* Resets all protected attributes to its initial values.
	 */
	inline void vReleaseMemory()
	{
		/* Allocated memory will be released!
		 */
		if ( pxSequenceRef != NULL ) 
		{
			free( pxSequenceRef );
		} // if
	} // protected method

	void vResetProtectedAttributes()
	{
		pxSequenceRef = NULL;
		uiSize = 0;
		uxCapacity = 0;
	} // protected method
	
	/* Tries to allocate the requested amount of memory and throws an exception if this process fails.
	 * uxRequestedSize is expressed in "number of requested elements.
	 */
	void vReserveMemory( size_t uxRequestedSize )
	{
		/* TO DO: This should be a bit more sophisticated ...
		 */
		kroundup32( uxRequestedSize );
		
		/* We try to reserve the requested memory.
		 * See: http://stackoverflow.com/questions/1986538/how-to-handle-realloc-when-it-fails-due-to-memory
		 */
		auto pxReallocRef = (ELEMENT_TYPE*)realloc( pxSequenceRef, uxRequestedSize * sizeof(ELEMENT_TYPE) );

		if ( pxReallocRef == NULL )
		{
			throw fasta_reader_exception( (std::string( "Memory Reallocation Failed for requested size ") + std::to_string( uxRequestedSize )).c_str() );
		} // if

		pxSequenceRef = pxReallocRef;
		uxCapacity = uxRequestedSize;
	} // method
	
public :
	PlainSequence() 
	{
		vResetProtectedAttributes();
	} // default constructor

	virtual ~PlainSequence()
	{
		/* Release all allocated memory.
		 */
		vReleaseMemory();
	} // destructor

	/* This moves the ownership of the protected attributes to another object.
	 * The receiver of pxSequenceRef is responsible for its deletion.
	 */
	void vTransferOwnership( PlainSequence &rReceivingSequence )
	{
		/* We transport the three protected attributes to the receiver ...
		 */
		rReceivingSequence.pxSequenceRef = this->pxSequenceRef;
		rReceivingSequence.uiSize = this->uiSize;
		rReceivingSequence.uxCapacity = this->uxCapacity;

		/* ... and delete the knowledge here
		 */
		vResetProtectedAttributes();
	} // protected method

	/* Clears the inner sequence, but does not deallocate the memory.
	 */
	inline void vClear()
	{
		uiSize = 0;
	} // method

	/* Returns whether the sequence is empty or not.
	 */
	inline bool bEmpty()
	{
		return uiSize == 0;
	} // method

	inline bool empty() const
	{
		return uiSize == 0;
	} // method

	/* Fast getter and setter for element access.
	 * If assertions activated we do a range check.
	 */
	inline ELEMENT_TYPE operator[]( size_t uiSubscript ) const
	{
		assert( uiSubscript < uiSize );
		return pxSequenceRef[uiSubscript];
	} // method (get)
	inline ELEMENT_TYPE & operator[]( size_t uiSubscript )
	{
		assert( uiSubscript < uiSize );
		return pxSequenceRef[uiSubscript];
	} // method (set)

	/* Resizes the internal buffer of the sequence to the requested value.
	 */
	inline void resize( size_t uiRequestedSize ) // throws exception
	{	/* Check, whether we have enough capacity, if not reserve memory
		 */
		if ( uxCapacity < uiRequestedSize )
		{
			vReserveMemory( uiRequestedSize );
		} // if
		
		uiSize = uiRequestedSize;
	} // method

	/* Because we want the reference to the sequence private we offer a getter method.
	 * WARNING! Here you can get a null-pointer.
	 */
	inline const ELEMENT_TYPE* const pGetSequenceRef() const
	{
		return this->pxSequenceRef;
	} // method

	/* Because we want to keep the size private we offer a getter method.
	 */
	inline const size_t uxGetSequenceSize() const
	{
		return this->uiSize;
	} // method

	inline const size_t length() const
	{
		return this->uiSize;
	} // method

	/* Reverse the elements of the plain sequence.
	 */
	inline void vReverse()
	{
		reverse( pxSequenceRef, uiSize );
	} // method
	
	/* WARNING: the inner string might not null-terminated after this operation.
	 */
	inline PlainSequence& vAppend( const ELEMENT_TYPE* pSequence, size_t uxNumberOfElements )
	{
		size_t uxRequestedSize = uxNumberOfElements + this->uiSize;

		if ( uxCapacity < uxRequestedSize )
		{
			vReserveMemory ( uxRequestedSize );
		} // if

		/* WARNING: If we work later with non 8-bit data we have to be careful here
		 */
		memcpy( this->pxSequenceRef + uiSize, pSequence, uxNumberOfElements * sizeof(ELEMENT_TYPE) );

		uiSize = uxRequestedSize;

		return *this;
	} // method

	/* Push back of a single symbol.
	 */
	inline void push_back( const ELEMENT_TYPE xElement )
	{
		if ( this->uiSize >= this->uxCapacity )
		{
			vReserveMemory( this->uiSize + 1 );
		} // if

		pxSequenceRef[uiSize + 1] = xElement;
		uiSize++;
	} // method

	/* Compares two sequences for equality
	 */
	inline bool equal(const PlainSequence &rOtherSequence)
	{
		if ( this->uiSize == rOtherSequence.uiSize )
		{
			return memcmp(this->pxSequenceRef, rOtherSequence.pxSequenceRef, sizeof(ELEMENT_TYPE) * uiSize ) == 0;
		} // if
		
		return false;
	} // method
}; // class PlainSequence

/* This call was exclusively build for the fasta-reader.
 * It shall boost performance for long inputs.
 */
class TextSequence : public PlainSequence<char>
{
public :
	TextSequence() : PlainSequence<char>()
	{
	} // default constructor

	TextSequence( const char* pcString ) : PlainSequence<char>()
	{
		vAppend(pcString );
	} // text constructor

	/* Terminates the inner string in the C-style using a null-character and 
	 * returns a reference to the location of the inner string.
	 */
	inline char* cString()
	{
		if ( uxCapacity < this->uiSize + 1 )
		{
			vReserveMemory ( this->uiSize + 1 );
		} // if

		this->pxSequenceRef[this->uiSize] = '\0';
		
		return pxSequenceRef;
	} // method

	inline void vAppend( const char &pcChar )
	{
		if ( uxCapacity < ( this->uiSize + 1 ) )
		{
			vReserveMemory ( this->uiSize + 1 );
		} // if

		this->pxSequenceRef[this->uiSize] = pcChar;
		this->uiSize++;
	} // method

	/* Appends the content of pcString to the current buffer
	 */
	inline void vAppend( const char* pcString )
	{
		PlainSequence<char>::vAppend( pcString, strlen( pcString ) );
	} // method
}; // class

/* Forward declarations
 */
class GeneticSequenceSlice;
class DoubleStrandedSequence;

/* Special Class for Genetic Sequences
 * IDEA: BioSequence objects use numbers instead of characters for sequence representation
 * Supports:
 *  - translation from textual representation to representation as sequence of numbers.
 *  - generation of the reverse strand 
 */

class GeneticSequence : public PlainSequence<uint8_t>
{
public :
	/* The type of elements represented by our sequence.
	 * (the type has to be decided in the context of the construction)
	 */
	// const SequenceType eContentType; 

	/* Disable the copy constructor.
	 */
	GeneticSequence( const GeneticSequence& ) = delete;

	GeneticSequence( ) // : eContentType( SEQUENCE_IS_NUCLEOTIDE )
	{
	} // default constructor

	GeneticSequence (GeneticSequence && g) // : eContentType( SEQUENCE_IS_NUCLEOTIDE )
	{
	} // default constructor

	/* Deprecated constructor, don't use it any longer!.
	 * Takes all elements of the interval for which the predicate is true.
	 * The predicate's input is the relative index within the interval.
	 */
	GeneticSequence( const PlainSequence &rSequence, const DeprecatedIntervalDescriptor &rInterval, const std::function<bool (size_t)>& predciateFunction )
		: PlainSequence( rSequence, rInterval, predciateFunction ) // forward call to the copy constructor of base class
		  // eContentType( SEQUENCE_IS_NUCLEOTIDE )
	{
	} // tailored copy constructor

	/* Copies a subsection of the current sequence on the foundation of the given interval.
	 * Used e.g. by sequence slices.
	 */
	GeneticSequence( const PlainSequence &rxSequence, const ZeroBasedIntervalDescriptor &rxInterval )
		: PlainSequence( rxSequence, rxInterval ) // forward call to the copy constructor of base class
		// eContentType( SEQUENCE_IS_NUCLEOTIDE )
	{} // tailored copy constructor

	/* The full sequence as slice.
	 */
	GeneticSequenceSlice fullSequenceAsSlice() const;

	/* A section of the sequence as slice
	 */
	GeneticSequenceSlice fromTo( size_t from, size_t end );

	/* Get a section by means of an interval.
	 */
	GeneticSequenceSlice getSliceByInterval( const DeprecatedIntervalDescriptor &interval ) const;

	/* TO DO: Make the 5 a class constant!
	 */
	inline uint8_t uxAlphabetSize() const
	{
		/* The Alphabet size for sequences of nucleotides is 5.
		 */
		return 5; // eContentType == SEQUENCE_IS_NUCLEOTIDE ? 5 : 20;
	} // method

	/* Appends a slice of some other sequence to our current sequence. 
	 * Externally defined.
	 */
	void vAppend( const GeneticSequenceSlice & );
}; // class

class NucleotideSequenceSliceHavingSelectedStrand;

/* Class for genetic sequence that consist of nucleotides. (A, C, G, T)
 * TO DO: Create a move constructor!
 */
class NucleotideSequence : public GeneticSequence, public Container
{
private :

public :
	/* The table used to translate from base pairs to numeric codes for nucleotides
	 */
	static const unsigned char xNucleotideTranslationTable[256];

	/* Forward call to the constructor of the superclass
	 */
	NucleotideSequence( const PlainSequence &rxSequence, const ZeroBasedIntervalDescriptor &rxInterval )
		: GeneticSequence( rxSequence, rxInterval )
	{ } // constructor

	/* Default constructor
	 */
	NucleotideSequence()
		: GeneticSequence()
	{ } // default constructor

	/* Constructor that get the initial content of the sequence in text form.
	 * FIX ME: This can be done a bit more efficient via the GeneticSequence class.
	 */
	NucleotideSequence( const std::string &rsInitialText )
		: GeneticSequence()
	{
		vAppend( rsInitialText.c_str() );
	} // constructor

	/* is implicitly deleted by geneticSequence but boost python needs to know */
	NucleotideSequence(const NucleotideSequence&) = delete;

	/*used to identify the nucleotide sequence datatype in the aligner pipeline*/
    ContainerType getType(){return ContainerType::nucSeq;}

	/* The full sequence as nucleotide sequence slice.
	 */
	NucleotideSequenceSlice fullSequenceAsSlice() const;

	/* The full sequence as slice having additionally selected forward or reverser strand.
	 */
	NucleotideSequenceSliceHavingSelectedStrand fullSequenceAsSliceHavingSelectedStrand
		( bool = false ) const;

	/* Overwritten form of the base class method
	 */
	NucleotideSequenceSlice getSliceByInterval( const DeprecatedIntervalDescriptor &interval ) const;

	NucleotideSequenceSliceHavingSelectedStrand getSliceByInterval( const DeprecatedDoubleStrandIntervalDescriptor &interval ) const;

	/* Move constructor on the foundation of text sequences.
	 * Reuses the space of the text-sequence! TO DO: move & to &&
	 */
	NucleotideSequence( TextSequence &rSequence ) 
	{
		/* We strip the given sequence of its content and move it to our new sequence.
		 * WARNING: Here we assume that the sizes for the types char and uint8_t are equal.
		 */
		rSequence.vTransferOwnership( (PlainSequence<char>&)(*this) );

		/* The given PlainSequence should be in textual, we have to translate it.
		 */
		vTranslateToNumericFormUsingTable( xNucleotideTranslationTable, 0 );
	} // constructor


	/* Delivers the complement of a single nucleotide.
	 */
	static inline char nucleotideComplement( char iNucleotide )
	{
		/* Complements of nucleotides
		 *							   0  1  2  3
		 */
		static const char chars[4] = { 3, 2, 1, 0 };

		return ( iNucleotide < 4 ) ? chars[(int)iNucleotide] : 5;
	} // static method

	/* A section of the sequence as slice.
	 * TO DO: This method should get a sisters fromToCopy, fromToReverseCopy
	 */
	NucleotideSequenceSlice fromTo(size_t from, size_t end) const;
	/* by markus
	*returns a sequence slice of the index from while checking if from >= 0
	*		with the length of uiSize while checking if the end of the slice is bigger than sequence length
	*	if the slice is out of bounds a slice of the size 0 is returned
	*	uses fromTo
	*/
	NucleotideSequenceSlice saveFromSize(long uiFrom, size_t uiSize) const;

	/* Iterates over all base pairs in the sequence and creates the complement. 
	 * (A -> T, T -> A, C -> G, G -> C)
	 */
	void vSwitchAllBasePairsToComplement()
	{
		for( size_t uxIterator = 0; uxIterator < uiSize; uxIterator++ )
		{
			pxSequenceRef[uxIterator] = nucleotideComplement( pxSequenceRef[uxIterator] );
		} // for
	} // function

/*	void vReverse()
	{
		PlainSequence.vReverse();
	}*/

	/* transforms the character representation into a representation on the foundation of digits.
	 */
	void vTranslateToNumericFormUsingTable( const unsigned char *alphabetTranslationTable,
											size_t uxStartIndex
										  )
	{
		for( size_t uxIterator = uxStartIndex; uxIterator < uiSize; uxIterator++ )
		{
			pxSequenceRef[uxIterator] = alphabetTranslationTable[ pxSequenceRef[uxIterator] ];
		} // for
	} // method

	/* Gives the textual representation for some numeric representation.
	 * Important: Keep this inline, so that it is not compiled into a function of its own. 
	 */
	static inline char translateACGTCodeToCharacter( uint8_t uiNucleotideCode )
	{
		static const char chars[4] = {'A', 'C', 'G', 'T'};
		if (uiNucleotideCode < 4)
		{
			return chars[uiNucleotideCode];
		} // if
		else
		{
			return 'N';
		} // else
	} // static method

	/* The symbol on some position in textual form.
	 * We count starting from 0.
	 */
	inline char charAt( size_t uxPosition )
	{
		if ( uxPosition >= uiSize)
		{
			throw fasta_reader_exception("Index out of range");
		} // if

		return translateACGTCodeToCharacter( pxSequenceRef[uxPosition] );
	} // method

	/* Appends a string containing nucleotides as text and automatically translates the symbols.
	 */
	void vAppend( const char* pcString )
	{
		size_t uxSizeBeforeAppendOperation = this->uiSize;
		
		/* WARNING! char and uint8_t must have the same size or we get a serious problem here!
		 */
		PlainSequence<uint8_t>::vAppend( (const uint8_t*)pcString, strlen( pcString ) );

		vTranslateToNumericFormUsingTable( xNucleotideTranslationTable, uxSizeBeforeAppendOperation );
	} // method

	/* wrapper for boost
	 */
	void vAppend_boost( const char* pcString )
	{
		vAppend(pcString);
	} // method

	/* Appends a slice of some other sequence to our current sequence. (Only Prototype)
	 */
	void vAppend( const NucleotideSequenceSlice & );

	/* Appends a slice that contains additional information regarding the strand. (Only Prototype)
	 */
	void vAppend( const NucleotideSequenceSliceHavingSelectedStrand & );

	/* Catenate all intervals of the interval vector.
	 * Write the outcome in the sequence given as second argument. 
	 */
	void catenate( const DeprecatedIntervalVector<DeprecatedDoubleStrandIntervalDescriptor> &rxIntervalVector,
				   NucleotideSequence &rxNucleotideSequence
				 ) const;
}; // class NucleotideSequence

/* Takes a sequence slice and creates a vector with virtual contigs
 * TYPE should be GeneticSequenceSlice, NucleotideSequenceSlice or NucleotideSequenceSliceHavingSelectedStrand
 * (Note: the trick here is the automatic selection of the appropriate form of getSubSliceFromTo by means of TYPE)
 * sequence digestion 
 */
template<typename TYPE>
std::unique_ptr< std::vector< TYPE > > genericCreateVectorOfPseudoContigs  
(
	const TYPE &rSequenceSlice,
	const unsigned int uiContigSize, // size of each contig
	const unsigned int uiShift, // shifting of the start position between two contigs
	const unsigned int uiAmplificationRate // number of copies create of some single contig
) 
{
	std::unique_ptr< std::vector< TYPE > > pxReturnedVectorRef ( new std::vector< TYPE >() );
	
	for ( unsigned int uiCounter = 0; uiCounter < uiAmplificationRate; uiCounter++ )
	{
		size_t uiContigStartIndex = 0;

		/* Iterate by continuously shifting.
		 */
		while ( uiContigStartIndex + uiContigSize + (uiContigSize / 2) < rSequenceSlice.xInterval.uxSize )
		{
			pxReturnedVectorRef->emplace_back( rSequenceSlice.getSubSliceFromTo( uiContigStartIndex, uiContigStartIndex + uiContigSize ) );

			uiContigStartIndex += uiShift;
		} // while

		pxReturnedVectorRef->emplace_back( rSequenceSlice.getSubSliceFromTo( uiContigStartIndex, rSequenceSlice.xInterval.uxSize ) );
	} // for 

	return pxReturnedVectorRef;
} // method

/* This function is, like asString() quite inefficient, because it can trigger repeated string copy operations.
 * We put this outside of the nucleotide class, in order to get correct static bindings on compile time.
 * FIX ME: Check, whether there is some trick for getting this back into the class.
 */
template<typename TYPE>
std::string genericSequenceAsFastaText( const TYPE &rSequenceSlice,
										const char* pcFastaNameRef 
									  )
{
	BOOST_LOG_TRIVIAL(trace) << "FASTA creation";
	/* We start with the symbol '>' ...
	 */
	std::string psFastaTextRef = ">";

	/* ... then we append the fasta name ...
	 */
	psFastaTextRef.append( pcFastaNameRef ).append("\n");

	/* ... the sequence itself ...
	 */
	psFastaTextRef.append( *rSequenceSlice.asSequenceOfACGT() );

	/* and some final carriage return.
	 */
	psFastaTextRef.append( "\n" );

	return psFastaTextRef;
} // method

/* Forward declaration
 */
class NucleotideSequenceSlice;


/* Objects of this class represent slices of sequences that contain genetic data.
 * For such slices are for reading purposes merely. 
 * TO DO: represent GeneticSequenceSlice on foundation of an interval (maybe as derived class.)
 */
template<typename TYPE>
class GenericGeneticSequenceSlice
{
public :
	/* Remark: In the moment this works only for genetic sequences
	 * WARNING!: This reference construct can be come dangerous in the context of the pooling.
	 * Here we should work with shared pointers in order to avoid problems.
	 */
	const GeneticSequence &hostObject;

	/* The interval that specifies begin and end of the slice.
	 */
	const ZeroBasedIntervalDescriptor xInterval;

public :
	friend class NCBIGeneLocusDescriptor;

	/* uxStartOffsetInHostSection is inclusive,
	 * uxEndOffsetInHostSection is exclusive, i.e. first element after the end of the sequence
	 */
	GenericGeneticSequenceSlice<TYPE>( const GeneticSequence &hostObject, 
									   size_t uxStartOffsetInHostSection,
									   size_t uxEndOffsetInHostSection
									 )
		 : hostObject( hostObject ), 
		   xInterval( uxStartOffsetInHostSection, uxEndOffsetInHostSection - uxStartOffsetInHostSection )
	{	/* Do some consistency checks for the input parameter.
		 */
		if ( uxEndOffsetInHostSection > hostObject.length( ) )
		{
			throw std::runtime_error( "GeneticSequenceSlice creation failed because slice end > hostObject.uiSize. (start : "
									  + std::to_string( uxStartOffsetInHostSection ) + " hostObject.uiSize : " + std::to_string( hostObject.length( ) ) + " )"
									  );
		} // if
		if ( uxStartOffsetInHostSection > uxEndOffsetInHostSection )
		{
			throw std::runtime_error( "GeneticSequenceSlice creation failed because start > end" );
		} // if
	} // constructor

	/* Temporary bridge to the new class NucleotideSequenceSlice.
	 * Should be later removed, if we switched to NucleotideSequence.
	 */
	// NucleotideSequenceSlice castToNucleotideSequenceSlice();

	/* Because we want the reference to the sequence private we offer a getter method.
	 * WARNING! Here you can get a null-pointer.
	 */
	inline const uint8_t* getSequenceRefInHost() const 
	{
		return hostObject.pGetSequenceRef() + xInterval.uxOffsetInHostSection;
	} // method

	/* Because we want to keep the size private we offer a getter method.
	 */
	inline const size_t getSize() const
	{
		return this->xInterval.uxSize;
	} // method

	/* Get a subsection of the current slice.
	 * The method is generic and adapts its return type.
	 */
	TYPE getSubSliceFromTo( size_t uxStartOffsetInSlice, size_t uxEndOffsetInSlice ) const
	{
		// std::cout << "Offset " << uxOffsetInHostSection << "\n";
		// std::cout << "Make slice " << uxOffsetInHostSection + uxStartOffsetInSlice << " to " << uxOffsetInHostSection + uxEndOffsetInSlice << "\n";
		return TYPE( hostObject, xInterval.uxOffsetInHostSection + uxStartOffsetInSlice, xInterval.uxOffsetInHostSection + uxEndOffsetInSlice );
	} // method

	/* Copy the slice as genetic sequence using forward order.
	 */
	std::shared_ptr<GeneticSequence> makeSequenceUsingForwardOrder() const
	{
		return std::make_shared<GeneticSequence>( hostObject, xInterval );
	} // method

	/* Copy the slice as genetic sequence using reverse order.
	 */
	std::shared_ptr<GeneticSequence> makeGenericSequenceUsingReverseOrder() const
	{
		auto xCopiedSequence = makeSequenceUsingForwardOrder();
		xCopiedSequence->vReverse();
		return xCopiedSequence;
	} // method

	/* Checks for equality of the content of two sequence slices.
	 */
	bool equals( const GenericGeneticSequenceSlice<TYPE> &rOtherSequenceSlice )
	{
		if ( this->xInterval.uxSize == rOtherSequenceSlice.xInterval.uxSize )
		{
			return std::memcmp( (void *)( hostObject.pGetSequenceRef() + xInterval.uxOffsetInHostSection ), 
								(void *)( rOtherSequenceSlice.hostObject.pGetSequenceRef() + rOtherSequenceSlice.xInterval.uxOffsetInHostSection ), 
								(size_t)(sizeof(uint8_t) * this->xInterval.uxSize) 
							  ) == 0;
		} // if 

		return false;
	} // method


	/* Creates a histogram of the occurrence of single symbols in a sequence.
	 * The output is delivered in form of some vector that has the size of the sequences alphabet
	 */
	std::shared_ptr<std::vector<unsigned int>> pCreateHistogramVector( /* bool (*predicateFunction) (size_t) */ ) const
	{
		/* We make a histogram vector with a size according to the size of the sequence alphabet.
		 */
		auto pHistogramVectorRef = std::make_shared<std::vector<unsigned int>>( hostObject.uxAlphabetSize() );

		/* We iterate over the sequence
		 */
		for( size_t uxIterator = 0; uxIterator < this->xInterval.uxSize; uxIterator++ )
		{
			(*pHistogramVectorRef)[ getSequenceRefInHost()[uxIterator] ] += 1;
		} // for

		return pHistogramVectorRef;
	} // method
}; // templated class

class GeneticSequenceSlice : public GenericGeneticSequenceSlice<GeneticSequenceSlice> 
{
public :
	/* Forward of constructor
	 */
	GeneticSequenceSlice( const GeneticSequence &hostObject,
						  size_t uxStartOffsetInHostSection,
					      size_t uxEndOffsetInHostSection
						)
		: GenericGeneticSequenceSlice<GeneticSequenceSlice>( hostObject, uxStartOffsetInHostSection, uxEndOffsetInHostSection)
	{ } // constructor
}; // class

/* General form of some slice(subsection) of a nucleotide sequence. 
 * So, we have the info that the slice belongs to a nucleotide sequence, but we know nothing with respect to the selected strand.
 */
template<typename TYPE>
class GenericNucleotideSequenceSlice : public GenericGeneticSequenceSlice<TYPE> // public GeneticSequenceSlice
{
private :
	/* subfunction of NucleotideStatisiticAsString()
	 */
	void vAppendValueToString( std::string &sOutputText, const std::string &sSymbol, unsigned int uxNumberOccurences ) const
	{
		sOutputText.append( sSymbol )
			.append( " : " )
			.append( numberToStringDeprecated( uxNumberOccurences ) )
			.append( " (% " )
			.append( numberToStringDeprecated( floor( (double)(uxNumberOccurences * 10000)
			/ (double)GenericGeneticSequenceSlice<TYPE>::xInterval.uxSize
			) / 100
			)  // numberToString
			) // append
			.append( ") " );
	} // method

public :
	/* Constructor
	 */
	GenericNucleotideSequenceSlice<TYPE>( const GeneticSequence &hostObject, 
										  size_t uxStartOffsetInHostSection, // inclusive starting counting with 0
										  size_t uxEndOffsetInHostSection // exclusive, i.e. first element after the end of the sequence
										)
		 : GenericGeneticSequenceSlice<TYPE>( hostObject, uxStartOffsetInHostSection, uxEndOffsetInHostSection )
	{
		/* TO DO: Check, whether the hostObject is really some nucleotide sequence.
		 */
	} // constructor

	/* SECTION WITH METHODS FOR SEQUENCE EXTRACTIONS
	 */

	/* Copy the slice as genetic sequence using forward order.
	 */
	std::shared_ptr<NucleotideSequence> makeSequenceUsingForwardOrder() const
	{
		return std::make_shared<NucleotideSequence>( GenericGeneticSequenceSlice<TYPE>::hostObject, GenericGeneticSequenceSlice<TYPE>::xInterval );
	} // method

	/* Copy the slice as genetic sequence using reverse order.
	 */
	std::shared_ptr<NucleotideSequence> makeGenericSequenceUsingReverseOrder() const
	{
		auto pxSequenceRef = makeSequenceUsingForwardOrder();
		pxSequenceRef->vReverse();
		return pxSequenceRef;
	} // method

	/* The slice as sequence from the reverse strand (complemented bp) in forward order.
	 */
	std::shared_ptr<NucleotideSequence> makeComplementSequenceUsingForwardOrder() const
	{
		auto pxSequenceRef = std::make_shared<NucleotideSequence>( GenericGeneticSequenceSlice<TYPE>::hostObject, GenericGeneticSequenceSlice<TYPE>::xInterval );
		pxSequenceRef->vSwitchAllBasePairsToComplement();
		return pxSequenceRef;
	} // method

	/*  The slice as sequence from the reverse strand (complemented bp) in reverse order.
	 */
	std::shared_ptr<NucleotideSequence> makeComplementSequenceUsingReverseOrder() const
	{
		auto pxSequenceRef = makeGenericSequenceUsingReverseOrder();
		pxSequenceRef->vSwitchAllBasePairsToComplement();
		return pxSequenceRef;
	} // method
	
	/* SECTION WITH METHODS FOR STATISTIC EVALUATIONS
	 */
	 
	/* Helper method for counting single nucleotides.
	 */
	std::shared_ptr<std::string> NucleotideStatisiticAsString() const
	{
		/* Create and get the histogram
		 */
		auto pHistogramVectorRef = GenericGeneticSequenceSlice<TYPE>::pCreateHistogramVector();

		std::shared_ptr<std::string> pTextRef = std::make_shared<std::string>();

		for ( uint8_t uxIterator = 0; uxIterator < GenericGeneticSequenceSlice<TYPE>::hostObject.uxAlphabetSize(); uxIterator++ )
		{
			vAppendValueToString( *pTextRef, 
								  std::string ( 1, NucleotideSequence::translateACGTCodeToCharacter( uxIterator ) ), 
								  (*pHistogramVectorRef)[uxIterator]
								);
		} // for
		
		/* WARNING!!! If the table changes, we have to change this code as well.
		 * A == 0; C == 1; G == 2; T == 3;
		 */
		vAppendValueToString( *pTextRef, std::string("GC"), (*pHistogramVectorRef)[1] + (*pHistogramVectorRef)[2] );
		vAppendValueToString( *pTextRef, std::string("AT"), (*pHistogramVectorRef)[0] + (*pHistogramVectorRef)[3] );

		return pTextRef;
	} // method

	/* SECTION WITH METHODS FOR TEXTUAL OUTPUT
	 */
	
	/* Delivers a string object that contains a textual representation of the referred genetic sequence. 
	 */
	std::unique_ptr<std::string> asSequenceOfACGT() const
	{
		std::unique_ptr<std::string> psSequenceAsTextRef ( new std::string() );

		/* We resize the string (reserve enough memory) for the iteration coming.
		 */
		psSequenceAsTextRef->resize( this->xInterval.uxSize );

		for( size_t uxIterator = 0; uxIterator < this->xInterval.uxSize; uxIterator++ )
		{
			(*psSequenceAsTextRef)[uxIterator] = NucleotideSequence::translateACGTCodeToCharacter( GenericGeneticSequenceSlice<TYPE>::getSequenceRefInHost()[uxIterator] );
		} // for

		return psSequenceAsTextRef;
	} // method
}; // templated class

class NucleotideSequenceSlice : public GenericNucleotideSequenceSlice<NucleotideSequenceSlice> 
{
public :
	/* Forward of constructor
	 */
	NucleotideSequenceSlice( const GeneticSequence &hostObject,
						     size_t uxStartOffsetInHostSection,
					         size_t uxEndOffsetInHostSection
						   )
		: GenericNucleotideSequenceSlice<NucleotideSequenceSlice>( hostObject, uxStartOffsetInHostSection, uxEndOffsetInHostSection)
	{} // constructor
}; // class

/* Generic front-end of NucleotideSequenceSliceHavingSelectedStrand
 */
template<typename TYPE>
class GenericNucleotideSequenceSliceHavingSelectedStrand : public GenericNucleotideSequenceSlice<TYPE> // public GeneticSequenceSlice
{
public :
	typedef enum
	{
		FORWARD_STRAND = 0,
		REVERSE_STRAND = 1
	} tSelectedStrandType;

	/* The strand, where the slice shall belong to
	 */
	const tSelectedStrandType tSelectedStrand;

	/* Constructor
	 */
	GenericNucleotideSequenceSliceHavingSelectedStrand<TYPE>( const GeneticSequence &hostObject, // host sequence 
															  size_t uxStartOffsetInHostSection, // inclusive starting counting with 0
															  size_t uxEndOffsetInHostSection, // exclusive, i.e. first element after the end of the sequence
															  tSelectedStrandType tSelectedStrand // forward or reverse strand
															  )
		 : GenericNucleotideSequenceSlice<TYPE>( hostObject, uxStartOffsetInHostSection, uxEndOffsetInHostSection ),
		   tSelectedStrand( tSelectedStrand )
	{
		/* TO DO: Check, whether the hostObject is really some nucleotide sequence.
		 */
	} // constructor

	/* Forwards to 
	 *	   makeSequenceUsingForwardOrder()
	 *  or makeComplementSequenceUsingReverseOrder()
	 *  depending on the initially select strand. 
	 */
	std::shared_ptr<NucleotideSequence> makeSequence() const
	{
		return tSelectedStrand == FORWARD_STRAND ? GenericNucleotideSequenceSlice<TYPE>::makeSequenceUsingForwardOrder()
												 : GenericNucleotideSequenceSlice<TYPE>::makeComplementSequenceUsingReverseOrder();
	} // method

	/* Get a subsection of the current slice. Has to be overloaded here, because we have to copy the strand as well.
	 * The method is generic and adapts its return type.
	 */
	TYPE getSubSliceFromTo( size_t uxStartOffsetInSlice, size_t uxEndOffsetInSlice ) const
	{
		return TYPE( GenericGeneticSequenceSlice<TYPE>::hostObject, 
					 GenericGeneticSequenceSlice<TYPE>::xInterval.uxOffsetInHostSection + uxStartOffsetInSlice,
					 GenericGeneticSequenceSlice<TYPE>::xInterval.uxOffsetInHostSection + uxEndOffsetInSlice,
					 tSelectedStrand
				   );
	} // method

	/* Delivers a string object that contains a textual representation of the referred genetic sequence. 
	 */
	std::unique_ptr<std::string> asSequenceOfACGT() const
	{
		if ( tSelectedStrand == FORWARD_STRAND )
		{
			/* We forward this call simply to the generic nucleotide sequence 
			 */
			return GenericNucleotideSequenceSlice<TYPE>::asSequenceOfACGT();
		} // if
		
		BOOST_LOG_TRIVIAL(trace) << "asSequenceOfACGT for reverse strand";

		/* Code for the case of the reverse strand.
		 * Alternative approach: first create a revers strand nucleotide sequence and then print this sequence (drawback: intermediate object)
		 */
		std::unique_ptr<std::string> psSequenceAsTextRef ( new std::string() );

		/* We resize the string (reserve enough memory) for the iteration coming.
		 */
		psSequenceAsTextRef->resize( this->xInterval.uxSize );

		for( size_t uiIterator = 0; uiIterator < this->xInterval.uxSize; uiIterator++ )
		{
			/* We behave, so that we get the reverse strand sequence.
			 */
			size_t uiInversePos = this->xInterval.uxSize - uiIterator - 1;
			(*psSequenceAsTextRef)[uiIterator] = NucleotideSequence::translateACGTCodeToCharacter
				( NucleotideSequence::nucleotideComplement( GenericGeneticSequenceSlice<TYPE>::getSequenceRefInHost()[uiInversePos] ) );
		} // for

		return psSequenceAsTextRef;
	} // method
}; // templated class

/* For slices of nucleotide sequences, where we have a decision with respect to the selected strand.
 * Class is for efficiently working with double stranded sequences.
 */
class NucleotideSequenceSliceHavingSelectedStrand : public GenericNucleotideSequenceSliceHavingSelectedStrand<NucleotideSequenceSliceHavingSelectedStrand> 
{
public :
	/* Forward of constructor
	 */
	NucleotideSequenceSliceHavingSelectedStrand( const GeneticSequence &hostObject,
												 size_t uxStartOffsetInHostSection,
												 size_t uxEndOffsetInHostSection,
												 tSelectedStrandType tSelectedStrand // forward or reverse strand
											   )
		: GenericNucleotideSequenceSliceHavingSelectedStrand<NucleotideSequenceSliceHavingSelectedStrand>
		  ( hostObject, uxStartOffsetInHostSection, uxEndOffsetInHostSection, tSelectedStrand )
	{ } // constructor

}; // class

void vForAllTranslationsDo( int64_t iGeneId,
							int64_t iTaxonomicId,
							const std::function< void( const NucleotideSequence& ) > functor
							// Functor &&functor
							);

/* export this module to boost python */
void exportSequence();