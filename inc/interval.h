#pragma once

#include <functional>
#include <iostream>
#include "exception.h"

class GeneticSequence;
class NucleotideSequence;
class GeneticSequenceSlice;
class NucleotideSequenceSlice;
class NucleotideSequenceSliceHavingSelectedStrand;

/* NCBI name consisting of a gi and a ref name.
 * FIX ME: Move this into some modul NCBI_general
 */
class NCBI_GiRefName
{
private :
	/* gi number and name of the ref-name
	 */
	std::string _sRefName;
	uint64_t _iGiNumber = 0;

public :
	friend class NCBIRefNamedOneBasedIntervalDescriptor;

	inline const std::string &sRefName()
	{
		return _sRefName;
	}

	NCBI_GiRefName( const std::string &sRefName, uint64_t iGiNumber = 0 ) 
	:	 _sRefName( sRefName ),
		 _iGiNumber( iGiNumber )
	{} // constructor

	/* Default constructor
	 */
	NCBI_GiRefName()
	{} // constructor

	NCBI_GiRefName( const std::string &rsGiRefName )
	{
		vParseUsingString( rsGiRefName );
	} // constructor

	/* We check, whether the names are equal.
	 */
	bool operator==( const NCBI_GiRefName &other ) const
	{
		return this->_sRefName == other._sRefName;
	} // operator
	
	/* Decomposes a named OneBasedIntervalDescriptor following the scheme:
	 * 1) text "gi"
	 * 2) gi number
	 * 3) text "ref"
	 * 4) ref name
	 * where all 4 fields occur separated by some bar and all is closed by some additional final bar.
	 * Example: gi|568815597|ref|NC_000001.11|
	 */
	void vParseUsingString( const std::string &rsGiRefName )
	{
		std::string sTextGi, sTextRef;

		/* Parse the string and set 4 variables according to the outcome.
		 */
		vSplitAndSet<false>		// false => everything after fourth token is ignored ...
		(	rsGiRefName,		// (in) input string
		 	'|',				// (in) separator, delimiter 
			sTextGi,			// (out) receives first token by reference
			this->_iGiNumber,	// (out) receives second token by reference
			sTextRef,			// (out) receives third token by reference
			this->_sRefName		// (out) receives fourth token by reference
		); // function call

		/* Check whether there is really a NCBI ref on the correct position.
		 */
		if ( sTextRef != "ref" )
		{
			throw DataReceiverException( "Parsed name has not the form gi|<number>|ref|<text>" );
		} // if
	} // method
}; // class

typedef enum
{
	FORWARD_STRAND = 0,
	REVERSE_STRAND = 1
} tSelectedStrandType;

/* Class for describing intervals within sequences using 0-based scheme.
 *  - counting start with 0
 *  - half-closed-half-open interval
 */
class ZeroBasedIntervalDescriptor
{
public :
	/* The offset from the start of the Host section.
	 */
	const size_t uxOffsetInHostSection;

	/* The size of the slice.
	 */
	const size_t uxSize;

	/* Preferred constructor setting offset and size.
	 */
	ZeroBasedIntervalDescriptor( size_t uxOffsetInHostSection, size_t uxSize ) : 
		uxOffsetInHostSection( uxOffsetInHostSection ), 
		uxSize( uxSize )
	{ } // constructor
}; // class

/* Class for describing intervals within sequences using 1-based scheme.
 *  - counting start with 1
 *  - region is specifieded by a closed interval
 */
class NCBIRefNamedOneBasedIntervalDescriptor
{
private :
	size_t _uxStartOffset = 0;

	size_t _uxEndOffset = 0;

	tSelectedStrandType _eSelectedStrand = FORWARD_STRAND;

	NCBI_GiRefName _sNCBIRefName;

public :
	inline const NCBI_GiRefName &xNCBIRefName()
	{
		return _sNCBIRefName;
	} // accessor

	inline const size_t &uxStartOffset()
	{
		return _uxStartOffset;
	} // accessor

	/* The default constructor keeps the primitive values uninitialized!
	 */
	NCBIRefNamedOneBasedIntervalDescriptor( const std::string &rsQNAMEAsString )
	{
		vParseFrom5Token( rsQNAMEAsString );
	} // constructor

	NCBIRefNamedOneBasedIntervalDescriptor()
	{} // default constructor

	/* Decomposes a named OneBasedIntervalDescriptor following the scheme:
	 * 1) reference sequence the contig originates from
	 * 2) start on reference sequence
	 * 3) end on reference sequence
	 * 4) the strand 0:forward strand, 1:backward strand
	 * 5) the interval number (if belonging to Exon or Intron)
	 * where all 5 fields occur separated by some bar.
	 */
	void vParseFrom5Token( const std::string &rsQNAMEAsString )
	{
		int iFowardReverse, iIntervalNumber;
		std::string sGeneName;
	
		/* Parse the string and set the 5 references according to the outcome.
		 */
		vSplitAndSet<false>					// false => everything after iIntervalNumber is ignored ...
		(	rsQNAMEAsString,				// input string
		 	'|',							// separator, delimiter 
			sGeneName,
			this->_sNCBIRefName._sRefName,	// receives first token by reference (for efficiency direct access to the private attribute)
			this->_uxStartOffset,			// receives second token by reference
			this->_uxEndOffset,				// receives third token by reference
			iFowardReverse,					// receives fourth token by reference
			iIntervalNumber					// receives sixth token by reference
		); // function call
		this->_eSelectedStrand = iFowardReverse == 0 ? FORWARD_STRAND : REVERSE_STRAND;
		this->_sNCBIRefName._iGiNumber = 0;
	} // method
}; // class

/************************************************************************/
/* DEPRECATED STUFF                                                     */
/************************************************************************/

/* DistanceDescriptor is for keeping data of ans start-end form
 * Counting start with 0, uxStart is the index of the first element, uxEnd the index of the first element after the end of the sequence.
 */
class DeprecatedIntervalDescriptor
{
public :
	const size_t uxStart;
	const size_t uxEnd;
	std::string sAccession; // The name of the interval (should be part of some sequence slice)

	DeprecatedIntervalDescriptor( size_t uxStart, size_t uxEnd ) : uxStart(uxStart), uxEnd(uxEnd)
	{ } // constructor

	void vStringAppendTextualDescription( std::string &string) const;

	void vStringAppendTextualDescription( std::string &string, const NucleotideSequenceSlice *pSequenceSliceRef ) const;

	void vStringAppendTextualDescription( std::string &string, const NucleotideSequenceSliceHavingSelectedStrand *pSequenceSliceRef ) const;

	/* the length of our interval
	 */
	size_t length()
	{
		return uxEnd - uxStart;
	} // method
}; // class

/* Intervals on chromosomes are additionally on the forward or reverse strand
 * This class keeps this additional information.
 */
class DeprecatedDoubleStrandIntervalDescriptor : public DeprecatedIntervalDescriptor
{
private :
	static const char *_namesForStrandTypes[3];

public :
	typedef enum 
	{
		UNKNOWN_STRAND		= 0,
		IS_FORWARD_STRAND	= 1,
		IS_REVERSE_STRAND	= 2
	} tStrandType;

public :
	/* Intervals for double stranded sequences keep the additional info whether we look on the forward strand
	 * or on the reverse strand.
	 * TO DO: make me eStrand const.
	 */
	tStrandType eStrand;

	DeprecatedDoubleStrandIntervalDescriptor( size_t uxStart, 
									size_t uxEnd, 
									tStrandType eStrand = UNKNOWN_STRAND 
								  ) 
		: DeprecatedIntervalDescriptor( uxStart, uxEnd ),
		  eStrand( eStrand )
	{ } // constructor
	
	inline void vStringAppendTextualDescription( std::string &string, const NucleotideSequenceSlice *pSequenceRef ) const
	{
		string.append(" strand : ").append( pcTextForStrandType( eStrand ) ).append(" ");
		DeprecatedIntervalDescriptor::vStringAppendTextualDescription( string, pSequenceRef );
	} // method

	inline void vStringAppendTextualDescription( std::string &string, const NucleotideSequenceSliceHavingSelectedStrand *pSequenceRef ) const
	{
		string.append( " strand : " ).append( pcTextForStrandType( eStrand ) ).append( " " );
		DeprecatedIntervalDescriptor::vStringAppendTextualDescription( string, pSequenceRef );
	} // method

	static const char *pcTextForStrandType( const tStrandType eStrand )
	{
		return _namesForStrandTypes[ eStrand <= IS_REVERSE_STRAND ? eStrand : 0 ];
	}; // static function
};

/* Abstract class for general Distance related data collections.
 * Base class should be either IntervalDescriptor for single stranded sequences 
 * or DoubleStrandIntervalDescriptor for double stranded sequences
 * TO DO: Change the design, so that vector becomes a private member!
 */
template <typename BASECLASS>
class DeprecatedIntervalVector : public std::vector<BASECLASS> 
{
private :
	/* IntervalDescriptorVector(const IntervalDescriptorVector&);
	 * Override this in a subclass, so that it delivers the appropriate value.
	 * This should give a textual description of the nature of data hold in the vector.
	 */
	virtual std::string sName()
	{
		return std::string( "Intervals" );
	} // method

public :
	template <typename TP_FUNC_APPLY>
	void vForAllIntervalsDoCounting( TP_FUNC_APPLY &&fApply ) const
	{
		unsigned int uxCounter = 1;
		for ( const BASECLASS &descriptor : *this )
		{
			fApply( descriptor, uxCounter );
			uxCounter++;
		} // for
#if 0
		std::for_each
			(	std::vector<BASECLASS>::begin(), 
				std::vector<BASECLASS>::end(), 
				[&uxCounter] ( BASECLASS descriptor ) 
					{ 
						fApply( descriptor, uxCounter ); 
						uxCounter++;
					} // lambda
			); // for_each
#endif
	} // generic method

	/* TO DO: Check whether there is a solution in C++ without the (unnecessary) lambda term
	 */
	template <typename TP_FUNC_APPLY>
	void vForAllIntervalsDo( TP_FUNC_APPLY &&fApply )
	{
		std::for_each( std::vector<BASECLASS>::begin(), std::vector<BASECLASS>::end(), fApply ); 
	} // generic method

	/* Calculates the statistical average (mean value) of all interval sizes
	 */
	size_t meanValueOfAllIntervalSizes()
	{
		size_t uxSizeOfAllIntervals = 0;

		vForAllIntervalsDo
		( 
			[&]( BASECLASS &interval )
			{
				uxSizeOfAllIntervals += interval.length();
			} // lambda
		); // function call

		return this->size() > 0 ? ( uxSizeOfAllIntervals / this->size() ) : 0;
	} // method

#if 0
	template <typename SEQUENCE_TP>
	void vWriteAllIntervalsAsFasta( const char *pcAccessionPrefix, // for the fasta reference sequence name
									const std::string &sGeneId, 
									std::ostream &rxOutputStream, // -> OUT the stream in which we write our output 
									CppSQLiteExtInsertStatementForTableInputs &xSQLInsertRowIntoTableInputs, // -> OUT The insert statement used for appending one row to table
									SEQUENCE_TP* pSequenceRef = NULL // the root sequence, from which we take the introns and exons
								  )
	{
		/* Strange GCC 4.7 bug: The replacement of sExonOrIntron in the ostream-statement leads to a segmentation fault.
		 * Only in GCC not MSVC.
		 */
		const std::string sExonOrIntron = sName();

		vForAllIntervalsDoCounting
		( 
			[&]( BASECLASS &interval, size_t uxCounter )
			{
				/* Counter for numbering.
				 * sName will resolve to Exon or Intron.
				 */
				const std::string sSequenceName( std::string( sGeneId ).append( "_" ).append( pcAccessionPrefix ).append( "_" ).append( sExonOrIntron ).append( std::to_string(uxCounter ) ) );
				
				rxOutputStream << ">" << sSequenceName << "\n";
#if 1
				/* The basic data of the current interval as, e.g. start, end, length
			     */
				if ( pSequenceRef != nullptr )
				{
					auto xSequenceSlice = pSequenceRef->getSliceByInterval( interval );
					rxOutputStream << *( xSequenceSlice.asString() );
				} // if
				rxOutputStream << "\n";
#endif			
				xSQLInsertRowIntoTableInputs( std::stol( sGeneId ),		// int,			 1. e.g. (729533 for FAM72A)
											  pcAccessionPrefix, 		// const char *, 2. e.g. NM?NP?
											  sExonOrIntron != "Intron"	// int,			 3. 1 -> Exon / 0 -> Intron
											  uxCounter,				// int,			 4. Number of the exon intro [1..n]
											  sSequenceName.c_str()		// const char *	 5. Name of the query sequence (fasta identification)
											);				
			} // lambda
		); // function call
	} // method
#endif

	/* For getting a textual representation of our sequence data
	 * If we get a pSequenceRef we cut the sections of this sequence and decompose it.
	 * SEQUENCE_TP can be single stranded as well as double stranded.
	 */
	// template <typename SEQUENCE_TP>
	void vStringAppendTextualDescription( std::string &sString, NucleotideSequence* pSequenceRef = NULL )
	{
		/* The name of our interval group. (e.g. exon)
		 */
		sString.append( "\n" )
			   .append( sName() )
			   .append( " (# ")
			   .append( numberToStringDeprecated( this->size() ) )
			   .append( ") average interval size: " )
			   .append( numberToStringDeprecated( meanValueOfAllIntervalSizes() ) )
			   .append( "\n" );

		vForAllIntervalsDoCounting
		( 
			[&]( const BASECLASS &interval, size_t uxCounter )
			{
				/* Counter for numbering.
				 */
				sString.append("(").append( numberToStringDeprecated(uxCounter) ).append( ") " );

				/* The basic data of the current interval as, e.g. start, end, length
			     */
				if ( pSequenceRef == NULL )
				{
					interval.vStringAppendTextualDescription( sString, (NucleotideSequenceSlice *)nullptr );
				} // if
				else
				{
					sString.append("FIX ME here has to come the interval as textual description");
					// auto xSequenceSlice = pSequenceRef->getSliceByInterval( interval ); 
					// interval.vStringAppendTextualDescription( sString, &xSequenceSlice );
				} // else
				sString.append( "\n" );
				
				sString.append( "\n" );
			} // lambda
		); // function call
	} // method

	/* Debug function for dumping an interval vector
	 */
	void vDumpToCout()
	{
		std::string string;
		vStringAppendTextualDescription( string );
		std::cout << string << std::endl;
	} // method
}; // class

/* Vector class for keeping intron related interval information
 */
class IntronVector : public DeprecatedIntervalVector<DeprecatedDoubleStrandIntervalDescriptor>
{
public :
	virtual std::string sName( )
	{
		return std::string( "Intron" );
	} // method
}; // class

/* Vector class for keeping exon related interval information
 * There are double-stranded and single-stranded intervals.
 */
template <typename BASECLASS>
class ExonVector : public DeprecatedIntervalVector<BASECLASS>
{
public :
	virtual std::string sName( )
	{
		return std::string( "Exon" );
	} // method
}; // class