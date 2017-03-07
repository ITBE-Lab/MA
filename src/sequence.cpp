#include "sequence.h"

/* Delivers a slice that spans over the complete sequence.
 */
GeneticSequenceSlice GeneticSequence::fullSequenceAsSlice() const
{
	return GeneticSequenceSlice( *this, 0, this->uiSize );
} // method


/* Documentation - see class definition.
 */
NucleotideSequenceSlice NucleotideSequence::fullSequenceAsSlice() const
{
	return NucleotideSequenceSlice( *this, 0, this->uiSize );
} // method

NucleotideSequenceSliceHavingSelectedStrand NucleotideSequence::fullSequenceAsSliceHavingSelectedStrand
	( bool bSelectReverseStrand) 
const
{
	return NucleotideSequenceSliceHavingSelectedStrand
	( *this, 
	  0, 
	  this->uiSize, 
	  bSelectReverseStrand ? NucleotideSequenceSliceHavingSelectedStrand::REVERSE_STRAND 
						   : NucleotideSequenceSliceHavingSelectedStrand::FORWARD_STRAND
	);
} // method

/* Extracts a part of a genetic sequence as sequence slice.
 */
GeneticSequenceSlice GeneticSequence::fromTo( size_t uxFrom, size_t uxEnd )
{	/* Range check is done as part of the sequence slice constructor
	 */
	return GeneticSequenceSlice( *this, uxFrom, uxEnd );
} // method

GeneticSequenceSlice GeneticSequence::getSliceByInterval( const DeprecatedIntervalDescriptor &interval ) const
{
	return GeneticSequenceSlice( *this, interval.uxStart, interval.uxEnd ); 
} // method

NucleotideSequenceSlice NucleotideSequence::getSliceByInterval( const DeprecatedIntervalDescriptor &interval ) const
{
	return NucleotideSequenceSlice( *this, interval.uxStart, interval.uxEnd );
} // method

NucleotideSequenceSliceHavingSelectedStrand NucleotideSequence::getSliceByInterval( const DeprecatedDoubleStrandIntervalDescriptor &interval ) const
{
	return NucleotideSequenceSliceHavingSelectedStrand( *this, 
														interval.uxStart, 
														interval.uxEnd, 
														interval.eStrand == DeprecatedDoubleStrandIntervalDescriptor::IS_FORWARD_STRAND 
															? NucleotideSequenceSliceHavingSelectedStrand::FORWARD_STRAND
															: NucleotideSequenceSliceHavingSelectedStrand::REVERSE_STRAND
													  );
} // method

#if 0
/* Temporary bridge function
 */
NucleotideSequenceSlice GeneticSequenceSlice::castToNucleotideSequenceSlice()
{
	/* TO DO: Check whether the host sequence is really of type nucleotide sequence
	 */
	return NucleotideSequenceSlice( hostObject, xInterval.uxOffsetInHostSection, xInterval.uxOffsetInHostSection + xInterval.uxSize );
}
#endif

void GeneticSequence::vAppend( const GeneticSequenceSlice &rSequenceSlice )
{
	/* We redirect the call to an append of our base class.
	 */
	PlainSequence<uint8_t>::vAppend( rSequenceSlice.getSequenceRefInHost(), rSequenceSlice.getSize() );
} // method

void NucleotideSequence::vAppend( const NucleotideSequenceSlice &rSequenceSlice )
{
	/* We redirect the call to an append of our base class.
	 */
	PlainSequence<uint8_t>::vAppend( rSequenceSlice.getSequenceRefInHost(), rSequenceSlice.getSize() );
} // method

/* Performance hurting code.
 */
void NucleotideSequence::vAppend( const NucleotideSequenceSliceHavingSelectedStrand &rSequenceSlice )
{
	if ( rSequenceSlice.tSelectedStrand == NucleotideSequenceSliceHavingSelectedStrand::REVERSE_STRAND )
	{	/* TO DO: This solution was an emergency fix and is quite inefficient.
		 * Quickly solve this TO DO !!
		 */
		vAppend( (*( rSequenceSlice.makeSequence() )).fullSequenceAsSlice() );
	} // if
	else
	{	/* We redirect the call to an append of our base class.
		 * The true branch line would work as well, but is less inefficient.
		 */
		PlainSequence<uint8_t>::vAppend( rSequenceSlice.getSequenceRefInHost(), rSequenceSlice.getSize() );
	} // else
} // method

/* Catenate all intervals of the interval vector.
 * Write the outcome in the sequence given as second argument. 
 */
void NucleotideSequence::catenate( const DeprecatedIntervalVector<DeprecatedDoubleStrandIntervalDescriptor> &rxIntervalVector,
								   NucleotideSequence &rxNucleotideSequence
								 ) const
{
	rxIntervalVector.vForAllIntervalsDoCounting
	(
		[&]( const DeprecatedDoubleStrandIntervalDescriptor &interval, size_t uxCounter )
		{	/* Transform the interval to a slice and append.
			 */
			rxNucleotideSequence.vAppend( getSliceByInterval( interval ) );
		} // lambda
	); // function call
} // method

/* by markus
*	returns a sequence slice of the index from while checking if from >= 0 
*		with the length of uiSize while checking if the end of the slice is bigger than the sequence length
*	if the slice is out of bounds a slice of the size 0 is returned
*/
NucleotideSequenceSlice NucleotideSequence::saveFromSize(long uiFrom, size_t uiSize) const
{
	if (uiFrom < 0)
	{
		uiSize += uiFrom;
		uiFrom = 0;
	}//if
	if (uiFrom + uiSize > length())
		uiSize = length() - uiFrom;
	if (uiSize < 0)
		uiSize = 0;
	return fromTo(uiFrom , uiFrom + uiSize);
} // method

/* TO DO: Check the index for correct range.
 */
NucleotideSequenceSlice NucleotideSequence::fromTo( size_t uiFrom, size_t uiEnd ) const
{	
	return NucleotideSequenceSlice( *this, uiFrom, uiEnd );
} // method

#if 0
/* Transforms the genetic sequence into some set of pseudo-contigs.
 */
std::unique_ptr<std::vector<GeneticSequenceSlice>> GeneticSequence::generateVectorOfPseudoContigs
(
	unsigned int uiContigSize, // size of each contig
	unsigned int uiShift, // shifting of the start position between two contigs
	unsigned int uiAmplificationRate // umber of copies create of some single contig
)
{
	size_t uiContigStartIndex = 0;
	std::unique_ptr<std::vector<GeneticSequenceSlice>> pxReturnedVectorRef ( new std::vector<GeneticSequenceSlice>() );

	/* Iterate by continuously shifting.
	 */
	while ( uiContigStartIndex + uiContigSize + (uiContigSize / 2) < this->length() )
	{
		pxReturnedVectorRef->emplace_back( GeneticSequenceSlice( *this, uiContigStartIndex, uiContigStartIndex + uiContigSize ));
		
		uiContigStartIndex += uiShift;
	} // while

	pxReturnedVectorRef->emplace_back( GeneticSequenceSlice( *this, uiContigStartIndex, this->length() ) );

	return pxReturnedVectorRef;
} // method
#endif

#if 0
NucleotideSequenceSliceHavingSelectedStrand DoubleStrandedSequence::getSliceByInterval( DeprecatedDoubleStrandIntervalDescriptor &interval )
{
	NucleotideSequenceSliceHavingSelectedStrand( *pForwardStrandRef, // host sequence 
												 interval.uxStart, // inclusive starting counting with 0
												 interval.uxEnd, // exclusive, i.e. first element after the end of the sequence
												 interval.eStrand == DeprecatedDoubleStrandIntervalDescriptor::IS_FORWARD_STRAND
													 ? NucleotideSequenceSliceHavingSelectedStrand::FORWARD_STRAND
													 : NucleotideSequenceSliceHavingSelectedStrand::REVERSE_STRAND 
											   );
#if 0
	/* Deprecated code segment, DELETE LATER
	 */
	switch ( interval.eStrand )
	{
	case DeprecatedDoubleStrandIntervalDescriptor::IS_FORWARD_STRAND :
		return 
			
			GeneticSequenceSlice( *pForwardStrandRef, // we take the forward strand
									 interval.uxStart, // see PPT slide 
									 interval.uxEnd // see PPT slide
								   );
	
	case DeprecatedDoubleStrandIntervalDescriptor::IS_REVERSE_STRAND :
		return GeneticSequenceSlice( *pReverseStrandRef, // we take the reverse strand
									 pReverseStrandRef->uxGetSequenceSize() - interval.uxEnd, // see PPT slide 
									 pReverseStrandRef->uxGetSequenceSize() - interval.uxStart  // see PPT slide 
								   );

	default :	
		/* We have to report some misbehavior.
		 */
		throw AlignerException( "getSliceByInterval does not work for UNKNOWN strands" );
	}
#endif
} // method
#endif

/* The translation table for columns.
 * Translates a single character into a 2-bit compressed code.
 */
const unsigned char NucleotideSequence::xNucleotideTranslationTable[256] = 
{
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,  // A == 0; C == 1; G == 2;
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  // T == 3;
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,  // a == 0; c == 1; g == 2;
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  // t == 3;
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
}; // predefined array

/* The following array is used in order to translate a AA code into the corresponding symbol for output
 */
const char AminoAcidSequence::xCodeToSymbolTranslationArray[23] =
{ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '$' };

/* Translation table for translating a codon (NT triple) in to a Amino Acid code
 */
const char AminoAcidSequence::xNTpackedTripleToAACodeTranslationArray[64] =
{ 
  8  ,11 ,8  ,11 ,17 ,17 ,17 ,17,
  15 ,16 ,15 ,16 ,7  ,7  ,10 ,7 ,
  14 ,6  ,14 ,6  ,13 ,13 ,13 ,13,
  15 ,15 ,15 ,15 ,9  ,9  ,9  ,9 ,
  3  ,2  ,3  ,2  ,0  ,0  ,0  ,0 ,
  5  ,5  ,5  ,5  ,19 ,19 ,19 ,19,
  22 ,21 ,12 ,21 ,16 ,16 ,16 ,16,
  18 ,1  ,20 ,1  ,9  ,4  ,9  ,4
};  // array

/* AA atom counting table.
 * Organzied as follows:
 * bit 0 - 4 : C, bit 5 - 9 : H, bit 10 - 14 : N, bit 15 - 19 : O, bit 20 - 24 : S
 * The lines follow the AA order
 */
const uint32_t AtomLevelCount::aCHNOS_CountingTable[23] =
{	
	66787  , // A	Alanine
	1115363, // C	Cysteine63
	132324 , // D	Aspartic acid4
	132389 , // E	Glutamic acid9
	66921  , // F	Phenylalanine
	66722  , // G	Glycine
	68902  , // H	Histidine
	66982  , // I	Isoleucine
	68038  , // K	Lysine
	66982  , // L	Leucine
	1115493, // M	Methionine93
	100612 , // N	Asparagine2
	0	   , // O	Pyrrolysine
	66853  , // P	Proline
	100677 , // Q	Glutamine7
	70086  , // R	Arginine
	99555  , // S	Serine
	99620  , // T	Threonine
	0	   , // U	Selenocysteine
	66917  , // V	Valine
	67979  , // W	Tryptophan
	99689  , // Y	Tyrosine
	0	   	 // $	Stop codon
}; // counting Table

/* Appending a NT-sequence to an AA-sequence.
 * TO DO: Improve me, so that is is a real append.
 */
void AminoAcidSequence::vAppendNucleoditeSequence( const NucleotideSequence &xNucleotideSequence ) // throws Exception
{
	if ( xNucleotideSequence.uxGetSequenceSize() % 3 != 0 )
	{
		throw std::runtime_error( "NT to AA translation works only for sequences with size multiple of 3" );
	} // if

	/* Resize the current sequence appropriately.
	 */
	this->resize( xNucleotideSequence.uxGetSequenceSize() / 3 );

	size_t uiNTSequenceIndex = 0;
	for ( size_t uiAASeqeuenceIndex = 0; uiAASeqeuenceIndex < this->uiSize; )
	{	/* Collect 3 NT symbols and pack these in into uiPackedNTtriple for translation.
		 * TO DO: check additionally, that the symbol is in A, C, G, T.
		 */
		assert( uiNTSequenceIndex < xNucleotideSequence.uxGetSequenceSize() - 2 );
		
		uint8_t	uiPackedNTtriple =					  ( ((xNucleotideSequence[uiNTSequenceIndex++] & 0x03) << 4) & 0xFF );
				uiPackedNTtriple = uiPackedNTtriple | ( ((xNucleotideSequence[uiNTSequenceIndex++] & 0x03) << 2) & 0xFF );
				uiPackedNTtriple = uiPackedNTtriple | (  (xNucleotideSequence[uiNTSequenceIndex++] & 0x03	   ) & 0xFF );
		
		/* We write directly into the sequence for performance reasons.
		 */
		this->pxSequenceRef[uiAASeqeuenceIndex++] = xNTpackedTripleToAACodeTranslationArray[uiPackedNTtriple];
	} // for

	assert( uiNTSequenceIndex == xNucleotideSequence.uxGetSequenceSize() );
} // method

/* Counts the 5 atoms occurring in an AA-sequence
 */
AtomLevelCount AminoAcidSequence::vCountAtoms( void )
{
	AtomLevelCount xAtomLevelCount;
	for ( size_t uiAASeqeuenceIndex = 0; uiAASeqeuenceIndex < this->uiSize; uiAASeqeuenceIndex++ )
	{
		xAtomLevelCount.addAminoAcid( this->pxSequenceRef[uiAASeqeuenceIndex] );
	} // for
	return xAtomLevelCount;
} // method

/* Adds the atoms of the given amino-acid code to the counters
 * We could create a faster type of addition by relying on parallel additions.
 */
void AtomLevelCount::addAminoAcid( uint8_t uiAACode )
{
	assert( uiAACode < 23);

	uint32_t uiPackedCount = aCHNOS_CountingTable[uiAACode];
	for ( uint8_t uiAtomIndex = 0; uiAtomIndex < 5; uiAtomIndex++ )
	{	/* Add atoms for current atom to counter
		 */
		aAtomCounters[uiAtomIndex] += uiPackedCount & 0x1F;
		uiPackedCount = uiPackedCount >> 5;
	} // for
} // method

void exportSequence()
{
	 //export the nucleotidesequence class
	boost::python::class_<NucSeqContainer, boost::python::bases<Container>, std::shared_ptr<NucSeqContainer>>("NucSeq")
		.def(boost::python::init<const std::string>())
        .def("at", &NucSeqContainer::charAt)
        .def("append", &NucSeqContainer::vAppend);

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< std::shared_ptr<NucSeqContainer>, std::shared_ptr<Container> >(); 
}//function