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

void exportSequence()
{
	 //export the nucleotidesequence class
	boost::python::class_<NucSeqContainer, boost::python::bases<Container>, std::shared_ptr<NucSeqContainer>>("NucSeq")
		.def(boost::python::init<const std::string>())
		.def(boost::python::init<const char*>())
        .def("at", &NucSeqContainer::charAt)
        .def("append", &NucSeqContainer::vAppend)
        .def("size", &NucSeqContainer::size);

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< std::shared_ptr<NucSeqContainer>, std::shared_ptr<Container> >(); 
}//function