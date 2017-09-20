#include "sequence.h"


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
	boost::python::class_<
			NucleotideSequence, 
			boost::noncopyable,
			boost::python::bases<Container>, 
			std::shared_ptr<NucleotideSequence>
		>(
			"NucSeq", 
			"Holds a single nucleotide sequence.\n",
			boost::python::init<const char*>(
				"arg1: self\n"
				"arg2: string to initialize the sequence from\n"
			)
		)
			.def(boost::python::init<const std::string>(
				"arg1: self\n"
				"arg2: string to initialize the sequence from\n"
			))
			.def(
					"at", 
					&NucleotideSequence::charAt,
					"arg1: self\n"
					"arg2: index at which to look\n"
					"returns: the char at the given index\n"
				)
			.def(
					"__getitem__", 
					&NucleotideSequence::charAt,
					"arg1: self\n"
					"arg2: index at which to look\n"
					"returns: the char at the given index\n"
				)
			.def(
					"append", 
					&NucleotideSequence::vAppend_boost,
					"arg1: self\n"
					"arg2: sequence to append\n"
					"returns: nil\n"
					"\n"
					"Appends the given string to the end of the sequence.\n"
					"Any character other than A,C,T,G,a,c,t and g "
					"will result in an N beeing appended.\n"
				)
			.def(
					"length", 
					&NucleotideSequence::length,
					"arg1: self\n"
					"returns: the length of the sequence\n"
				)
			.def(
					"__len__", 
					&NucleotideSequence::length,
					"arg1: self\n"
					"returns: the length of the sequence\n"
				)
			.def(
					"__str__", 
					&NucleotideSequence::toString,
					"arg1: self\n"
					"returns: the sequence as string\n"
				)
			.def(
					"reverse", 
					&NucleotideSequence::vReverse,
					"arg1: self\n"
					"returns: nil\n"
					"\n"
					"Reverses the sequence.\n"
				);

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<NucleotideSequence>, 
			std::shared_ptr<Container>
		>(); 

	//register return values of vectors of nucseqs
	boost::python::class_<std::vector<std::shared_ptr<NucleotideSequence>>>("VecRetNuc")
	.def(boost::python::vector_indexing_suite<
			std::vector<std::shared_ptr<NucleotideSequence>>,
			/*
			*	true = noproxy this means that the content of the vector is already exposed by
			*	boost python. 
			*	if this is kept as false, StripOfConsideration would be exposed a second time.
			*	the two StripOfConsiderations would be different and not intercastable.
			*	=> keep this as true
			*/
			true
		>());

}//function