// SwithWaterman.cpp : Defines the entry point for the console application.
//

#include <string>
#include <iostream>
#include <algorithm>
#include "sw_sse_avx.h"
#include "smith_waterman.h"


std::vector<ContainerType> SMW::getInputType()
{
	return std::vector<ContainerType>{
			//the query sequence
			ContainerType::nucSeq,
			//the ref
			ContainerType::packedNucSeq
		};
}//function


ContainerType SMW::getOutputType()
{
	return ContainerType::alignment;
}//function

/* Random nucleotide sequence of length uiLen
 */
std::string randomNucSeq(const int uiLen) {
	static const char nucleotides[] = "ACGT";
	
	std::string sNucSeq(uiLen, ' ');
	for (int i = 0; i < uiLen; ++i) 
	{
		sNucSeq[i] = nucleotides[rand() % (sizeof(nucleotides) - 1)];
	} // for
	
	return sNucSeq;
} // function


/* SIMD (SSE, AVX2) boosted Smith-Waterman alignments. 
 */
int16_t alignSW_SIMD( const NucleotideSequence &rQuerySequence, // query sequence
					  const NucleotideSequence &rReferenceSequence, // reference sequence
					  SmithWatermanParamaterSet<int16_t> &rSWparameterSet, // Smith Waterman alignment parameter
					  std::vector<size_t> &rvMaxScorePositions // vector will recieve positions, where we have a max score
					)
{
	int16_t iMaxScore = 0;

	try
	{	// Create aligner object on the foundation of the query
		SW_SIMD_Aligner<int16_t, size_t> xSIMDAligner( rQuerySequence, rSWparameterSet );
		
		/* Do the actual alignment using a reference.
		 * iMaxScore is the maximum score computed in the context of the alignment.
		 * rvMaxScorePositions is a vector that contains all positions having maxscore.
		 * Position counting starts with zero. (First row has index zero)
		 */
		iMaxScore = xSIMDAligner.align( rReferenceSequence, rvMaxScorePositions );
	} // try
	catch (std::exception &e)
	{	/* Aligner can throw expcetions of type AlignerException
		 */
		std::cout << "Catched exception: " << e.what();
		exit(0);
	} // catch

	return iMaxScore;
} // function


std::shared_ptr<Container> SMW::execute(
		std::vector<std::shared_ptr<Container>> vpInput
	)
{
	std::shared_ptr<NucleotideSequence> pQuerySeq = 
		std::static_pointer_cast<NucleotideSequence>(vpInput[0]);
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefPack = 
		std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[1]);


	// 1. Prepare the SW parameter set ...
	SmithWatermanParamaterSet<int16_t> xSWparameterSet
	(	1, // score for match (must be positive)
		-1, // score for mismatch (must be negative)
		20, // penalty for gap open ("gap" means insertion or deletion)
		1, // penalty for gap extension
		pQuerySeq->uxAlphabetSize() // alphabet size of input
	);
	
	// 2. Do the alignment ...
	std::vector<size_t> vMaxScorePositions;
	/*int16_t imaxScore =*/ alignSW_SIMD( *pQuerySeq, // query sequence
									  *pRefPack->vColletionAsNucleotideSequence(), // reference sequence
									  xSWparameterSet, // Smith Waterman alignment parameter
									  vMaxScorePositions // vector will recieve positions, where we have a max score
									);

	nucSeqIndex start = vMaxScorePositions[0];
	if(vMaxScorePositions.size() != 1)
		std::cout << "WARNING SMW found "<< vMaxScorePositions.size() <<" positions!" << std::endl;
	return std::shared_ptr<Alignment>(new Alignment(start, 0));
}//function


void demoCode()
{
	srand( (unsigned)time( NULL ) );	
	std::string sRandNucSeq = randomNucSeq( 10000000 );
	NucleotideSequence xReference( sRandNucSeq );

	size_t uiStart = 10000;
	size_t uiMatchSize = 1000;
	NucleotideSequence xQuery( sRandNucSeq.substr( uiStart, uiMatchSize ) );

	// 1. Prepare the SW parameter set ...
	SmithWatermanParamaterSet<int16_t> xSWparameterSet
	(	10, // score for match (must be positive)
		-3, // score for mismatch (must be negative)
		4, // penalty for gap open ("gap" means insertion or deletion)
		1, // penalty for gap extension
		xReference.uxAlphabetSize() // alphabet size of input
	);
	
	// 2. Do the alignment ...
	std::vector<size_t> vMaxScorePositions;
	int16_t imaxScore = alignSW_SIMD( xQuery, // query sequence
									  xReference, // reference sequence
									  xSWparameterSet, // Smith Waterman alignment parameter
									  vMaxScorePositions // vector will recieve positions, where we have a max score
									);

	// 3. Print info ...
	std::cout << "Max score is: " << imaxScore << "\n"
			  << "Max positions are: ";
	for( auto uiPosition : vMaxScorePositions )
	{
		std::cout << uiPosition << " ";
	} // for
	std::cout << std::endl;

	// 4. Selfchecks ...
	//size_t uiMatchEndIndex = uiStart + uiMatchSize - 1;
	if(    (  (size_t)imaxScore != xSWparameterSet.iWeightMatch * uiMatchSize )
		|| ( std::find(vMaxScorePositions.begin(), vMaxScorePositions.end(), imaxScore ) != vMaxScorePositions.end() ) )
	{
		std::cout << "Panic. Selfcheck failed 1!";
		exit(0);
	} // if

	for( auto uiPos : vMaxScorePositions )
	{	// Sequence at max score position must be equal to query
		if( sRandNucSeq.substr( uiPos - uiMatchSize + 1, uiMatchSize ) != sRandNucSeq.substr( uiStart, uiMatchSize ) )
		{
			std::cout << "Panic. Selfcheck 2 failed for position: " << uiPos << std::endl;
		} // if
	} // for 
} // function

void exportSMW()
{
	boost::python::def("demoSMW", &demoCode);
    boost::python::class_<
			SMW, 
			boost::python::bases<CppModule>,
        	std::shared_ptr<SMW>
			>(
			"SMW",
            "SMW alignment.\n"
        );
	boost::python::implicitly_convertible< 
		std::shared_ptr<SMW>,
		std::shared_ptr<CppModule> 
	>();
}