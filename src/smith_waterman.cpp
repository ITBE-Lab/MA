// SwithWaterman.cpp : Defines the entry point for the console application.
//

#include <string>
#include <iostream>
#include <algorithm>
#include "sw_sse_avx.h"
#include "smith_waterman.h"


ContainerVector SMW::getInputType() const
{
	return ContainerVector{
			//the query sequence
			std::shared_ptr<Container>(new NucSeq()),
			//the ref
			std::shared_ptr<Container>(new Pack())
		};
}//function


std::shared_ptr<Container> SMW::getOutputType() const
{
	return std::shared_ptr<Container>(
			new ContainerVector(std::shared_ptr<Container>(new Alignment()))
		);
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
int16_t alignSW_SIMD( const NucSeq &rQuerySequence, // query sequence
					  const NucSeq &rReferenceSequence, // reference sequence
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
		ContainerVector vpInput
	)
{
	std::shared_ptr<NucSeq> pQuerySeq = 
		std::static_pointer_cast<NucSeq>(vpInput[0]);
	std::shared_ptr<Pack> pRefPack = 
		std::static_pointer_cast<Pack>(vpInput[1]);


	// 1. Prepare the SW parameter set ...
	SmithWatermanParamaterSet<int16_t> xSWparameterSet
	(	1, // score for match (must be positive)
		-1, // score for mismatch (must be negative)
		20, // penalty for gap open ("gap" means insertion or deletion)
		1, // penalty for gap extension
		pQuerySeq->uxAlphabetSize() // alphabet size of input
	);

	std::shared_ptr<NucSeq> pReference = pRefPack->vColletionAsNucSeq();

	//reversing query and reference in order to obtain start instead of end (markus)
	pQuerySeq->vReverse();
	pReference->vReverse();
	
	// 2. Do the alignment ...
	std::vector<size_t> vMaxScorePositions;
	/*int16_t imaxScore =*/ alignSW_SIMD( *pQuerySeq, // query sequence
									  *pReference, // reference sequence
									  xSWparameterSet, // Smith Waterman alignment parameter
									  vMaxScorePositions // vector will recieve positions, where we have a max score
									);

	// 3. collect the results ...
	auto pvRet = std::make_shared<ContainerVector>();

	for(nucSeqIndex revStart : vMaxScorePositions)
	{
		//undo the reversion by substracting from absolute length
		std::shared_ptr<Container> pAlignment = std::shared_ptr<Alignment>(new Alignment(pRefPack->uiUnpackedSizeForwardPlusReverse() - revStart));
		pvRet->push_back(pAlignment);
	}//for

	return std::shared_ptr<ContainerVector>(new ContainerVector(pvRet));
}//function


void demoCode()
{
	srand( (unsigned)time( NULL ) );	
	std::string sRandNucSeq = randomNucSeq( 10000000 );
	NucSeq xReference( sRandNucSeq );

	size_t uiStart = 10000;
	size_t uiMatchSize = 1000;
	NucSeq xQuery( sRandNucSeq.substr( uiStart, uiMatchSize ) );

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