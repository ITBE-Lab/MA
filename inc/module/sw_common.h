/**
 * @file sw_common.h
 * @brief Implements the smith waterman algorithm.
 * @author Arne Kutzner
 */

#pragma once
/* Common code used in Smith-Waterman aligners.
 */

/* As described in the papers:
 * TO DO: insert references
 */
template<class SCORE_TP>
class SimilarityMatrix
{
public:
	const SCORE_TP xWeightMatch;
	const SCORE_TP xWeightMismatch;
	const uint8_t uxAlphabetSize;

private:
	/* We avoid the copy of Objects
	* In C++11: FastaString(const FastaString& ) = delete;
	*/
	SimilarityMatrix(const SimilarityMatrix&);

	/* the scoring table must have size 5 X 5 of type SCORE_TP
	*/
	void initialzeSimilarityMatrix()
	{
		/* initialize scoring matrix
		*/
		size_t matrixIterator = 0;
		for (uint8_t i = 0; i < uxAlphabetSize - 1; ++i)
		{
			for (uint8_t j = 0; j < uxAlphabetSize - 1; ++j)
				pSimilarityMatrixRef[matrixIterator++] = i == j ? this->xWeightMatch
				: this->xWeightMismatch;

			pSimilarityMatrixRef[matrixIterator++] = 0; // ambiguous base
		}
		/* Fill the final row of the scoring matrix with zeros.
		*/
		for (int j = 0; j < uxAlphabetSize; ++j)
		{
			pSimilarityMatrixRef[matrixIterator++] = 0;
		}
	} // method

public:
	/* Quick lookup table for scoring. \
	* At the moment we allow direct access because of performance.
	* Here we should integrate a friend declaration.
	*/
	SCORE_TP* pSimilarityMatrixRef;

	SimilarityMatrix(SCORE_TP xWeightMatch, SCORE_TP xWeightMismatch, uint8_t uxAlphabetSize)
		: xWeightMatch(xWeightMatch),
		xWeightMismatch(xWeightMismatch),
		uxAlphabetSize(uxAlphabetSize)
	{
		/* We reserve memory for the scoring table and initialize the table.
		*/
		pSimilarityMatrixRef = new SCORE_TP[uxAlphabetSize * uxAlphabetSize];
		initialzeSimilarityMatrix();
	} // constructor

	  /* The destructor frees the allocated memory automatically
	  */
	~SimilarityMatrix()
	{
		delete[] pSimilarityMatrixRef;
	} // destructor

	  /* Gets the maximum of the similarity matrix and multiplies it with the sequence size
	  */
	SCORE_TP xCalculateMaxScore( size_t xSequenceLength)
	{
		return xSequenceLength * xWeightMatch;
	} // method
}; // class

template<class SCORE_TP>
class SmithWatermanParamaterSet
{
private:
	/* We avoid the copy of SmithWatermanParamaterSet Objects, for efficiency reasons.
	*/
	SmithWatermanParamaterSet(const SmithWatermanParamaterSet&);

public:
	/* Weights for Match and Mismatch
	*/
	const SCORE_TP iWeightMatch;
	const SCORE_TP iWeightMismatch;
	/* This size of the alphabet for out parameter set.
	*/
	const unsigned int uiAlphabetSize;

	const SimilarityMatrix<SCORE_TP> similarityMatrix;

	const SCORE_TP iGapOpen;
	const SCORE_TP iGapExtend;


	SmithWatermanParamaterSet(SCORE_TP xWeightMatch, SCORE_TP xWeightMismatch, SCORE_TP xGapo, SCORE_TP xGape, uint8_t uiAlphabetSize)
		: iWeightMatch(xWeightMatch),
		iWeightMismatch(xWeightMismatch),
		uiAlphabetSize(uiAlphabetSize),
		similarityMatrix(xWeightMatch, xWeightMismatch, uiAlphabetSize),
		iGapOpen(xGapo),
		iGapExtend(xGape)
	{ } // constructor

	SCORE_TP xGetApproximatedBackTrackDistance(size_t uxSequenceSize, SCORE_TP xScore)
	{
		/* WARNING 'size_t' to 'int16_t', possible loss of data
		*/
		SCORE_TP xMaximallyPossibleScore = (SCORE_TP)(uxSequenceSize * iWeightMatch);

		SCORE_TP xScoreDifference = xMaximallyPossibleScore - xScore;

		return (SCORE_TP)uxSequenceSize + (xScoreDifference / iGapExtend);
	} // method

	/* Delivers the maximal possible score for the given query len.
	 */
	long iPossibleMaxScore(size_t uiQueryLen) const
	{
		return (long)(iWeightMatch * uiQueryLen);
	} // public method

	  /* Delivers a textual description for all parameter of the parameter set
	  */
	std::string getTextualDescription()
	{
		std::string sDescripton = "SW-";
		sDescripton.append(std::to_string(iWeightMatch)).append("-")
			.append(std::to_string(iWeightMismatch)).append("-")
			.append(std::to_string(iGapOpen)).append("-")
			.append(std::to_string(iGapExtend));

		return sDescripton;
	} // method
}; // class
