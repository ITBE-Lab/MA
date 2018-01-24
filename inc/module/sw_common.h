/**
 * @file sw_common.h
 * @brief Implements the smith waterman algorithm.
 * @author Arne Kutzner
 */

#pragma once
/* Common code used in Smith-Waterman aligners.
 */

#ifdef _WIN32   
typedef unsigned __int8 uint8_t;
#endif

/* As described in the papers:
 * TO DO: insert references
 */
template<class T_scoring, // type for scoring 
         class T_size_t      // size_t type
>
class SimilarityMatrix
{
public:
    const T_scoring xWeightMatch;
    const T_scoring xWeightMismatch;
    const uint8_t uxAlphabetSize;

private:
    /* We avoid the copy of Objects
    * In C++11: FastaString(const FastaString& ) = delete;
    */
    SimilarityMatrix(const SimilarityMatrix&);

    /* the scoring table must have size 5 X 5 of type T_scoring
    */
    void initialzeSimilarityMatrix()
    {
        /* initialize scoring matrix
        */
        T_size_t matrixIterator = 0;
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
    T_scoring* pSimilarityMatrixRef;

    SimilarityMatrix(T_scoring xWeightMatch, T_scoring xWeightMismatch, uint8_t uxAlphabetSize)
        : xWeightMatch(xWeightMatch),
        xWeightMismatch(xWeightMismatch),
        uxAlphabetSize(uxAlphabetSize)
    {
        /* We reserve memory for the scoring table and initialize the table.
        */
        pSimilarityMatrixRef = new T_scoring[uxAlphabetSize * uxAlphabetSize];
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
    T_scoring xCalculateMaxScore(T_size_t xSequenceLength)
    {
        return xSequenceLength * xWeightMatch;
    } // method
}; // class

template<class T_scoring>
class SmithWatermanParamaterSet
{
private:
    /* We avoid the copy of SmithWatermanParamaterSet Objects, for efficiency reasons.
    */
    SmithWatermanParamaterSet(const SmithWatermanParamaterSet&);

public:
    /* Weights for Match and Mismatch
    */
    const T_scoring iWeightMatch;
    const T_scoring iWeightMismatch;

    /* This size of the alphabet for out parameter set.
    */
    const unsigned int uiAlphabetSize;
    
    const SimilarityMatrix<T_scoring, size_t> similarityMatrix;

    const T_scoring iGapOpen;
    const T_scoring iGapExtend;


    SmithWatermanParamaterSet(T_scoring xWeightMatch, T_scoring xWeightMismatch, T_scoring xGapo, T_scoring xGape, uint8_t uiAlphabetSize)
        : iWeightMatch(xWeightMatch),
        iWeightMismatch(xWeightMismatch),
        uiAlphabetSize(uiAlphabetSize),
        similarityMatrix(xWeightMatch, xWeightMismatch, uiAlphabetSize),
        iGapOpen(xGapo),
        iGapExtend(xGape)
    { } // constructor

    T_scoring xGetApproximatedBackTrackDistance(size_t uxSequenceSize, T_scoring xScore)
    {
        /* WARNING 'size_t' to 'int16_t', possible loss of data
        */
        T_scoring xMaximallyPossibleScore = (T_scoring)(uxSequenceSize * iWeightMatch);

        T_scoring xScoreDifference = xMaximallyPossibleScore - xScore;

        return (T_scoring)uxSequenceSize + (xScoreDifference / iGapExtend);
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
