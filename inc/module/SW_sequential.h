#pragma once
/* Naive implementation of Smith-Waterman for debugging purposes */
#include "container/nucSeq.h"
#include "module/sw_common.h"
#include <vector>
#include <functional>

typedef enum {
	STOP	= 0,
	LEFT_UP = 1,
	UP		= 2,
	LEFT	= 3
} sw_direction_t;

typedef enum {
	EQUAL_PAIR				 = 0,
	INSERTION_AT_ROW_SIDE	 = 1,
	INSERTION_AT_COLUMN_SIDE = 2,
	UNEQUAL_PAIR			 = 3
} e_alignment_element_kind_t;

typedef enum 
{
	ROW_SIDE	= 0,
	COLUMN_SIDE	= 1,
	BAR			= 2
} eSelectionType;

template<class SCORE_TP>
struct generic_eh_type {
	SCORE_TP h, e;
};

template<class SCORE_TP>
struct sw_alignment_result_t
{
	/* The absolute index.
	* It has to interpreted as if the matrix were a plain data structure
	*/
	size_t uxIndex;
	SCORE_TP score;
};

/* The compare function used in the below qsort
*/
template<class SCORE_TP>
int compare_sw_alignment_result_t(const void* p1, const void* p2) 
{
	const sw_alignment_result_t<SCORE_TP> f1 = *(const sw_alignment_result_t<SCORE_TP>*) p1;
	const sw_alignment_result_t<SCORE_TP> f2 = *(const sw_alignment_result_t<SCORE_TP>*) p2;
	return f2.score - f1.score;
} // compare_sw_alignment_result_t


  /* class for the description of alignments
  */
template<class T>
class alignment_description_element
{
private:
	char cMarkedElement( const char c ) const
	{
		return (eElementKind == UNEQUAL_PAIR) ? tolower(c) : c;
	} // private method

public :
	/* Data elements of single alignment.
	*/
	e_alignment_element_kind_t eElementKind;
	T entryOnRowSide;
	T entryOnColumnSide;

	alignment_description_element( e_alignment_element_kind_t eElementKind, 
								   T entryOnRowSide,
								   T entryOnColumnSide ) : eElementKind( eElementKind ), 
		entryOnRowSide( entryOnRowSide ), 
		entryOnColumnSide( entryOnColumnSide )
	{ } // constructor
}; // class

template<class T>
class alignment_description : public std::vector< alignment_description_element<T> >
{
private:
	/* Explanation for the typename below: (MSVC2013 accepts the code without typename)
	* http://www.parashift.com/c++-faq-lite/nondependent-name-lookup-types.html
	*/
	typedef typename std::vector< alignment_description_element<T> >::iterator vectorIterator;

	void appendSingleRow( typename std::vector< alignment_description_element<T> >::iterator begin, 
						  typename std::vector< alignment_description_element<T> >::iterator end, 
						  std::string &text, 
						  eSelectionType selectedSide )
	{
		std::for_each( begin, end, [&]( alignment_description_element<T> &element )
		{
			text.append( element.sHtmlText( selectedSide ) );
		} // lambda
		); // for_each
	} // method

public:
	/* Increments rForwardIterator as long as the predicate is true
	*/
	vectorIterator vIterateForwardWhile ( vectorIterator forwardIterator, 
										  const std::function<bool ( alignment_description_element<T> )>& predciateFunction )
	{
		vectorIterator resultIterator = forwardIterator;

		while (   ( resultIterator < this->end() ) 
				&& ( predciateFunction( *resultIterator ) == true )
				)
		{
			resultIterator++;
		};

		return resultIterator;
	} // method


	std::shared_ptr<std::stringstream> createTextForPositionInformation(   
		vectorIterator &rSectionStartIterator,
		vectorIterator &rSectionEndIterator,
		const size_t uxStartForAbsolutePosition,
		const size_t uxReverseStrandSize )
	{
		std::shared_ptr<std::stringstream> pStringStreamRef = std::make_shared<std::stringstream>();
		if ( uxReverseStrandSize == 0 )
		{

			( *pStringStreamRef ) << ( uxStartForAbsolutePosition + rSectionStartIterator ) - this->begin() + 1 
				<< " - " 
				<< ( uxStartForAbsolutePosition + rSectionEndIterator ) - this->begin();
		} // if
		else
		{
			( *pStringStreamRef ) << uxReverseStrandSize - ( uxStartForAbsolutePosition + rSectionEndIterator - this->begin() ) + 1 
				<< " - " 
				<< uxReverseStrandSize - ( uxStartForAbsolutePosition + rSectionStartIterator - this->begin() );
		} // else

		return pStringStreamRef;
	} // method


	std::string & rsAppendStringRepesentation( std::string &string ) 
	{
		int stride = 100;
		vectorIterator frontIterator = this->begin();

		while ( frontIterator < this->end() ) 
		{
			vectorIterator tailIterator = ( this->end() - frontIterator ) >= stride ? frontIterator + stride
				: this->end();
			appendSingleRow( frontIterator, tailIterator, string, COLUMN_SIDE );
			string.append( "\n" );
			appendSingleRow( frontIterator, tailIterator, string, BAR );
			string.append( "\n" );
			appendSingleRow( frontIterator, tailIterator, string, ROW_SIDE );
			string.append( "\n" );

			frontIterator = tailIterator;
		} // while

		return string;
	} // method
}; // class


/* Class for the description of scoring matrices.
 * Scoring matrices can be used by the serial aligner for backtracking purposes.
 */
template<class SCORE_TP, class T_size_t>
struct AlignmentOutcomeMatrix
{
	SCORE_TP* scoringOutcomeMatrix;
	sw_direction_t*	backtrackMatrix;

#if 0
	SequenceString &rRowSequenceString;
	SequenceString &rColumnSequenceString;
#endif

	const T_size_t numberOfColumns;
	const uint8_t* puxColumnSequenceRef;
	const T_size_t numberOfRows;
	const uint8_t* puxRowSequenceRef;


	/* Initializes the first column and the first row.
	*/
	void initializeFristColumnAndFirstRow() 
	{

		for(T_size_t uxIterator = 0; uxIterator < numberOfColumns; uxIterator++ )
		{
			scoringOutcomeMatrix[ uxIterator ] = 0;
			backtrackMatrix[ uxIterator ] = STOP;
		} // for

		for(T_size_t uxIterator = 0; uxIterator < numberOfRows; uxIterator++ )
		{
			scoringOutcomeMatrix[ uxIterator * numberOfColumns ] = 0;
			backtrackMatrix[ uxIterator * numberOfColumns ] = STOP;
		} // for
	}  // method

	   /* Construction of an alignment outcome matrix.
	   */
	AlignmentOutcomeMatrix( T_size_t numberOfColumns, const uint8_t*columnSequence, 
							T_size_t numberOfRows, const uint8_t*rowSequence ) 
		: numberOfColumns( numberOfColumns + 1 ), // + 1 because the matrix gets a virtual first column 
		puxColumnSequenceRef( columnSequence ),
		numberOfRows( numberOfRows + 1 ), // + 1 because the matrix gets a virtual first row
		puxRowSequenceRef( rowSequence )
	{
		const T_size_t sizesOfOutcomeAndBacktrackingMatrix = (this->numberOfColumns) * (this->numberOfRows);

		/* Reserve memory for the matrices with scoring and backtracking information
		*/
		scoringOutcomeMatrix = new SCORE_TP[sizesOfOutcomeAndBacktrackingMatrix];
		backtrackMatrix = new sw_direction_t[sizesOfOutcomeAndBacktrackingMatrix];

		initializeFristColumnAndFirstRow();
	} // constructor

	  /* Destructor.
	  */
	~AlignmentOutcomeMatrix()
	{
		delete[] scoringOutcomeMatrix;
		delete[] backtrackMatrix;
	} // destructor

	  /* Dumps the matrix to cout for debugging purposes
	  */
	void dump()
	{
		char directionSymbols[4] = {'x', '\\', '^', '<'};

		std::cout << "**********************************************" << std::endl;
		std::cout << "The scoring matrix is given by  " << std::endl << std::endl;

		std::cout << "    - \t";

		for( size_t uxIteratorColum = 1; uxIteratorColum < numberOfColumns; uxIteratorColum++ )
		{
			std::cout << NucSeq::translateACGTCodeToCharacter(puxColumnSequenceRef[ uxIteratorColum - 1 ]) << "\t";
		}

		std::cout << std::endl;

		for( size_t uxIteratorRow = 0; uxIteratorRow < numberOfRows; uxIteratorRow++ )
		{
			if (uxIteratorRow > 0)
				std::cout << NucSeq::translateACGTCodeToCharacter(puxRowSequenceRef [uxIteratorRow - 1]) << " : ";
			else
				std::cout << "- : ";

			for( size_t uxIteratorColum = 0; uxIteratorColum < numberOfColumns; uxIteratorColum++ )
			{
				std::cout << directionSymbols[backtrackMatrix[uxIteratorRow * numberOfColumns + uxIteratorColum]]
					<< scoringOutcomeMatrix[(uxIteratorRow * numberOfColumns + uxIteratorColum)]
					<< " "; // << "\t";
			} // for uxIteratorColum
			std::cout << std::endl;
		} // for uxIteratorRow
	} // method

	  /* Performs a backtrack within the scoring matrix.
	  * The caller has to deliver enough space in resultingSequenceForRow as well as resultingSequenceForColumn.
	  */
	void backtrackFromIndexText( size_t startIndex,
								 char *resultingSequenceForRow,
								 char *resultingSequenceForColumn,
								 size_t &ruiNumberMatches,
								 size_t &ruiNumberMismatches,
								 size_t &ruiNumberInsertions,
								 size_t &ruiNumberDeletions
	)
	{
		size_t index = startIndex;
		size_t endIndex = startIndex;
		ruiNumberMatches = ruiNumberMismatches = ruiNumberInsertions = ruiNumberDeletions = 0;

		/* We store the initial references for later reversing
		*/
		char *rowSequenceIterator = resultingSequenceForRow;
		char *columnSequenceIterator = resultingSequenceForColumn;

		while( backtrackMatrix[index] != STOP && (scoringOutcomeMatrix[index] > 0) )
		{
			/* size_t should be an integer type !
			*/

			size_t uxColumn = index % numberOfColumns;
			size_t uxRow = index / numberOfColumns;
			assert( uxColumn > 0 && uxRow > 0 );

			/* We check whether the two symbols match
			* WARNING: The scoring matrix is currently shifted into X as well Y direction by one position
			*/
			if ( puxRowSequenceRef[uxRow - 1] == puxColumnSequenceRef[uxColumn - 1] )
			{
				*(rowSequenceIterator++) = NucSeq::translateACGTCodeToCharacter(puxRowSequenceRef[uxRow - 1]);
				*(columnSequenceIterator++) = NucSeq::translateACGTCodeToCharacter(puxRowSequenceRef[uxRow - 1]);
				ruiNumberMatches++;
			} // if
			else
			{	
				switch ( backtrackMatrix[index] )
				{
				case LEFT_UP :
					ruiNumberMismatches++;
					*(rowSequenceIterator++) = 'x';
					break;
				case LEFT :
					ruiNumberDeletions++;
					*(rowSequenceIterator++) = '+';
					break;
				case UP :
					ruiNumberInsertions++;
					*(rowSequenceIterator++) = NucSeq::translateACGTCodeToCharacter(puxRowSequenceRef[uxRow - 1]);
				} // switch

				switch ( backtrackMatrix[index] )
				{
				case LEFT_UP :
					*(columnSequenceIterator++) = 'x';
					break;
				case LEFT :
					*(columnSequenceIterator++) = NucSeq::translateACGTCodeToCharacter(puxColumnSequenceRef[uxColumn - 1]);
					break;
				case UP :
					*(columnSequenceIterator++) = '+';
				} // switch
			} // else

			  /* Now we do the real backtracking.
			  */
			endIndex= index;
			switch ( backtrackMatrix[index] )
			{
			case LEFT_UP :
				index = ((uxRow - 1) * numberOfColumns) + (uxColumn - 1);
				break;
			case LEFT :
				index--;
				break;
			case UP :
				index = ((uxRow - 1) * numberOfColumns) + uxColumn;
			} // switch
		}
	} // method

	  /* Performs a backtrack within the scoring matrix. The outcome is here an STL vector. 
	  * The caller has to deliver an empty vector and is responsible for the memory allocation and deallocation.
	  */
	void backtrackFromIndex( T_size_t startIndex,
							 alignment_description<char> &alignmentOutcomeVector,
							 T_size_t &startPositionInColumn, 
							 T_size_t &endPositionInColumn,
							 T_size_t &startPositionInRow, 
							 T_size_t &endPositionInRow
	)
	{
		T_size_t currentIndex = startIndex;
		T_size_t previousIndex = startIndex;

		while( backtrackMatrix[currentIndex] != STOP && (scoringOutcomeMatrix[currentIndex] > 0) )
		{
			/* We calculate the current column and row on the foundation of the current index
			*/
			T_size_t uxColumn = currentIndex % numberOfColumns;
			T_size_t uxRow = currentIndex / numberOfColumns;
			assert( uxColumn > 0 && uxRow > 0 );

			/* We check whether the two symbols match
			* WARNING: The scoring matrix is shifted into X as well Y direction by one position
			*/
            switch ( backtrackMatrix[currentIndex] )
            {
            case LEFT_UP :
                if ( puxRowSequenceRef[uxRow - 1] == puxColumnSequenceRef[uxColumn - 1] )
                {
                    /* They match, so we insert a match into the vector. (push_back using &&)
                    */
                    alignmentOutcomeVector.push_back( alignment_description_element<char> ( EQUAL_PAIR, 
                                                    NucSeq::translateACGTCodeToCharacter( puxRowSequenceRef[uxRow - 1] ), 
                                                    NucSeq::translateACGTCodeToCharacter( puxColumnSequenceRef[uxColumn - 1] ) ) );
                } // if
                else
                {	
                    alignmentOutcomeVector.push_back( alignment_description_element<char> ( UNEQUAL_PAIR, 
                                                    NucSeq::translateACGTCodeToCharacter( puxRowSequenceRef[uxRow - 1] ), 
                                                    NucSeq::translateACGTCodeToCharacter( puxColumnSequenceRef[uxColumn - 1] ) ) );
                } // else
                break;
            case LEFT :
                alignmentOutcomeVector.push_back( alignment_description_element<char> ( INSERTION_AT_ROW_SIDE, 
                                                    '+', 
                                                    NucSeq::translateACGTCodeToCharacter( puxColumnSequenceRef[uxColumn - 1] ) ) );
                break;
            case UP :
                alignmentOutcomeVector.push_back( alignment_description_element<char> ( INSERTION_AT_COLUMN_SIDE, 
                                                    NucSeq::translateACGTCodeToCharacter( puxRowSequenceRef[uxRow - 1] ), 
                                                    '+' ) );
            case STOP :
                ;
            } // switch

			  /* Now we do the real backtracking.
			  */
			previousIndex = currentIndex;
			switch ( backtrackMatrix[currentIndex] )
			{
			case LEFT_UP :
				currentIndex = ((uxRow - 1) * numberOfColumns) + (uxColumn - 1);
				break;
			case LEFT :
				currentIndex--;
				break;
			case UP :
				currentIndex = ((uxRow - 1) * numberOfColumns) + uxColumn;
			case STOP :
				;
			}
		}

		/* Finally we inform about the begin and end of our sequences
		* Here we take 1 way, so we get values according to a counting starting with 0 instead of 1
		*/
        assert(previousIndex % numberOfColumns > 0);
		startPositionInColumn = (previousIndex % numberOfColumns) - 1;
		endPositionInColumn = (startIndex % numberOfColumns) - 1;

		startPositionInRow = (previousIndex / numberOfColumns) - 1;
		endPositionInRow = (startIndex / numberOfColumns) - 1;

		/* !!!WARNING!!! alignmentOutcomeVector was created in REVERSED order.
		* So, we reverse the alignmentOutcomeVector for getting the correct output
		*/
		std::reverse( alignmentOutcomeVector.begin(), alignmentOutcomeVector.end() );
	}

	void dumpAlignmentFromRowColumn( SCORE_TP score, T_size_t row, T_size_t column )
	{
		//// std::cout << "row is here " << row << " and column is here " << column << std::endl;
		T_size_t index = ((row + 1) * numberOfColumns) + (column + 1);
		assert( scoringOutcomeMatrix[index] == score );
		//// std::cout << "index is here " << index;
		dumpAlignmentFromIndex( index, scoringOutcomeMatrix[index] );
	}

	void dumpAlignmentFromIndex( T_size_t index, SCORE_TP score )
	{
		auto bufferSize = numberOfColumns + numberOfRows + 1;

		//// char *textBufferforRow = (char *)malloc( bufferSize * sizeof(char) );
		//// char *textBufferforColumn = (char *)malloc( bufferSize * sizeof(char) );

		std::vector<char> textBufferforRow( bufferSize );
		std::vector<char> textBufferforColumn( bufferSize );

		size_t uiNumberMatches = 0;
		size_t uiNumberMismatches = 0;
		size_t uiNumberInsertions = 0;
		size_t uiNumberDeletions = 0;

		//// std::cout << "index " << index << " has score " << score << " :: " << std::endl;
		backtrackFromIndexText( index,
								&textBufferforRow[0],
								&textBufferforColumn[0],
								uiNumberMatches,
								uiNumberMismatches,
								uiNumberInsertions,
								uiNumberDeletions);
		//// startColumn, endColumn, startRow, endRow);
		std::cout << "Query: " << (char*)(&textBufferforColumn[0]) << "|\n"
			<< "Refer: " << (char*)(&textBufferforRow[0]) << "|\n" 
			<< "Statistic: Matches " << uiNumberMatches << "  Mismatches " << uiNumberMismatches 
			<< "  Insertions " << uiNumberInsertions << "  Deletions " << uiNumberDeletions << std::endl;
	}

	/* Dumps all alignments for the current matrix
	*/
	void dumpAllAlignments() 
	{
		T_size_t matrixSize = numberOfRows * numberOfColumns;
		sw_alignment_result_t<SCORE_TP>* swAlignmentOutcomes = new sw_alignment_result_t<SCORE_TP>[matrixSize];	

		for(T_size_t uxIterator = 0; uxIterator < matrixSize; uxIterator++ )
		{
			swAlignmentOutcomes[uxIterator].score = scoringOutcomeMatrix[ uxIterator ];
			swAlignmentOutcomes[uxIterator].uxIndex = uxIterator;
		}

		/* We take the C like quick sort because it is faster than the stl implementation
		*/
		qsort( swAlignmentOutcomes, matrixSize, sizeof(sw_alignment_result_t<SCORE_TP>), compare_sw_alignment_result_t<SCORE_TP> );

		for(T_size_t uxIterator = 0; uxIterator < matrixSize; uxIterator++ )
		{
			if ( swAlignmentOutcomes[uxIterator].score > 0 )
			{
				dumpAlignmentFromIndex( swAlignmentOutcomes[uxIterator].uxIndex, swAlignmentOutcomes[uxIterator].score );
			}
		} // for

		delete[] swAlignmentOutcomes;
	}
};

/* queryOutcomeMatrix must have the size (numberOfRows + 1) x (numberOfColumns + 1)
* CONF_FILL_OUTCOME_MATRIX if you deliver true here, then the code for filling the outcome matrix is included
*/
template<bool CONF_FILL_OUTCOME_MATRIX, // Compile switch, that decides about filling outcome matrix
	class SCORE_TP				// the type that shall be used for the scoring
>
struct SW_align_type
{
	size_t numberOfColumns; const uint8_t *columnSequence;
	size_t numberOfRows; const uint8_t *rowSequence;

	/* number of symbols in your alphabet (don't forget +1 -> A,C,G,T => 5)
	*/
	const int alphabetSize; 

	/* Local sub-object that is initialized in the context of the default constructor.
	*/
	AlignmentOutcomeMatrix<SCORE_TP, size_t> alignmentOutcomeMatrix;

	/* Matrix that keeps the similarity weights.
	* In the moment it is initialized in the context of the constructor.
	*/
	SimilarityMatrix<SCORE_TP> similarityMatrix;

	SmithWatermanParamaterSet<SCORE_TP> &pSWparameterSetRef;

	/* The parameter
	* columnSequence : already translated sequence of numbers, representing the query sequence
	* rowSequence	  : already translated sequence of numbers, representing the database sequence
	* gapo			  : gap penalty for introducing a gap (together with gape)
	* gape			  : gap penalty for extending a gap
	*/
	SW_align_type( size_t numberOfColumns, 
				   const uint8_t *columnSequence, 
				   size_t numberOfRows, 
				   const uint8_t *rowSequence, 
				   SmithWatermanParamaterSet<SCORE_TP> &SWparameterSet
	) : 
		numberOfColumns(numberOfColumns),
		columnSequence(columnSequence),
		numberOfRows(numberOfRows),
		rowSequence(rowSequence),
		alphabetSize( SWparameterSet.uiAlphabetSize ),
		/* Initialize the scoring matrix
		*/
		alignmentOutcomeMatrix(numberOfColumns, columnSequence, numberOfRows, rowSequence ),
		similarityMatrix(SWparameterSet.iWeightMatch, SWparameterSet.iWeightMismatch, alphabetSize),
        pSWparameterSetRef( SWparameterSet )
	{ } // constructor

		/* The central alignment method
		*/
	SCORE_TP swAlign( 
#if (CONF_BAND_LIMITATION == 1)
		int w, 
#endif

#if (CONF_SET_INITIAL_SCORE_VALUE == 1)
		int initalScoreValue,
#endif
		size_t *pColumnIndexOfMaxScore, 
		size_t *pRowIndexOfMaxScore,
		size_t &indexOfMaxElementInScoringTable,
		std::vector<std::pair<size_t, size_t>> &rvMaxScorePositions // vector that receives the max-positions
	)
	{
		size_t uxIteratorRow; 
		size_t uxIteratorColumn, start, end;
		SCORE_TP maxScoreValue;
		long long rowIndexOfMaxScore, columnIndexOfMaxScore; 

		SCORE_TP gapoe = pSWparameterSetRef.iGapOpen + pSWparameterSetRef.iGapExtend;

		int scoreTableRowIterator;
		int queryProfileIterator;

#if 0
		std::cout << "One element of scoring values has " << sizeof(T_scoring) << " bytes" << std::endl;
#endif

#if (CONF_BAND_LIMITATION == 1)
		int sizeOfMatrix, upperLimitOfGapValue;
#endif

#if (CONF_SET_INITIAL_SCORE_VALUE == 1)
		if (initalScoreValue < 0)
			initalScoreValue = 0;
#endif

		/* allocate memory
		* 1. Query Profile (see ppt-file)
		* 2. h_and_e_Columns : Two one dimensional array for intermediate data storage (see ppt-file)
		*/
		SCORE_TP* queryProfile = (SCORE_TP *)malloc( numberOfColumns * alphabetSize * sizeof(SCORE_TP));
		generic_eh_type<SCORE_TP> *h_and_e_Columns = (generic_eh_type<SCORE_TP> *)calloc( numberOfColumns + 1, sizeof(generic_eh_type<SCORE_TP>) ); // memory for the score array

																																					/* Generate the query profile by using the scoring table
																																					*/
		queryProfileIterator = 0;
		for (scoreTableRowIterator = 0; scoreTableRowIterator < alphabetSize; ++scoreTableRowIterator) {
			SCORE_TP *referenceToStartOfRow = &similarityMatrix.pSimilarityMatrixRef[scoreTableRowIterator * alphabetSize];

			for (uxIteratorColumn = 0; uxIteratorColumn < numberOfColumns; ++uxIteratorColumn) {
				queryProfile[queryProfileIterator] = referenceToStartOfRow[ columnSequence[uxIteratorColumn] ];
				queryProfileIterator++;
			} // for
		} // for

		  /* Initialize the h-elements of the scoreArray.
		  * The e-elements are set 0 by the initialization
		  */
#if (CONF_SET_INITIAL_SCORE_VALUE == 1)
		h_and_e_Columns[0].h = initalScoreValue; 
		h_and_e_Columns[1].h = initalScoreValue > gapoe ? initalScoreValue - gapoe 
			: 0;
#else
		h_and_e_Columns[0].h = 0;
		h_and_e_Columns[1].h = 0;
#endif

		for( uxIteratorColumn = 2; 

			 uxIteratorColumn <= numberOfColumns
			 && h_and_e_Columns[uxIteratorColumn-1].h > pSWparameterSetRef.iGapExtend;

			 ++uxIteratorColumn)
		{
			h_and_e_Columns[uxIteratorColumn].h = h_and_e_Columns[uxIteratorColumn-1].h - pSWparameterSetRef.iGapExtend;
		} // for


#if (CONF_BAND_LIMITATION == 1)
		  /* adjust $w if it is too large
		  * gape = 1; gapo = 3; gapoe = 4;
		  */
		sizeOfMatrix = alphabetSize * alphabetSize;
		maxScoreValue = 0;
		for (uxIteratorRow = 0;  uxIteratorRow < sizeOfMatrix; ++uxIteratorRow) // get the max score
			maxScoreValue = maxScoreValue > similarityMatrix[uxIteratorRow] ? maxScoreValue 
			: similarityMatrix[uxIteratorRow];

		upperLimitOfGapValue = (int)((double)(numberOfColumns * maxScoreValue - gapo) / gape + 1.);

		upperLimitOfGapValue = upperLimitOfGapValue > 1 ? upperLimitOfGapValue 
			: 1;

		w = w < upperLimitOfGapValue ? w 
			: upperLimitOfGapValue;
#endif

		/* computation of H matrix
		*/
#if (CONF_SET_INITIAL_SCORE_VALUE == 1)
		maxScoreValue = initalScoreValue;
#else
		maxScoreValue = 0;
#endif

		/* rowIndexOfMaxScore and columnIndexOfMaxScore were originally initialized to -1 !
		* they should be of type int instead of size_t
		*/
		rowIndexOfMaxScore = 0; 
		columnIndexOfMaxScore = 0;
		start = 0;
		end = numberOfColumns;

		for (uxIteratorRow = 0; uxIteratorRow < numberOfRows; ++uxIteratorRow) {
			SCORE_TP f = 0;
			SCORE_TP hInNextRound;
			long long indexOfMaxScoreInCurrentRow = -1;
			SCORE_TP maxScoreInCurrentRow = 0;

			/* We select the row in our query profile according to the row-sequence symbol that belongs to the current row
			*/
			SCORE_TP *referenceToQueryProfileRow = &queryProfile[ rowSequence[uxIteratorRow] * numberOfColumns ];

			/* compute the first column
			*/
#if (CONF_SET_INITIAL_SCORE_VALUE == 1)
			hInNextRound = initalScoreValue - (gapo + gape * (uxIteratorRow + 1));
			if (hInNextRound < 0) 
				hInNextRound = 0;
#else
			hInNextRound = 0;
#endif

#if (CONF_BAND_LIMITATION == 1)
			/* apply the band and the constraint (if provided)
			* w seem to a 'band' that limits the computational space to +- w elements
			*/
			if (start < uxIteratorRow - w) 
				start = uxIteratorRow - w;

			if (end > uxIteratorRow + w + 1) 
				end = uxIteratorRow + w + 1;

			if (end > numberOfColumns) 
				end = numberOfColumns;
#endif

			/* The outcome Matrix has size (numberOfRows + 1) x (numberOfColumns + 1)
			*/
			auto uxMatrixStartShift = CONF_FILL_OUTCOME_MATRIX ? ( ( (uxIteratorRow + 1) * (numberOfColumns + 1) ) + (start + 1) )
															   : 0;

			SCORE_TP *queryOutcomeMatrixIterator = CONF_FILL_OUTCOME_MATRIX ? alignmentOutcomeMatrix.scoringOutcomeMatrix + uxMatrixStartShift
																			: NULL;

			sw_direction_t *backtrackingMatrixIterator = CONF_FILL_OUTCOME_MATRIX ? alignmentOutcomeMatrix.backtrackMatrix + uxMatrixStartShift
																				  : NULL;

			for ( uxIteratorColumn = start; uxIteratorColumn < end; ++uxIteratorColumn ) {
				/* At the beginning of the loop: scoreArray[j] = { H(i - 1, j - 1), E(i, j) }, 
				*								 f = F(i, j) 
				*                               hInNextRound = new value for H(i, j - 1)
				*/
				generic_eh_type<SCORE_TP> *refIntoScoreArrayForColumn = &h_and_e_Columns[uxIteratorColumn];
				SCORE_TP h = refIntoScoreArrayForColumn->h; // get H(i-1,j-1)
				SCORE_TP e = refIntoScoreArrayForColumn->e; // get E(i,j)

															/* 1. add match or mismatch distance
															* 2. h = H(i,j)= max{H(i-1, j-1) + Score(i, j), E(i, j), F(i, j)}
															* Here we can imply the entry in the H matrix as well as the backtracking direction !
															*/
				sw_direction_t eDirection;
				if ( CONF_FILL_OUTCOME_MATRIX )
				{
					eDirection = LEFT_UP;
				}

				h += referenceToQueryProfileRow[uxIteratorColumn]; // H(i-1, j-1) + Score(i, j)

                /*
                 * NOTE: the >= is very important for the backtracking
                 */
				if (e >= h)
				{
					h = e;
					if ( CONF_FILL_OUTCOME_MATRIX )
					{
						eDirection = UP;
					}
				}

                /*
                 * NOTE: the >= is very important for the backtracking
                 */
				if (f >= h)
				{
					h = f;
					if ( CONF_FILL_OUTCOME_MATRIX )
					{
						eDirection = LEFT;
					}
				}

				/* We store h in the Matrix for later analysis
				*/
				if ( CONF_FILL_OUTCOME_MATRIX )
				{
					*(queryOutcomeMatrixIterator++) = h;
					*(backtrackingMatrixIterator++) = eDirection;
				}

				refIntoScoreArrayForColumn->h = hInNextRound; // set H(i,j-1) for the next row
				hInNextRound = h; // save H(i,j) to hInNextRound for storage in the next iteration.

				indexOfMaxScoreInCurrentRow = maxScoreInCurrentRow > h ? indexOfMaxScoreInCurrentRow 
					: uxIteratorColumn;

				maxScoreInCurrentRow = maxScoreInCurrentRow > h ? maxScoreInCurrentRow 
					: h;   // m is stored at eh[mj+1]

						   /* Compute new e value and store this value for the next round
						   * E(i + 1,j) = max{H(i, j) - gapo, E(i,j)} - gape
						   */
				h -= gapoe;
				h = h > 0 ? h 
					: 0;
				e -= pSWparameterSetRef.iGapExtend;
				e = e > h ? e 
					: h;   
				refIntoScoreArrayForColumn->e = e; // save E(i + 1, j) for the next row

												   /* Compute new f value and store this for the next iteration
												   * F(i,j+1) = max{H(i, j) - gapo, F(i,j)} - gape
												   */
				f -= pSWparameterSetRef.iGapExtend;
				f = f > h ? f 
					: h;  
			} // inner for

			  /* Save the final hInNextRound, because there are no more iterations.
			  * Seems to be necessary for begin and end management
			  */
#if (CONF_FOLLOW_MAX_PATH_OPTIMIZATION == 1)
			h_and_e_Columns[end].h = hInNextRound; 
			h_and_e_Columns[end].e = 0;
#endif
#if 0
			if (maxScoreInCurrentRow == 0) 
			{	/* There is no reason to continue, because all scores were zero
				*/
				break;
			}
#endif
			/* Logging of the maximum scores */
			if( maxScoreInCurrentRow >= maxScoreValue )
			{
				if( maxScoreInCurrentRow > maxScoreValue )
				{	// fresh overall maximum detected 
					rvMaxScorePositions.clear();
					maxScoreValue = maxScoreInCurrentRow;
                    columnIndexOfMaxScore = indexOfMaxScoreInCurrentRow;//added by markus
                    rowIndexOfMaxScore = uxIteratorRow;//added by markus
				} // if
				  /* Log pairs (row, column (max-pos within row) )*/
				rvMaxScorePositions.push_back( std::pair<size_t, size_t>( uxIteratorRow, indexOfMaxScoreInCurrentRow) );
			} // if

#if (CONF_FOLLOW_MAX_PATH_OPTIMIZATION == 1)
			  /* We count down from the maximum position until we find score 0 or we reach the first column
			  */
			for (uxIteratorColumn = indexOfMaxScoreInCurrentRow; uxIteratorColumn >= start && h_and_e_Columns[uxIteratorColumn].h; --uxIteratorColumn)
				;
			start = uxIteratorColumn + 1;

			/* We do something similar for the end
			*/
			for (uxIteratorColumn = indexOfMaxScoreInCurrentRow + 2; uxIteratorColumn <= end && h_and_e_Columns[uxIteratorColumn].h; ++uxIteratorColumn)
				;
			end = uxIteratorColumn;
#endif

			//beg = 0; end = qlen; // uncomment this line for debugging
		} // outer for

		  /* Free all allocated memory
		  */
		free(h_and_e_Columns); 
		free(queryProfile);

		/* We return the position (row and column index) of the maximum score value.
		* The values are relative to the matrix. (+1 applied)
		*/
		if ( pColumnIndexOfMaxScore )
		{
			*pColumnIndexOfMaxScore = columnIndexOfMaxScore + 1;
		} // if
		if ( pRowIndexOfMaxScore )
		{
			*pRowIndexOfMaxScore = rowIndexOfMaxScore + 1;
		} // if

        
        DEBUG(
            if( pColumnIndexOfMaxScore && pRowIndexOfMaxScore )
            {
                std::cout << "Max score is: " << maxScoreValue << std::endl;
            } // if
        )

		  /* The index of the maximum within the table.
		  */
		indexOfMaxElementInScoringTable = (rowIndexOfMaxScore + 1) * alignmentOutcomeMatrix.numberOfColumns + (columnIndexOfMaxScore + 1);
		return maxScoreValue;
	}
}; // struct