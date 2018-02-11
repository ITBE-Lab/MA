// SwithWaterman.cpp : Defines the entry point for the console application.
//

#include <string>
#include <iostream>
#include <algorithm>
#include "container/nucSeq.h" // sequence slices
using namespace libMA;
#include "module/sw_sse_avx.h"
#include "module/smith_waterman.h"
#include "module/SW_sequential.h"


extern int iGap;
extern int iExtend;
extern int iMatch;
extern int iMissMatch;

ContainerVector SMW::getInputType() const
{
    return ContainerVector{
            //the query sequence
            std::shared_ptr<Container>(new NucSeq()),
            //the ref
            std::shared_ptr<Container>(new NucSeq())
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
int16_t alignSW_SIMD( NucSeq &rQuerySequence, // query sequence
                      NucSeq &rReferenceSequence, // reference sequence
                      SmithWatermanParamaterSet<int16_t> &rSWparameterSet, // Smith Waterman alignment parameter
                      std::vector<size_t> &rvMaxScorePositions // vector will recieve positions, where we have a max score
                    )
{
    //reversing query and reference in order to obtain start instead of end (markus)
    rQuerySequence.vReverse();
    rReferenceSequence.vReverse();

    int16_t iMaxScore = 0;

    try
    {    // Create aligner object on the foundation of the query
        SW_SIMD_Aligner<int16_t, size_t> xSIMDAligner( rQuerySequence, rSWparameterSet );
        
        /* Do the actual alignment using a reference.
         * iMaxScore is the maximum score computed in the context of the alignment.
         * rvMaxScorePositions is a vector that contains all positions having maxscore.
         * Position counting starts with zero. (First row has index zero)
         */
        iMaxScore = xSIMDAligner.align( rReferenceSequence, rvMaxScorePositions );
    } // try
    catch (std::exception &e)
    {    /* Aligner can throw expcetions of type AlignerException
         */
        std::cout << "Catched exception: " << e.what();
        exit(0);
    } // catch

    return iMaxScore;
} // function

extern int iGap;
extern int iExtend;
extern int iMatch;
extern int iMissMatch;

std::shared_ptr<Container> SMW::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<NucSeq> pQuerySeq = 
        std::static_pointer_cast<NucSeq>((*vpInput)[0]);
    std::shared_ptr<NucSeq> pReference = 
        std::static_pointer_cast<NucSeq>((*vpInput)[1]);
    //std::shared_ptr<NucSeq> pReference = pRefPack->vColletionAsNucSeq();

    if(pQuerySeq->length() == 0 || pReference->length() == 0)
        return std::shared_ptr<ContainerVector>(new ContainerVector());

    DEBUG(
        std::cout << pQuerySeq->length() << "x" << pReference->length() << std::endl;
    )//debug

    // 1. Prepare the SW parameter set ...
    SmithWatermanParamaterSet<int16_t> xSWparameterSet
    (   iMatch, // score for match (must be positive)
        -iMissMatch, // score for mismatch (must be negative)
        iGap, // penalty for gap open ("gap" means insertion or deletion)
        iExtend, // penalty for gap extension
        pQuerySeq->uxAlphabetSize() // alphabet size of input
    );

    if(bBacktrack)
    {
        /* Create aligner object */
        SW_align_type <
            true, // create data for backtracking
            int16_t> // forward type for scoring
            xAligner( pQuerySeq->length(), // length query
                    pQuerySeq->pGetSequenceRef(), // address first query symbol
                    pReference->length(), // length reference
                    pReference->pGetSequenceRef(), // address first reference symbol
                    xSWparameterSet ); // SW-parameter

        /* Prepare variables that will receive aligner output */
        size_t uiColumnIndexOfMaxScore; // Receives column of max-score
        size_t uiRowIndexOfMaxScore; // Receives row of max-score
        size_t indexOfMaxElementInScoringTable; // linear index for max-position in matrix
        std::vector<std::pair<size_t, size_t>> vMaxScorePositions; // receives the maximum positions as tuples (row, col)

        /* Do the actual alignment... */
        xAligner.swAlign(
            &uiColumnIndexOfMaxScore,
            &uiRowIndexOfMaxScore,
            indexOfMaxElementInScoringTable,
            vMaxScorePositions
        );
        // collect the results ...
        auto pvRet = std::shared_ptr<ContainerVector>(
            new ContainerVector(std::shared_ptr<Alignment>(new Alignment()))
        );

        alignment_description<char> vTemp;
        size_t startPositionQuery = 0;
        size_t endPositionQuery = 0;
        size_t startPositionRef = 0;
        size_t endPositionRef = 0;
        xAligner.alignmentOutcomeMatrix.backtrackFromIndex(
            indexOfMaxElementInScoringTable,
            vTemp,
            startPositionQuery,
            endPositionQuery,
            startPositionRef,
            endPositionRef
        );
        std::shared_ptr<Alignment> pAlignment = std::shared_ptr<Alignment>(
            new Alignment(startPositionRef, endPositionRef+1, startPositionQuery, endPositionQuery+1));
        for(auto pos : vTemp)
            switch (pos.eElementKind){
                case EQUAL_PAIR:
                    pAlignment->append(MatchType::match);
                    break;
                case UNEQUAL_PAIR:
                    pAlignment->append(MatchType::missmatch);
                    break;
                case INSERTION_AT_ROW_SIDE://ref -> deletion
                    pAlignment->append(MatchType::insertion);
                    break;
                case INSERTION_AT_COLUMN_SIDE://query -> insertion
                    pAlignment->append(MatchType::deletion);
                    break;
                default:
                    std::cout << "unkown symbol in alignment" << std::endl;
            }//switch
        //for
        pvRet->push_back(pAlignment);
        return pvRet;
    }//if
    else
    {
        // 2. Do the alignment ...
        std::vector<size_t> vMaxScorePositions;
        int16_t imaxScore = alignSW_SIMD( *pQuerySeq, // query sequence
                                        *pReference, // reference sequence
                                        xSWparameterSet, // Smith Waterman alignment parameter
                                        vMaxScorePositions // vector will recieve positions, where we have a max score
                                        );

        // 3. collect the results ...
        auto pvRet = std::make_shared<ContainerVector>();

        for(nucSeqIndex revStart : vMaxScorePositions)
        {
            //undo the reversion by substracting from absolute length
            std::shared_ptr<Alignment> pAlignment = std::shared_ptr<Alignment>(new Alignment(
                pReference->length() - revStart));
            pAlignment->iScore = imaxScore;
            pvRet->push_back(pAlignment);
        }//for

        return pvRet;
    }//else
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
    (   iMatch, // score for match (must be positive)
        -iMissMatch, // score for mismatch (must be negative)
        iGap, // penalty for gap open ("gap" means insertion or deletion)
        iExtend, // penalty for gap extension
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
    {    // Sequence at max score position must be equal to query
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
            boost::python::bases<Module>,
            std::shared_ptr<SMW>
            >("SMW",boost::python::init<bool>())
        ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<SMW>,
        std::shared_ptr<Module> 
    >();
}