/** 
 * @file needlemanWunsch.cpp
 * @author Markus Schmidt
 */
#include "module/needlemanWunsch.h"
#include <bitset>
#include "ksw/ksw2.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>


using namespace libMA;

using namespace libMA::defaults;
extern int libMA::defaults::iGap;
extern int libMA::defaults::iExtend;
extern int libMA::defaults::iGap2;
extern int libMA::defaults::iExtend2;
extern int libMA::defaults::iMatch;
extern int libMA::defaults::iMissMatch;
extern nucSeqIndex libMA::defaults::uiMaxGapArea;
extern nucSeqIndex libMA::defaults::uiPadding;

DEBUG(
    bool bAnalyzeHeuristics = false;
)// DEBUG

//match missmatch matrix for ksw
int8_t mat[25];

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = 0;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = 0;
}// function

void ksw_simplified(
        int qlen, const uint8_t *query,
        int tlen, const uint8_t *target,
        int8_t q, int8_t e, int8_t q2, int8_t e2, int& w,
        ksw_extz_t *ez
    )
{
    int minAddBandwidth = 10; // must be >= 0 otherwise ksw will not align till the end
    if( std::abs(tlen - qlen) + minAddBandwidth > w)
        w = std::abs(tlen - qlen) + minAddBandwidth;
    ksw_extd2_sse(nullptr, qlen, query, tlen, target, 5, mat, q, e, q2, e2, w, -1, -1, 0, ez);
}// function


//small wrapper that takes care of deallocation
class Wrapper_ksw_extz_t
{
public:
    ksw_extz_t* ez;

    Wrapper_ksw_extz_t()
    {
        ez = new ksw_extz_t {};// {} forces zero initialization

    }//default constructor

    ~Wrapper_ksw_extz_t()
    {
        delete ez->cigar;
        delete ez;
    }//default constructor
};//class

// banded global NW
void ksw(
        std::shared_ptr<NucSeq> pQuery,
        std::shared_ptr<NucSeq> pRef,
        nucSeqIndex fromQuery,
        nucSeqIndex toQuery,
        nucSeqIndex fromRef,
        nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment
    )
{
    //sanity checks
    if(toRef <= fromRef)
        if(toQuery <= fromQuery)
            return;
    if(toQuery <= fromQuery)
    {
        pAlignment->append(MatchType::deletion, toRef-fromRef);
        return;
    }//if
    if(toRef <= fromRef)
    {
        pAlignment->append(MatchType::insertion, toQuery-fromQuery);
        return;
    }//if

    Wrapper_ksw_extz_t ez;

    int uiBandwidth = 100;

    assert(toQuery < pQuery->length());
    assert(toRef < pRef->length());
    ksw_simplified(
        toQuery - fromQuery, 
        pQuery->pGetSequenceRef() + fromQuery,
        toRef - fromRef, 
        pRef->pGetSequenceRef() + fromRef, 
        iGap,
        iExtend,
        iGap2,
        iExtend2,
        uiBandwidth,
        ez.ez // return value
    );

    uint32_t qPos = fromQuery;
    uint32_t rPos = fromRef;
    for (int i = 0; i < ez.ez->n_cigar; ++i)
    {
        uint32_t uiSymbol = ez.ez->cigar[i]&0xf;
        uint32_t uiAmount = ez.ez->cigar[i]>>4;
        switch (uiSymbol)
        {
            case 0:
                for(uint32_t uiPos = 0; uiPos < uiAmount; uiPos++)
                {
                    if( (*pQuery)[uiPos + qPos] == (*pRef)[uiPos + rPos] )
                        pAlignment->append(MatchType::match);
                    else
                        pAlignment->append(MatchType::missmatch);
                }// for
                qPos+=uiAmount;
                rPos+=uiAmount;
                break;
            case 1:
                pAlignment->append(MatchType::insertion, uiAmount);
                qPos+=uiAmount;
                break;
            case 2:
                pAlignment->append(MatchType::deletion, uiAmount);
                rPos+=uiAmount;
                break;
            default:
                std::cerr << "obtained wierd symbol from ksw: " << uiSymbol << std::endl;
                assert(false);
                break;
        }//switch
    }//for
#if DEBUG_LEVEL >= 1
        const char vMIDN[] = {'M', 'I', 'D', 'N'};
        if(qPos != (uint32_t)toQuery && rPos != (uint32_t)toRef)
        {
            std::cerr << "ksw did neither extend till end of query nor ref " 
                << (int)toQuery - qPos << "; " << (int)toRef -  rPos
                << " bandwidth: " << uiBandwidth << std::endl;
            std::cout << pQuery->fromTo(fromQuery, toQuery) << std::endl;
            std::cout << pRef->fromTo(fromRef, toRef) << std::endl;
            std::cout << "CIGAR:";
            for (int i = 0; i < ez.ez->n_cigar; ++i)
            {
                uint32_t uiSymb = ez.ez->cigar[i]&0xf;
                uint32_t uiAmount = ez.ez->cigar[i]>>4;
                std::cout << uiAmount <<  vMIDN[uiSymb] << " ";
            }//for
            std::cout << std::endl;
            assert(false);
        }// if
#endif
    //ensure we do not append something negative.. (underflow)
    assert(toQuery >= qPos);
    assert(toRef >= rPos);
    //Append the remaining insertion or deletion. 
    //Remaining lengths of zero are taken care of in append
    pAlignment->append(MatchType::deletion, toQuery - qPos);
    pAlignment->append(MatchType::insertion, toRef - rPos);
}//function


/**
 * @brief the NW dynamic programming algorithm
 * @details
 * This is a naive very slow version,
 * but it computes the correct score and it is capable of doing semi-global alignments.
 *
 * if bNoGapAtBeginning || bNoGapAtEnd :
 *      returns the gap at the beginning or end
 * otherwise : returns 0 
 */
void NeedlemanWunsch::naiveNeedlemanWunsch(
        std::shared_ptr<NucSeq> pQuery, 
        std::shared_ptr<NucSeq> pRef,
        nucSeqIndex fromQuery,
        nucSeqIndex toQuery,
        nucSeqIndex fromRef,
        nucSeqIndex toRef,
        bool bNoGapAtBeginning,
        bool bNoGapAtEnd,
        std::shared_ptr<Alignment> pAlignment
    )
{
    /*
    * break conditions for actually empty areas
    */
    DEBUG(
        if(toQuery > pQuery->length())
        {
            std::cerr << toQuery << " " << pQuery->length() << std::endl;
            assert(false);
        }//if
        if(toRef > pRef->length())
        {
            std::cerr << toRef << " " << pRef->length() << std::endl;
            assert(false);
        }//if
    )// DEBUG
    if(toRef <= fromRef)
        if(toQuery <= fromQuery)
            return;
    DEBUG_2(
        std::cout << toQuery-fromQuery << std::endl;
        for(nucSeqIndex i = fromQuery; i < toQuery; i++)
            std::cout << pQuery->charAt(i);
        std::cout << std::endl;
        std::cout << toRef-fromRef << std::endl;
        for(nucSeqIndex i = fromRef; i < toRef; i++)
            std::cout << pRef->charAt(i);
        std::cout << std::endl;
    )//DEBUG
    if(toQuery <= fromQuery)
    {
        pAlignment->append(MatchType::deletion, toRef-fromRef);
        return;
    }//if
    if(toRef <= fromRef)
    {
        pAlignment->append(MatchType::insertion, toQuery-fromQuery);
        return;
    }//if

    /*
    * beginning of the actual NW
    */
#if ALLOCATE_ONCE == 0
    //allocate memory for the naive approach
    std::vector<std::vector<std::vector<int>>> s(
        3,
        std::vector<std::vector<int>>(
            toQuery-fromQuery+1,
            std::vector<int>(toRef-fromRef+1)
        )
    );
    std::vector<std::vector<std::vector<char>>> dir(
        3,
        std::vector<std::vector<char>>(
            toQuery-fromQuery+1,
            std::vector<char>(toRef-fromRef+1)
        )
    );
#endif

    /*
    * initialization:
    *      this part sets the scores for the last row and column (reverse order)
    * 
    * Note:
    *      if we do not want a gap at the end since the alignment ends there we need to 
    *      set the initial values along the reference to 0.
    *      we do not want a complete global alignment,
    *      merely a global alignment with respect to the query
    */
    //                  BINARY       DECIMAL
    #define DIA         /*000001*/   1
    #define INS         /*000010*/   2
    #define DEL         /*000100*/   4
    #define DIR_0_NEXT  /*001000*/   8
    #define DIR_1_NEXT  /*010000*/   16
    #define DIR_2_NEXT  /*100000*/   32

    // used to prevent the DP to make extensions from positions where this is set...
    // -iGap*10 should be more than enough
    #define LOWER -iGap*1000

    s[0][0][0] = 0;
    s[1][0][0] = LOWER;
    s[2][0][0] = LOWER;
    dir[0][0][0] = 0;// this position will throw an error is the backtracker tries to use it
    dir[1][0][0] = 0;// this position will throw an error is the backtracker tries to use it
    dir[2][0][0] = 0;// this position will throw an error is the backtracker tries to use it

    s[0][1][0] = LOWER;
    dir[0][1][0] = 0;// this position will throw an error is the backtracker tries to use it
    s[1][1][0] = - (iGap + iExtend);
    dir[1][1][0] = INS | DIR_0_NEXT;
    s[2][1][0] = LOWER;
    dir[2][1][0] = 0;// this position will throw an error is the backtracker tries to use it
    
    s[0][0][1] = LOWER;
    dir[0][0][1] = 0;// this position will throw an error is the backtracker tries to use it
    s[1][0][1] = LOWER;
    dir[1][0][1] = 0;// this position will throw an error is the backtracker tries to use it
    if(bNoGapAtEnd)//see note above
        s[2][0][1] = 0;
    else
        s[2][0][1] = - (iGap + iExtend);
    dir[2][0][1] = DEL | DIR_0_NEXT;
    for(unsigned int x=0; x < 3; x++)
    {
        for(nucSeqIndex uiI = 2; uiI < toQuery-fromQuery+1; uiI++)
        {
            if(bNoGapAtEnd)//see note above
                s[x][0][uiI] = 0;
            else
                s[x][uiI][0] = s[x][uiI - 1][0] - iExtend;
            dir[1][uiI][0] = INS | DIR_1_NEXT;
        }//for
        for(nucSeqIndex uiI = 2; uiI < toRef-fromRef+1; uiI++)
        {
            s[x][0][uiI] = s[x][0][uiI - 1] - iExtend;
            dir[2][0][uiI] = DEL | DIR_2_NEXT;
        }//for
    }//for
    /*
    * dynamic programming loop
    * Note:
    *      we iterate in the reverse order on reference and query
    *      so that the backtracking can be done in forward order
    *      this saves us the work to reverse the result
    *
    * This works as follows:
    *      for each cell compute the scores if resuling from an insertion deletion match/missmatch
    *      in this order. Store the score from the insertion and overwrite the score with the del
    *      match of missmatch score if any of them is higher. Also keep track of which direction
    *      we came from in the dir matrix.
    */
    int a, b;
    char c;
    for(nucSeqIndex uiI = 1; uiI < (toQuery-fromQuery)+1; uiI++)
    {
        for(nucSeqIndex uiJ = 1; uiJ < (toRef-fromRef)+1; uiJ++)
        {
            //match / missmatch
            a = s[0][uiI - 1][uiJ - 1];
            c = DIA | DIR_0_NEXT;
            b = s[1][uiI - 1][uiJ - 1];
            if(b > a)
            {
                a = b;
                c = DIA | DIR_1_NEXT;
            }//if
            b = s[2][uiI - 1][uiJ - 1];
            if(b > a)
            {
                a = b;
                c = DIA | DIR_2_NEXT;
            }//if
            if( (*pQuery)[toQuery - uiI] == (*pRef)[toRef - uiJ] )
                a += iMatch;
            else
                a -= iMissMatch;
            dir[0][uiI][uiJ] = c;
            s[0][uiI][uiJ] = a;

            //insertion
            a = s[1][uiI - 1][uiJ] - iExtend;
            c = INS | DIR_1_NEXT;
            b = s[0][uiI - 1][uiJ] - (iGap + iExtend);
            if(b >= a)
            {
                a = b;
                c = INS | DIR_0_NEXT;
            }//if
            b = s[2][uiI - 1][uiJ] - (iGap + iExtend);
            if(b > a)
            {
                a = b;
                c = INS | DIR_2_NEXT;
            }//if
            dir[1][uiI][uiJ] = c;
            s[1][uiI][uiJ] = a;

            //deletion
            a = s[2][uiI][uiJ - 1] - iExtend;
            c = DEL | DIR_2_NEXT;
            b = s[0][uiI][uiJ - 1] - (iGap + iExtend);
            if(b >= a)
            {
                a = b;
                c = DEL | DIR_0_NEXT;
            }//if
            b = s[1][uiI][uiJ - 1] - (iGap + iExtend);
            if(b > a)
            {
                a = b;
                c = DEL | DIR_1_NEXT;
            }//if
            dir[2][uiI][uiJ] = c;
            s[2][uiI][uiJ] = a;
        }//for
    }//for

    DEBUG_3(
        /*
        * sanity prints
        */
        for(nucSeqIndex uiI = 0; uiI < toRef-fromRef+1; uiI++)
        {
            if(uiI == 0)
                std::cout << " \t \t";
            else
                std::cout << pRef->charAt(toRef - uiI) << "\t";
        }//for
        std::cout << std::endl;
        for(nucSeqIndex uiI = 0; uiI < toQuery-fromQuery+1; uiI++)
        {
            if(uiI == 0)
                std::cout << " \t";
            else
                std::cout << pQuery->charAt(toQuery - uiI) << "\t";
            for(nucSeqIndex uiJ = 0; uiJ < toRef-fromRef+1; uiJ++)
                std::cout << s[2][uiI][uiJ] << "\t";
            std::cout << std::endl;
        }//for
    )//DEBUG

    /*
    * backtracking
    */
    nucSeqIndex iX = toQuery-fromQuery;
    nucSeqIndex iY = toRef-fromRef;
    
    char cLastDir = DIR_0_NEXT;
    /*
    * if there is no gap cost for the beginning 
    * we should start backtracking where the score is maximal
    * along the reference
    * also: in this case the last direction must be a match
    */
    if(bNoGapAtBeginning)
    {
        for(nucSeqIndex uiJ = 1; uiJ < toRef-fromRef; uiJ++)
            if(s[0][iX][uiJ] > s[0][iX][iY])
                iY = uiJ;
        DEBUG_2(
            std::cout << (toRef-fromRef) - iY << "D";
        )//DEBUG
        pAlignment->shiftOnRef((toRef-fromRef) - iY);
    }//if
    else // in this case the first direction might be an insertion or deletion
    {
        int a = s[0][iX][iY];
        int b = s[1][iX][iY];
        if(b > a)
        {
            a = b;
            cLastDir = DIR_1_NEXT;
        }//if
        b = s[2][iX][iY];
        if(b > a)
            cLastDir = DIR_2_NEXT;
    }// else
    while(iX > 0 || iY > 0)
    {
        //load the direction value from the correct matrix
        if(cLastDir & DIR_0_NEXT)
            cLastDir = dir[0][iX][iY];
        else if(cLastDir & DIR_1_NEXT)
            cLastDir = dir[1][iX][iY];
        else if(cLastDir & DIR_2_NEXT)
            cLastDir = dir[2][iX][iY];
        else
            std::cerr << "WARNING: no next pointer set in dynamic programming" << std::endl;
        //do the backtracking
        if(cLastDir & DIA)
        {
            if( (*pQuery)[toQuery - iX] == (*pRef)[toRef - iY] )
            {
                pAlignment->append(MatchType::match);
                DEBUG_2(
                    std::cout << "M";
                )//DEBUG
            }//if
            else
            {
                pAlignment->append(MatchType::missmatch);
                DEBUG_2(
                    std::cout << "W";
                )//DEBUG
            }//else
            iX--;
            iY--;
        }//if
        else if(cLastDir & INS)
        {
            pAlignment->append(MatchType::insertion);
            iX--;
            DEBUG_2(
                std::cout << "I";
            )//DEBUG
        }//if
        else if(cLastDir & DEL)
        {
            pAlignment->append(MatchType::deletion);
            iY--;
            DEBUG_2(
                std::cout << "D";
            )//DEBUG        
        }//if
        else{
            std::cerr << "WARNING: no direction set in dynamic programming" << std::endl;
        }//else
        /*
        * if there is no gap cost for the end 
        * we should stop backtracking once we reached the end of the query
        */
        if(bNoGapAtEnd && iX <= 0)
            return;
    }//while

    DEBUG_2(
        std::cout << std::endl;
    )//DEBUG

    //print the entire matrix if necessary
    DEBUG_3(
        std::cout << "\t";
        for(auto i = toRef; i > fromRef; i--)
            std::cout << "\t" << NucSeq::translateACGTCodeToCharacter((*pRef)[i - 1]);
        for(auto j = fromQuery; j <= toQuery; j++)
        {
            std::cout << "\n";
            if(j > fromQuery)
                std::cout << NucSeq::translateACGTCodeToCharacter((*pQuery)[toQuery - j]);
            for(auto i = fromRef; i <= toRef; i++)
                std::cout
                    << "\t"
                    << s[0][j - fromQuery][i - fromRef]
                    << ","
                    << s[1][j - fromQuery][i - fromRef]
                    << ","
                    << s[2][j - fromQuery][i - fromRef]
                    << " ("
                    << std::bitset<6>(dir[0][j - fromQuery][i - fromRef])
                    << ","
                    << std::bitset<6>(dir[1][j - fromQuery][i - fromRef])
                    << ","
                    << std::bitset<6>(dir[2][j - fromQuery][i - fromRef])
                    << ")"
                    ;
        }//for
        std::cout << std::endl;
    )//DEBUG
    return;
}//function

class MyPrinterMemory
{
public:
    std::shared_ptr<NucSeq> pQuery;
    std::shared_ptr<NucSeq> pRef;
    std::shared_ptr<Alignment> pAlignment;
    nucSeqIndex uiQPos, uiRPos;
    DEBUG(
        nucSeqIndex uiQueryExtensionSize = 0;
        nucSeqIndex uiRefExtensionSize = 0;
    )

    MyPrinterMemory(
            std::shared_ptr<NucSeq> pQuery,
            std::shared_ptr<NucSeq> pRef,
            std::shared_ptr<Alignment> pAlignment,
            nucSeqIndex uiQueryStart,
            nucSeqIndex uiRefStart
        )
            :
        pQuery(pQuery),
        pRef(pRef),
        pAlignment(pAlignment),
        uiQPos(uiQueryStart),
        uiRPos(uiRefStart)
    {}//constructor
};// class


char vTrans[9] = {'?', 'A', 'C', '?', 'G', '?', '?', '?', 'T'};

int printer(void* pVoid, uint64_t uiLen, char c)
{
    MyPrinterMemory* pM = (MyPrinterMemory*) pVoid;
    assert(pM != nullptr);
    assert(pM->pAlignment != nullptr);
    //std::cout << pM->pAlignment->length() << std::endl;
    //std::cout << pM->pAlignment->data.size() << std::endl;
    DEBUG_3(std::cout << c << uiLen << " ";)
    switch (c)
    {
        case 'M':
            for(size_t i = 0; i < uiLen; i++)
            {
                if(
                        (*pM->pQuery)[pM->uiQPos + i] 
                            ==
                        (*pM->pRef)[pM->uiRPos + i]
                    )
                    pM->pAlignment->append(MatchType::match);
                else
                    pM->pAlignment->append(MatchType::missmatch);
            }// for
            pM->uiQPos += uiLen;
            pM->uiRPos += uiLen;
            DEBUG(
                pM->uiQueryExtensionSize += uiLen;
                pM->uiRefExtensionSize += uiLen;
            )
            break;
        case 'I':
            pM->pAlignment->append(MatchType::insertion, uiLen);
            pM->uiQPos += uiLen;
            DEBUG(
                pM->uiQueryExtensionSize += uiLen;
            )
            break;
        case 'D':
            pM->pAlignment->append(MatchType::deletion, uiLen);
            pM->uiRPos += uiLen;
            DEBUG(
                pM->uiRefExtensionSize += uiLen;
            )
            break;
        default:
            std::cout << "GABA cigar contains unknown symbol: " << c << std::endl;
            assert(false);
    }// switch
    return 0;
}// function

void NeedlemanWunsch::dynPrg(
        const std::shared_ptr<NucSeq> pQuery, 
        const std::shared_ptr<NucSeq> pRef,
        const nucSeqIndex fromQuery, const nucSeqIndex toQuery,
        const nucSeqIndex fromRef, const nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment, // in & output
        const bool bLocalBeginning,
        const bool bLocalEnd
    )
{
// use NW naive backend for all cases
#if 0
    //temp
    naiveNeedlemanWunsch(
        pQuery, pRef,
        fromQuery, toQuery,
        fromRef, toRef,
        bLocalBeginning, bLocalEnd,
        pAlignment
    );
    return;
    //temp
#endif

    DEBUG_3(std::cout << "dynProg begin" << std::endl;)
    if(!bLocalBeginning && !bLocalEnd)
    {
        if(
                toQuery - fromQuery > NAIVE_MAX_SIZE ||
                toRef - fromRef > NAIVE_MAX_SIZE ||
                true
            )
            ksw(
                pQuery, pRef,
                fromQuery, toQuery,
                fromRef, toRef,
                pAlignment
            );
        else
            naiveNeedlemanWunsch(
                pQuery, pRef,
                fromQuery, toQuery,
                fromRef, toRef,
                false, false,
                pAlignment
            );
        DEBUG_3(std::cout << "dynProg end" << std::endl;)
        return;
    }// if

    DEBUG_3(std::cout << "sw1" << std::endl;)

    assert(! (bLocalBeginning && bLocalEnd) );
    assert( bLocalBeginning || bLocalEnd );

    const bool bReverse = bLocalBeginning;

    const unsigned int uiBandWidth = 512;

    // in all other cases we use libGaba
    // do some checking for empty sequences though since libGaba does not offer that
    if( toRef <= fromRef )
        if( toQuery <= fromQuery )
        {
            DEBUG_3(std::cout << "dynProg end" << std::endl;)
            return;
        }// if
    if( toQuery <= fromQuery )
    {
        pAlignment->append(MatchType::deletion, toRef - fromRef );
        DEBUG_3(std::cout << "dynProg end" << std::endl;)
        return;
    }//if
    if( toRef <= fromRef )
    {
        pAlignment->append(MatchType::insertion, toQuery - fromQuery );
        DEBUG_3(std::cout << "dynProg end" << std::endl;)
        return;
    }//if

    // if we reached this point we actually have to align something
    DEBUG_2(
        std::cout << pQuery->toString() << std::endl;
        std::cout << pRef->toString() << std::endl;
    )

    DEBUG_3(std::cout << "sw2" << std::endl;)
    /*
    * do the SW alignment
    */
    std::vector<uint8_t> vQuery4Bit = pQuery->as4Bit(fromQuery, toQuery+1, bReverse);
    std::vector<uint8_t>   vRef4Bit =   pRef->as4Bit(  fromRef,   toRef+1, bReverse);
    DEBUG_3(
        if(bLocalBeginning)
        {
            for(auto i : vQuery4Bit)
                std::cout << (int)i; 
            std::cout << std::endl;
            for(auto i : vRef4Bit)
                std::cout << (int)i; 
            std::cout << std::endl;
        }// if
        std::cout << "sw3" << std::endl;
    )// DEBUG_3

    uint8_t const t[ uiBandWidth ] = { 0 }; // tail array

	struct gaba_section_s asec = { // gaba_build_section(0,  &vQuery4Bit[0], vQuery4Bit.size());
        id: 0,
        len: (uint32_t) vRef4Bit.size(),
        base: &vRef4Bit[0]
    };//struct 
	struct gaba_section_s bsec = { // gaba_build_section(2, &vRef4Bit[0], vRef4Bit.size());
        id: 2,
        len: (uint32_t) vQuery4Bit.size(),
        base: &vQuery4Bit[0]
    };//struct 
	struct gaba_section_s tail = { // gaba_build_section(4, t, 64);
        id: 4,
        len: uiBandWidth,
        base: t
    };//struct 
    DEBUG_3(std::cout << "sw4" << std::endl;)

    assert(pGabaScoring->pContext != nullptr);
    
    Gaba_dp_tWrapper xDb( gaba_dp_init(pGabaScoring->pContext) );
    assert(xDb.pDp != nullptr);
    DEBUG_3(std::cout << "sw4.1" << std::endl;)

    //Note f does not need to be freed apparently
	struct gaba_section_s const *ap = &asec, *bp = &bsec;
    
    struct gaba_fill_s const *f = gaba_dp_fill_root(
        xDb.pDp,	    /* dp -> &dp[_dp_ctx_index(band_width)] makes the band width selectable */
		ap, 0,		/* a-side (reference side) sequence and start position */
		bp, 0,		/* b-side (query) */
		UINT32_MAX		/* max extension length */
	);
    DEBUG_3(std::cout << "sw5" << std::endl;)

	/* until X-drop condition is detected */
	struct gaba_fill_s const *m = f;
    /* track max */
	while((f->status & GABA_TERM) == 0)
    {
        /* substitute the pointer by the tail section's if it reached the end */
		if(f->status & GABA_UPDATE_A) { ap = &tail; }
		if(f->status & GABA_UPDATE_B) { bp = &tail; }

        /* extend the banded matrix */
		f = gaba_dp_fill(xDb.pDp, f, ap, bp, UINT32_MAX);
        /* swap if maximum score was updated */
		m = f->max > m->max ? f : m;
	}// while

	xDb.pR = gaba_dp_trace(
        xDb.pDp,
		m,		/* section with the max */
		NULL	/* custom allocator: see struct gaba_alloc_s in gaba.h */
	);//struct
    if(xDb.pR == nullptr)
    {
        DEBUG(
            std::cerr << "WARNING libGaba delivered nullptr on alignment "
                << pAlignment->xStats.sName << std::endl
                << "Maybe the area was too large?" << std::endl
                << "Query: " << fromQuery << " - " << toQuery << std::endl
                << "Ref: " << fromRef << " - " << toRef << std::endl;
        )// DEBUG
        //great what to do other than give up here...?
        if(bReverse)
        {
            pAlignment->shiftOnQuery( toQuery - fromQuery );
            pAlignment->shiftOnRef( toRef - fromRef );
        }// if
        return;
    }// if

    DEBUG_3(std::cout << "sw6" << std::endl;)

    // used int the lambda function
    MyPrinterMemory xPrinter(
        pQuery,
        pRef,
        pAlignment,
        // if we have a reverse alignment libGaba might not have fully aligned till the end
        // so we need to figure out where it ended
        bReverse ? toQuery - (xDb.pR->bgcnt + xDb.pR->dcnt) : fromQuery,
        bReverse ? toRef - (xDb.pR->agcnt + xDb.pR->dcnt) : fromRef
    );

	////printf("score(%" PRId64 "), path length(%" PRIu64 ")\n", xDb.pR->score, xDb.pR->plen);
    /*
     * Having duplicate code here, but need the different functions
     * Maybe template all that?
     */
    if(bReverse)
        gaba_print_cigar_reverse(
            printer,                        /* printer function */
            (void*) &xPrinter,              /* printer function input */
            xDb.pR->path,		            /* bit-encoded path array */
            0,					            /* offset is always zero */
            xDb.pR->plen		            /* path length */
        );
    else
        gaba_print_cigar_forward(
            printer,                        /* printer function */
            (void*) &xPrinter,              /* printer function input */
            xDb.pR->path,		            /* bit-encoded path array */
            0,					            /* offset is always zero */
            xDb.pR->plen		            /* path length */
        );
    DEBUG_3(std::cout << std::endl;)

    assert(xPrinter.uiQueryExtensionSize == xDb.pR->bgcnt + xDb.pR->dcnt);
    assert(xPrinter.uiRefExtensionSize == xDb.pR->agcnt + xDb.pR->dcnt);

    /*
     * Warning:
     * Order is important: the shifting needs to be done after the cigar extraction
     */
    if(bReverse)
    {
        pAlignment->shiftOnRef( vRef4Bit.size() - (xDb.pR->agcnt + xDb.pR->dcnt) );
        pAlignment->shiftOnQuery( vQuery4Bit.size() - (xDb.pR->bgcnt + xDb.pR->dcnt) );
    }// if

    DEBUG_3(std::cout << "dynProg end" << std::endl;)
    return;
}//function

NeedlemanWunsch::NeedlemanWunsch()
#if ALLOCATE_ONCE == 1
        :
    //allocate memory for the naive approach
    s(
        3,
        std::vector<std::vector<int>>(
            NAIVE_MAX_SIZE+1,
            std::vector<int>(NAIVE_MAX_SIZE+1)
        )
    ),
    dir(
        3,
        std::vector<std::vector<char>>(
            NAIVE_MAX_SIZE+1,
            std::vector<char>(NAIVE_MAX_SIZE+1)
        )
    )
#endif
{
    //Gaba context initialization
    /*
     *    GABA_PARAMS(
     *            .xdrop = 100,  
     *            //match award, 
     *            //mismatch penalty, 
     *            //gap open penalty (G_i), and 
     *            //gap extension penalty (G_e) 
     *            GABA_SCORE_SIMPLE(iMatch, iMissMatch, iGap, iExtend)
     *    )
     */
    gaba_params_s xParams = {
        score_matrix: {
                 (int8_t) iMatch, (int8_t) -iMissMatch, (int8_t) -iMissMatch, (int8_t) -iMissMatch,
            (int8_t) -iMissMatch,      (int8_t) iMatch, (int8_t) -iMissMatch, (int8_t) -iMissMatch,
            (int8_t) -iMissMatch, (int8_t) -iMissMatch,      (int8_t) iMatch, (int8_t) -iMissMatch,
            (int8_t) -iMissMatch, (int8_t) -iMissMatch, (int8_t) -iMissMatch,      (int8_t) iMatch
        },
        gi: (int8_t) iGap,
        ge: (int8_t) iExtend,
        gfa: 0,
        gfb: 0,
        xdrop: 100
    };// struct

    pGabaScoring = std::make_shared<Gaba_tWrapper>(xParams);

    //create match/missmatch matrix for ksw
    ksw_gen_simple_mat(5, mat, iMatch, iMissMatch);
}// constructor

std::string NeedlemanWunsch::getFullDesc() const
{
    return std::string("NeedlemanWunsch(") + 
        std::to_string(iMatch) + "," + 
        std::to_string(iMissMatch) + "," + 
        std::to_string(iGap) + "," + 
        std::to_string(iExtend) + "," + 
        std::to_string(uiMaxGapArea) + "," + 
        std::to_string(bLocal) + "," + 
        std::to_string(uiPadding) + ")"
        ;
}//function

ContainerVector NeedlemanWunsch::getInputType() const
{
    return ContainerVector{
        //the sound strip of consideration
        std::shared_ptr<Container>(new Seeds()),
        //the query sequence
        std::shared_ptr<Container>(new NucSeq()),
        //the reference sequence
        std::shared_ptr<Container>(new Pack())
    };
}//function

std::shared_ptr<Container> NeedlemanWunsch::getOutputType() const
{
    return std::shared_ptr<Container>(new Alignment());
}//function


std::shared_ptr<Container> NeedlemanWunsch::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    const auto& pSeeds   = std::static_pointer_cast<Seeds >((*vpInput)[0]);
    const auto& pQuery   = std::static_pointer_cast<NucSeq>((*vpInput)[1]);
    const auto& pRefPack = std::static_pointer_cast<Pack  >((*vpInput)[2]);

    if(pSeeds == nullptr)
        return std::shared_ptr<Alignment>(new Alignment());

    //no seeds => no spot found at all...
    if(pSeeds->empty())
    {
        std::shared_ptr<Alignment> pRet(new Alignment());
        pRet->xStats = pSeeds->xStats;
        pRet->xStats.sName = pQuery->sName;
        return pRet;
    }//if

    DEBUG_2(
        std::cout << "seedlist: (start_ref, end_ref; start_query, end_query)" << std::endl;
        for(Seed& rSeed : *pSeeds)
        {
            std::cout << rSeed.start_ref() << ", " << rSeed.end_ref() << "; "
                << rSeed.start() << ", " << rSeed.end() << std::endl;
        }//for
    )// DEBUG

    // Determine the query and reverence coverage of the seeds

    nucSeqIndex beginRef = pSeeds->front().start_ref();
    nucSeqIndex endRef = pSeeds->back().end_ref();
    //seeds are sorted by ther startpos so we 
    //actually need to check all seeds to get the proper end
    nucSeqIndex endQuery = pSeeds->back().end();
    nucSeqIndex beginQuery = pSeeds->front().start();
    for (auto xSeed : *pSeeds)
    {
        if(endRef < xSeed.end_ref())
            endRef = xSeed.end_ref();
        if(beginRef > xSeed.start_ref())
            beginRef = xSeed.start_ref();
        if(endQuery < xSeed.end())
            endQuery = xSeed.end();
        if(beginQuery > xSeed.end())
            beginQuery = xSeed.start();
        assert(xSeed.start() <= xSeed.end());
    }//for
    DEBUG_2(
        std::cout << beginRef << ", " << endRef << "; " << beginQuery << ", " << endQuery << std::endl;
    )// DEEBUG

    // here we have enough query coverage to attemt to fill in the gaps merely
    DEBUG_2(std::cout << "filling in gaps" << std::endl;)

    std::shared_ptr<Alignment> pRet;

    if(!bLocal)
    {
        beginRef -= uiPadding + beginQuery;
        if(beginRef > endRef)//check for underflow
            beginRef = 0;
        endRef += uiPadding + (pQuery->length() - endQuery);
        if(beginRef > endRef)//check for overflow
            endRef = pRefPack->uiUnpackedSizeForwardPlusReverse();
        endQuery = pQuery->length();
        beginQuery = 0;
    }//if
    assert(endQuery <= pQuery->length());
    pRet = std::shared_ptr<Alignment>(
            new Alignment(beginRef, beginQuery)
        );

    //save the strip of consideration stats in the alignment
    pRet->xStats = pSeeds->xStats;
    pRet->xStats.sName = pQuery->sName;

    DEBUG_2(
        std::cout << beginRef << " " << endRef << std::endl;
    )
    std::shared_ptr<NucSeq> pRef;
    try
    {
        pRef = pRefPack->vExtract(beginRef, endRef);
    } catch(std::runtime_error e)
    {
        std::shared_ptr<Alignment> pRet(new Alignment());
        pRet->xStats = pSeeds->xStats;
        pRet->xStats.sName = pQuery->sName;
        return pRet;
    }// catch

    //create the actual alignment

    if(!bLocal)
    {
        dynPrg(
                pQuery,
                pRef,
                0,
                pSeeds->front().start(),
                0,
                pSeeds->front().start_ref() - beginRef,
                pRet,
                true,
                false
            );
    }// if

    nucSeqIndex endOfLastSeedQuery = pSeeds->front().end();
    nucSeqIndex endOfLastSeedReference = pSeeds->front().end_ref() - beginRef;

    DEBUG(
        if(pRet->uiEndOnQuery != pSeeds->front().start())
        {
            std::cout << pRet->uiEndOnQuery << " ?= " << pSeeds->front().start() 
                << std::endl;
            std::cout << pRet->uiEndOnRef << " ?= " << pSeeds->front().start_ref() 
                << std::endl;
            assert(false);
        }// if
        if(pRet->uiEndOnRef != pSeeds->front().start_ref())
        {
            std::cout << pRet->uiEndOnQuery << " ?= " << pSeeds->front().start() 
                << std::endl;
            std::cout << pRet->uiEndOnRef << " ?= " << pSeeds->front().start_ref() 
                << std::endl;
            assert(false);
        }// if
    )//DEBUG

    pRet->append(MatchType::seed, pSeeds->front().size());
    bool bSkip = true;
    for(Seed& rSeed : *pSeeds)
    {
        // skip the first seed
        // we do this since the seed has already been appended before the loop
        // this makes the loop structure easier since this way 
        // we can always first compute the NW and then append a seed
        if(bSkip)
        {
            bSkip = false;
            continue;
        }//if
        if(rSeed.size() == 0)
            continue;
        nucSeqIndex ovQ = endOfLastSeedQuery - rSeed.start();
        if(rSeed.start() > endOfLastSeedQuery)
            ovQ = 0;
        nucSeqIndex ovR = endOfLastSeedReference - (rSeed.start_ref() - beginRef);
        if(rSeed.start_ref() > endOfLastSeedReference + beginRef)
            ovR = 0;
        nucSeqIndex len = rSeed.size();
        nucSeqIndex overlap = std::max(ovQ, ovR);
        DEBUG_2(
            std::cout << "overlap: " << overlap << std::endl;
        )//DEBUG
        if(len > overlap)
        {
            dynPrg(
                    pQuery,
                    pRef,
                    endOfLastSeedQuery,
                    rSeed.start(),
                    endOfLastSeedReference,
                    rSeed.start_ref() - beginRef,
                    pRet,
                    false,
                    false
                );
            DEBUG(
                //std::cout << pRet->vGapsScatter.size() << std::endl;
                pRet->vGapsScatter.push_back(std::make_pair(
                    rSeed.start() - endOfLastSeedQuery, 
                    rSeed.start_ref() - beginRef - endOfLastSeedReference
                ));
            )//DEBUG
            if(ovQ > ovR)
                pRet->append(MatchType::deletion, ovQ - ovR);
            DEBUG_2(
                for(nucSeqIndex i = ovR; i < ovQ; i++)
                    std::cout << "d";
            )
            if(ovR > ovQ)
                pRet->append(MatchType::insertion, ovR - ovQ);
            DEBUG_2(
                for(nucSeqIndex i = ovQ; i < ovR; i++)
                    std::cout << "i";
                std::cout << std::endl;
            )//DEBUG
            pRet->append(MatchType::seed, len - overlap);
            DEBUG_2(
                std::cout << len - overlap << std::endl;
            )//DEBUG_2
            DEBUG_2(
                for(nucSeqIndex i = overlap; i < len; i++)
                    std::cout << pQuery->charAt(i + rSeed.start());
                std::cout << std::endl;
                for(nucSeqIndex i = overlap; i < len; i++)
                    std::cout << pRef->charAt(i + rSeed.start_ref() - beginRef);
                std::cout << std::endl;
            )//DEBUG
            DEBUG_2(
                for(nucSeqIndex i = 0; i < len - overlap; i++)
                    std::cout << "m";
            )//DEBUG_2
            if(rSeed.end() > endOfLastSeedQuery)
                endOfLastSeedQuery = rSeed.end();
            if(rSeed.end_ref() > endOfLastSeedReference + beginRef)
                endOfLastSeedReference = rSeed.end_ref() - beginRef;
        }//if
    }//for

    if(bLocal)
        assert(std::get<0>(pRet->data.front()) == MatchType::seed);
    assert(std::get<0>(pRet->data.back()) == MatchType::seed);

    DEBUG_2(
        std::cout << std::endl;
    )
    if(bLocal)
        pRet->makeLocal();
    else
    {
        dynPrg(
                pQuery,
                pRef,
                endOfLastSeedQuery,
                endQuery-1,
                endOfLastSeedReference,
                endRef-beginRef-1,
                pRet,
                false,
                true
            );
        //there should never be dangeling deletions with libGaba
        pRet->removeDangeling();
    }//else
    return pRet;
}//function


/* Random nucleotide sequence of length uiLen, represented as codes.
 */
std::vector<char> randomNucSeq( const size_t uiLen )
{
	static const char nucleotides[] = { 0, 1, 2, 3 };

	std::vector<char> vNucSeq( uiLen );
	for (size_t i = 0; i < uiLen; ++i)
	{
		vNucSeq[i] = nucleotides[rand() % (sizeof( nucleotides ) - 1)];
	} // for

	return vNucSeq;
} // function

void testKsw()
{
	/* Seed the random number generator */
	auto uiSeed = (unsigned)time( NULL ); 
	//// std::cout << "SEED:" << uiSeed << std::endl;
	srand(uiSeed);

    size_t uiRefSize = 1;
    size_t uiQuerySize = 10000;

    int8_t iMatch = 1;
    int8_t iMissMatch = 4;
    int8_t iGap = 10;
    int8_t iExtend = 1;
    int iBandwidth = 100;

#if 0 // adaptive bandwidth
    int minAddBandwidth = 10; // must be >= 0 otherwise ksw will not align till the end
    if( std::abs(uiRefSize - uiQuerySize) > iBandwidth + minAddBandwidth)
        iBandwidth = std::abs(uiRefSize - uiQuerySize) + minAddBandwidth;
#endif

    Wrapper_ksw_extz_t ez;

	for (int uiCounter = 0; uiCounter < 100; ++uiCounter)
	{
        std::cout << "score=" << ez.ez->score << " ciglen= " << ez.ez->n_cigar << " cigar=" << std::endl;
        //create match/missmatch matrix
	    int8_t mat[25];
        ksw_gen_simple_mat(5, mat, iMatch, iMissMatch);

	    auto vRefSeq = randomNucSeq( uiRefSize ); // nucleotide sequence (codes)
	    auto vQuerySeq = randomNucSeq( uiQuerySize ); // nucleotide sequence (codes)

        const char vMIDN[] = {'M', 'I', 'D', 'N'};

        std::cout << "dp" << std::endl;
#if 0
		ksw_extz2_sse(
            nullptr, // used in kmalloc so we do not need this
            vQuerySeq.size(), 
            (const uint8_t*)&vQuerySeq[0],
            vRefSeq.size(), 
            (const uint8_t*)&vRefSeq[0], 
            5, //number of input symbols
            mat, // 5x5 matrix match missmatch
            iGap,
            iExtend,
            iBandwidth,
            -1, // break alignment if score drops to fast (-1 == disabled)
            -1, // bonus(?) for reaching the end of the query (-1 == disabled)
            0, // flags (0 == no flags set)
            ez.ez // return value
        );
#else
		ksw_simplified(
            vQuerySeq.size(), 
            (const uint8_t*)&vQuerySeq[0],
            vRefSeq.size(), 
            (const uint8_t*)&vRefSeq[0], 
            iGap,
            iExtend,
            iGap2,
            iExtend2,
            iBandwidth,
            ez.ez // return value
        );
#endif

        std::cout << "score=" << ez.ez->score << " ciglen= " << ez.ez->n_cigar << " cigar=" << std::endl;
        uint32_t qPos = 0;
        uint32_t rPos = 0;
		for (int i = 0; i < ez.ez->n_cigar; ++i)
        {
            uint32_t uiSymb = ez.ez->cigar[i]&0xf;
            uint32_t uiAmount = ez.ez->cigar[i]>>4;
			std::cout << uiAmount <<  vMIDN[uiSymb] << " ";
            if(uiSymb == 0 || uiSymb == 1)
                qPos+=uiAmount;
            if(uiSymb == 0 || uiSymb == 2)
                rPos+=uiAmount;
        }//for
        std::cout << std::endl;
        std::cout << "qPos=" << qPos << " rPos= " << rPos << std::endl;
        assert(qPos == (uint32_t)uiQuerySize);
        assert(rPos == (uint32_t)uiRefSize);
	} // for
}//function

#ifdef WITH_PYTHON
void exportNeedlemanWunsch()
{
    DEBUG(
        //broken debug function
        //boost::python::def("debugNW", &debugNW);

        //test ksw function
        boost::python::def("testKsw", &testKsw);
    )//DEBUG

     //export the segmentation class
    boost::python::class_<
        NeedlemanWunsch,
        boost::python::bases<Module>,
        std::shared_ptr<NeedlemanWunsch>
    >(
        "NeedlemanWunsch"
    )
#if DEBUG_LEVEL >= 1
        .def_readwrite("analyze_heuristics", &bAnalyzeHeuristics)
#endif //DEBUG_LEVEL
    ;
    boost::python::implicitly_convertible<
        std::shared_ptr<NeedlemanWunsch>,
        std::shared_ptr<Module>
    >();

}//function
#endif