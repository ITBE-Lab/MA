/** 
 * @file needlemanWunsch.cpp
 * @author Markus Schmidt
 */
#include "module/needlemanWunsch.h"
#include <bitset>


using namespace libMA;

int iMatch = 2;
int iMissMatch = 4;
int iGap = 6;
int iExtend = 1;
/// @brief the maximal allowed area for a gap between seeds (caps the NW runtime maximum)
//accuracy drops if parameter is set smaller than 10^6
nucSeqIndex uiMaxGapArea = 1000000;
/// @brief the padding on the left and right end of each alignment
nucSeqIndex uiPadding = 100;

/*
 * @todo: fix this problem
 * 
 * arghh this is really ugly...
 * At the moment there is only one sequence that sets all NW and SW parameters correctly:
 * 1) set the parameters.
 * 2) create a new NW module.
 * 3) Then create all other modules that you want to use.....
 */
//the match missmatch matrix
std::shared_ptr<Gaba_tWrapper> pGabaScoring;


/**
 * @todo: this should not be in a .h file
 * @brief the NW dynamic programming algorithm
 * @details
 * This is a naive very slow version,
 * but it computes the correct score and it is capable of doing semi-global alignments.
 *
 * if bNoGapAtBeginning || bNoGapAtEnd :
 *      returns the gap at the beginning or end
 * otherwise : returns 0 
 */
void naiveNeedlemanWunsch(
        std::shared_ptr<NucSeq> pQuery, 
        std::shared_ptr<NucSeq> pRef,
        nucSeqIndex fromQuery,
        nucSeqIndex toQuery,
        nucSeqIndex fromRef,
        nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment
    )
{
    DEBUG_3(std::cout << "naive" << std::endl;)
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
        int iY = toRef-fromRef;
        while(iY > 0)
        {
            pAlignment->append(MatchType::deletion);
            DEBUG_2(
                std::cout << "D";
            )//DEBUG
            iY--;
        }//while
        return;
    }//if
    if(toRef <= fromRef)
    {
        int iX = toQuery-fromQuery;
        while(iX > 0)
        {
            pAlignment->append(MatchType::insertion);
            DEBUG_2(  
                std::cout << "I";
            )//DEBUG
            iX--;
        }//while
        return;
    }//if

    /*
    * beginning of the actual NW
    */
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
    s[2][0][1] = - (iGap + iExtend);
    dir[2][0][1] = DEL | DIR_0_NEXT;
    for(unsigned int x=0; x < 3; x++)
    {
        for(nucSeqIndex uiI = 2; uiI < toQuery-fromQuery+1; uiI++)
        {
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
    }// score for int a and b
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

// @todo: this should not be in a .h file
// @todo: we should make two calls: start in the center and work our way towards both ends
void libMA::dynPrg(
        const std::shared_ptr<NucSeq> pQuery, 
        const std::shared_ptr<NucSeq> pRef,
        const nucSeqIndex fromQuery, const nucSeqIndex toQuery,
        const nucSeqIndex fromRef, const nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment, // in & output
        const bool bLocalBeginning,
        const bool bLocalEnd
    )
{
    DEBUG_3(std::cout << "dynProg begin" << std::endl;)
    if(!bLocalBeginning && !bLocalEnd)
    {
        naiveNeedlemanWunsch(
            pQuery, pRef,
            fromQuery, toQuery,
            fromRef, toRef,
            pAlignment
        );
        DEBUG_3(std::cout << "dynProg end" << std::endl;)
        return;
    }// if

    DEBUG_3(std::cout << "sw1" << std::endl;)

    assert(! (bLocalBeginning && bLocalEnd) );
    assert( bLocalBeginning || bLocalEnd );

    const bool bReverse = bLocalBeginning;

    const unsigned int uiBandWidth = 64;

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

NeedlemanWunsch::NeedlemanWunsch(bool bLocal)
        :
    bLocal(bLocal)
{
    //Gaba context initialization
    /*
     *    GABA_PARAMS(
     *            .xdrop = 100,  // @todo: ?
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
        xdrop: 100 // @todo: ?
    };// struct

    pGabaScoring = std::make_shared<Gaba_tWrapper>(xParams);
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


//@todo currently disabled
#if 0
nucSeqIndex NeedlemanWunsch::needlemanWunsch(
        std::shared_ptr<NucSeq> pQuery, 
        std::shared_ptr<NucSeq> pRef,
        nucSeqIndex fromQuery,
        nucSeqIndex toQuery,
        nucSeqIndex fromRef,
        nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment,
        bool bNoGapAtBeginning,
        bool bNoGapAtEnd
        DEBUG_PARAM(bool bPrintMatrix)
    )
{
    // in these cases parasail does not offer the wanted alignment approach
    // so we use our own slow implementation
    if(bNoGapAtBeginning || bNoGapAtEnd)
        return naiveNeedlemanWunsch(
            pQuery, pRef,
            fromQuery, toQuery,
            fromRef, toRef,
            pAlignment,
            bNoGapAtBeginning, bNoGapAtEnd
            DEBUG_PARAM(bPrintMatrix)
        );
    // in all other cases we use parasail
    // do some checking for empty sequences though since parasail does not offer that
    if(toRef <= fromRef)
        if(toQuery <= fromQuery)
            return 0;
    if(toQuery <= fromQuery)
    {
        pAlignment->append(MatchType::deletion, toRef-fromRef);
        return 0;
    }//if
    if(toRef <= fromRef)
    {
        pAlignment->append(MatchType::insertion, toQuery-fromQuery);
        return 0;
    }//if

    // okay if we reached here we actually have to align something

    /*
     * do the NW alignment
     */

    // Note: parasail does not follow the usual theme where for opening a gap 
   ParsailResultWrapper pResult(parasail_nw_trace_scan_32(
        (const char*)pQuery->pGetSequenceRef() + fromQuery, toQuery - fromQuery,
        (const char*)pRef->pGetSequenceRef() + fromRef, toRef - fromRef,
        iGap + iExtend, iExtend, &matrix
    ));

    //get the cigar
    ParsailCigarWrapper pCigar(parasail_result_get_cigar(
        pResult.get(),
        (const char*)pQuery->pGetSequenceRef() + fromQuery, toQuery - fromQuery,
        (const char*)pRef->pGetSequenceRef() + fromRef, toRef - fromRef,
        &matrix
    ));

    DEBUG_2(
        std::cout << "cigar length: " << pCigar->len << std::endl;
        std::cout << pCigar->beg_query << ", " << pCigar->beg_ref << std::endl;
    )

    //Due to using NW this should hold:
    assert(pCigar->beg_query == 0);
    assert(pCigar->beg_ref == 0);

    /*
     * Decode the cigar
     */

    nucSeqIndex uiQPos = 0;
    nucSeqIndex uiRPos = 0;
    //decode the cigar
    for(int i = 0; i < pCigar->len; i++)
    {
        char c = parasail_cigar_decode_op(pCigar->seq[i]);
        uint32_t uiLen = parasail_cigar_decode_len(pCigar->seq[i]);
        DEBUG_2(std::cout << c << " x" << uiLen << std::endl;)
        switch (c)
        {
            case '=':
                pAlignment->append(MatchType::match, uiLen);
                uiQPos += uiLen;
                uiRPos += uiLen;
                break;
            case 'I':
                pAlignment->append(MatchType::insertion, uiLen);
                uiQPos += uiLen;
                break;
            case 'D':
                pAlignment->append(MatchType::deletion, uiLen);
                uiRPos += uiLen;
                break;
            case 'X':
                pAlignment->append(MatchType::missmatch, uiLen);
                uiQPos += uiLen;
                uiRPos += uiLen;
                break;
            default:
                // there are different CIGAR symbols allowed in the SAM format
                // but parasail should never generate any of them
                assert(false);
                break;
        }//switch
    }//for
    
    //Due to using NW this should hold:
    assert(uiQPos == toQuery - fromQuery);
    assert(uiRPos == toRef - fromRef);
    
    DEBUG_2(
        std::cout << pQuery->length() - uiQPos << ", " << pRef->length() - uiRPos << std::endl;
    )

    //since we use the naive implementation for semi-global alignments we can always return 0 here
    return 0;
}//function
#endif

// Broken debug code for NW
#if 0 
DEBUG(
    /*
     * small debug function so that we can call a raw NW if we want to.
     */
    std::shared_ptr<Alignment> debugNW(
            std::shared_ptr<NucSeq> q, std::shared_ptr<NucSeq> r,
            bool bNoGapAtBeginning, bool bNoGapAtEnd
        )
    {
        auto pAlignment = std::make_shared<Alignment>(0, r->length(), 0, q->length());
        auto uiRes = NeedlemanWunsch(true).needlemanWunsch(
            q,
            r,
            0,
            q->length(),
            0,
            r->length(),
            pAlignment,
            bNoGapAtBeginning,
            bNoGapAtEnd
        );
        if(bNoGapAtBeginning)
            pAlignment->uiBeginOnRef += uiRes;
        if(bNoGapAtEnd)
            pAlignment->uiEndOnRef -= uiRes;
        return pAlignment;
    }//function
)//DEBUG
#endif

//@todo currently disabled
std::shared_ptr<Container> NeedlemanWunsch::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    assert(false);
    return nullptr;
#if 0
    std::shared_ptr<Seeds> pSeeds = std::static_pointer_cast<Seeds>((*vpInput)[0]);
    std::shared_ptr<NucSeq> pQuery 
        = std::static_pointer_cast<NucSeq>((*vpInput)[1]);
    std::shared_ptr<Pack> pRefPack = 
        std::static_pointer_cast<Pack>((*vpInput)[2]);

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

    /*
     * Here we decide weather to actually perform NW in the gaps between seeds or 
     * if we use SW to align the entire thing.
     */
    if(endQuery - beginQuery < pQuery->length() * fMinimalQueryCoverage)
    {
        DEBUG_2(std::cout << "computing SW for entire area" << std::endl;)
        //figure out the correct reference location
        nucSeqIndex refStart = beginRef;
        if(refStart >= uiPadding)
            refStart -=  uiPadding;
        else
            refStart = 0;
        nucSeqIndex refEnd = endRef;
        assert(endQuery <= pQuery->length());
        refEnd += uiPadding;

        assert(refEnd > refStart);
        nucSeqIndex refWidth = refEnd - refStart;
        if(refStart + refWidth >= pRefPack->uiUnpackedSizeForwardPlusReverse())
            refWidth = pRefPack->uiUnpackedSizeForwardPlusReverse() - refStart - 1;
        assert(refWidth > 0);
        if(pRefPack->bridgingSubsection(refStart, refWidth))
        {
            DEBUG_2(std::cout << "Un-bridging from " << refStart << ", " << refWidth;)
            pRefPack->unBridgeSubsection(refStart, refWidth);
            DEBUG_2(std::cout << " to: " << refStart << ", " << refWidth << std::endl;)
        }
        //extract the reference
        std::shared_ptr<NucSeq> pRef = pRefPack->vExtract(
            refStart,
            refStart + refWidth
        );
        //compute the SW alignment
        auto pRet = smithWaterman(pQuery, pRef, refStart);
        //copy the stats
        pRet->xStats = pSeeds->xStats;
        pRet->xStats.sName = pQuery->sName;
        // return the SW alignment
        return pRet;
    }//if

    // here we have enough query coverage to attemt to fill in the gaps merely
    DEBUG_2(std::cout << "filling in gaps" << std::endl;)

    std::shared_ptr<Alignment> pRet;

    if(!bLocal)
    {
        beginRef -= uiPadding;
        if(beginRef > endRef)//check for underflow
            beginRef = 0;
        endRef += uiPadding;
        if(beginRef > endRef)//check for overflow
            endRef = pRefPack->uiUnpackedSizeForwardPlusReverse()-1;
        endQuery = pQuery->length();
        beginQuery = 0;
    }//if
    assert(endQuery <= pQuery->length());
    pRet = std::shared_ptr<Alignment>(
        new Alignment(beginRef, endRef, beginQuery, endQuery)
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
    }

    //create the actual alignment

    if(!bLocal)
    {
        pRet->uiBeginOnRef += needlemanWunsch(
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
    }//else

    nucSeqIndex endOfLastSeedQuery = pSeeds->front().end();
    nucSeqIndex endOfLastSeedReference = pSeeds->front().end_ref() - beginRef;

    pRet->append(MatchType::seed, pSeeds->front().size());
    bool bSkip = true;
    for(Seed& rSeed : *pSeeds)
    {
        //skip the first seed
        if(bSkip)
        {
            bSkip = false;
            continue;
        }//if
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
            needlemanWunsch(
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
        pRet->uiEndOnRef -= needlemanWunsch(
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
        pRet->removeDangeling();
    }//else


    return pRet;

#endif
}//function


std::string LocalToGlobal::getFullDesc() const
{
    return std::string("LocalToGlobal(") + 
        std::to_string(fMappingQualMin) + ")"
        ;
}//function

ContainerVector LocalToGlobal::getInputType() const
{
    return ContainerVector{
        //the local alignment
        std::make_shared<ContainerVector>( std::make_shared<Alignment>() ),
        //the query sequence
        std::shared_ptr<Container>(new NucSeq()),
        //the reference sequence
        std::shared_ptr<Container>(new Pack())
    };
}//function

std::shared_ptr<Container> LocalToGlobal::getOutputType() const
{
    return std::make_shared<ContainerVector>( std::make_shared<Alignment>() );
}//function

//@todo currently disabled
std::shared_ptr<Container> LocalToGlobal::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    assert(false);
    return nullptr;
#if 0
    const auto& pAlignments = std::static_pointer_cast<ContainerVector>((*vpInput)[0]);
    const std::shared_ptr<NucSeq>& pQuery = std::static_pointer_cast<NucSeq>((*vpInput)[1]);
    const std::shared_ptr<Pack>& pRefPack = std::static_pointer_cast<Pack>((*vpInput)[2]);

    if(pAlignments->size() < 2)
        return pAlignments;

    // check if the mapping quality is lower than the threshold
    auto pFirst = std::static_pointer_cast<Alignment>((*pAlignments)[0]);
    auto pSecond = std::static_pointer_cast<Alignment>((*pAlignments)[1]);
    float fQual =
                ( pFirst->score() - pSecond->score() )
                    /
                (double)( pFirst->score() )
            ;

    if(fQual > fMappingQualMin)
        return std::make_shared<ContainerVector>( pAlignments );

    auto pRet = std::make_shared<ContainerVector>( std::make_shared<Alignment>() );

    for( const auto& pElement : *pAlignments )
    {
        const auto& pAlignment = std::static_pointer_cast<Alignment>(pElement);

        nucSeqIndex beginRef =
            pAlignment->uiBeginOnRef - uiPadding;
        //check for underflow
        if(uiPadding > pAlignment->uiBeginOnRef)
            beginRef = 0;
        DEBUG(
            if(pAlignment->uiBeginOnRef < beginRef)
            {
                std::cerr << pAlignment->uiBeginOnRef << " " << beginRef << " " << pAlignment->uiBeginOnQuery << " " << uiPadding << std::endl;
                assert(false);
            }// if
        )// DEBUG
        nucSeqIndex endRef =
            pAlignment->uiEndOnRef + uiPadding;

        std::shared_ptr<Alignment> pAppend(new Alignment(beginRef, endRef, 0, pQuery->length()));
        pAppend->xStats = pAlignment->xStats;

        std::shared_ptr<NucSeq> pRef;
        try
        {
            pRef = pRefPack->vExtract(beginRef, endRef);
        } catch(std::runtime_error e)
        {
            pRet->push_back(pAlignment);
            continue;
        }// catch

        // limit the maximal area that NW computes; 0 means unlimited
        if(uiMaxGapArea != 0 && pAlignment->uiBeginOnQuery * (pAlignment->uiBeginOnRef - beginRef) > uiMaxGapArea )
        {
            pAppend->uiBeginOnRef = pAlignment->uiBeginOnRef;
            pAppend->uiBeginOnQuery = pAlignment->uiBeginOnQuery;
        }//if
        else
            pAppend->uiBeginOnRef += naiveNeedlemanWunsch(
                pQuery,
                pRef,
                0,
                pAlignment->uiBeginOnQuery,
                0,
                pAlignment->uiBeginOnRef - beginRef,
                pAppend,
                true,
                false
            );

        pAppend->append(*pAlignment);

        // limit the maximal area that NW computes; 0 means unlimited
        if(uiMaxGapArea != 0 && (pQuery->length() - pAlignment->uiEndOnQuery) * 
            (pRef->length() - pAlignment->uiEndOnRef - beginRef) > uiMaxGapArea )
        {
            pAppend->uiEndOnQuery = pAlignment->uiEndOnQuery;
            pAppend->uiEndOnRef = pAlignment->uiEndOnRef;
        }//if
        else
            pAppend->uiEndOnRef -= naiveNeedlemanWunsch(
                pQuery,
                pRef,
                pAlignment->uiEndOnQuery,
                pQuery->length(),
                pAlignment->uiEndOnRef - beginRef,
                pRef->length(),
                pAppend,
                false,
                true
            );

        pRet->push_back(pAppend);
    }//for

    //sort ascending
    std::sort(
        pRet->begin(), pRet->end(),
        []
        (std::shared_ptr<Container> a, std::shared_ptr<Container> b)
        {
            return a->larger(b);
        }//lambda
    );//sort function call
    assert(pRet->size() <= 1 || !pRet->back()->larger(pRet->front()));

    return pRet;
#endif
}//function



ContainerVector CombatRepetitively::getInputType() const
{
    return ContainerVector{
        //the local alignment
        std::make_shared<ContainerVector>( std::make_shared<Alignment>() ),
        //the query sequence
        std::shared_ptr<Container>(new NucSeq()),
        //the reference sequence
        std::shared_ptr<Container>(new Pack())
    };
}//function

std::shared_ptr<Container> CombatRepetitively::getOutputType() const
{
    return std::make_shared<ContainerVector>( std::make_shared<Alignment>() );
}//function


//DEPRECATED
std::shared_ptr<Container> CombatRepetitively::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    assert(false);
    return nullptr;
#if 0
    const auto& pAlignments = std::static_pointer_cast<ContainerVector>((*vpInput)[0]);
    const std::shared_ptr<NucSeq>& pQuery = std::static_pointer_cast<NucSeq>((*vpInput)[1]);
    const std::shared_ptr<Pack>& pRefPack = std::static_pointer_cast<Pack>((*vpInput)[2]);

    if(pAlignments->size() < 2)
        return pAlignments;

    // check if the mapping quality is lower than the threshold
    auto pFirst = std::static_pointer_cast<Alignment>((*pAlignments)[0]);
    auto pSecond = std::static_pointer_cast<Alignment>((*pAlignments)[1]);
    float fQual =
                ( pFirst->score() - pSecond->score() )
                    /
                (double)( pFirst->score() )
            ;

    if(fQual > fMappingQualMax)
        return pAlignments;

    // if we reach this point we have to go through all alignments and check the are that they cover
    nucSeqIndex refStart = pFirst->uiBeginOnRef;
    nucSeqIndex refEnd = pFirst->uiEndOnRef;
    for( const auto& pElement : *pAlignments )
    {
        const auto& pAlignment = std::static_pointer_cast<Alignment>(pElement);
        if(refStart > pAlignment->uiBeginOnRef)
            refStart = pAlignment->uiBeginOnRef;
        if(refEnd < pAlignment->uiEndOnRef)
            refEnd = pAlignment->uiEndOnRef;
    }// for

    if(refEnd - refStart > uiRegionLength)
        return pAlignments;

    //add padding to the reference and make sure we do not extract something bridging
    if(refStart >= uiPadding)
        refStart -= uiPadding;
    else
        refStart = 0;
    refEnd += uiPadding;

    assert(refEnd > refStart);
    nucSeqIndex refWidth = refEnd - refStart;
    if(refStart + refWidth >= pRefPack->uiUnpackedSizeForwardPlusReverse())
        refWidth = pRefPack->uiUnpackedSizeForwardPlusReverse() - refStart - 1;
    assert(refWidth > 0);
    if(pRefPack->bridgingSubsection(refStart, refWidth))
    {
        DEBUG_2(std::cout << "Un-bridging from " << refStart << ", " << refWidth;)
        pRefPack->unBridgeSubsection(refStart, refWidth);
        DEBUG_2(std::cout << " to: " << refStart << ", " << refWidth << std::endl;)
    }
    //extract the reference
    std::shared_ptr<NucSeq> pRef = pRefPack->vExtract(
        refStart,
        refStart + refWidth
    );

    auto pAlign = smithWaterman(pQuery, pRef, refStart);
    //copy the stats
    pAlign->xStats = pFirst->xStats;
    pAlignments->push_back(pAlign);

    //sort alignments
    std::sort(
        pAlignments->begin(), pAlignments->end(),
        []
        (std::shared_ptr<Container> a, std::shared_ptr<Container> b)
        {
            return a->larger(b);
        }//lambda
    );//sort function call

    return pAlignments;
#endif
}// function

void exportNeedlemanWunsch()
{
    DEBUG(
        //broken debug function
        //boost::python::def("debugNW", &debugNW);
    )//DEBUG
     //export the segmentation class
    boost::python::class_<
        NeedlemanWunsch,
        boost::python::bases<Module>,
        std::shared_ptr<NeedlemanWunsch>
    >(
        "NeedlemanWunsch",
        boost::python::init<bool>()
    )
        // These are constants among the entire code...
        // We set them using the NW class for simplicity
        // @todo: this should be changed ...
        .def_readwrite("penalty_gap_open", &iGap)
        .def_readwrite("penalty_gap_extend", &iExtend)
        .def_readwrite("score_match", &iMatch)
        .def_readwrite("penalty_missmatch", &iMissMatch)
        .def_readwrite("max_gap_Area", &uiMaxGapArea)
        // actual parameters of NW
        .def_readwrite("local", &NeedlemanWunsch::bLocal)
        .def_readwrite("min_coverage", &NeedlemanWunsch::fMinimalQueryCoverage)
    ;
    boost::python::implicitly_convertible<
        std::shared_ptr<NeedlemanWunsch>,
        std::shared_ptr<Module>
    >();

     //export the LocalToGlobal class
    boost::python::class_<
        LocalToGlobal,
        boost::python::bases<Module>,
        std::shared_ptr<LocalToGlobal>
    >(
        "LocalToGlobal",
        boost::python::init<double>()
    );
    boost::python::implicitly_convertible<
        std::shared_ptr<LocalToGlobal>,
        std::shared_ptr<Module>
    >();

     //export the CombatRepetitively class
    boost::python::class_<
        CombatRepetitively,
        boost::python::bases<Module>,
        std::shared_ptr<CombatRepetitively>
    >(
        "CombatRepetitively",
        boost::python::init<double, nucSeqIndex>()
    );
    boost::python::implicitly_convertible<
        std::shared_ptr<CombatRepetitively>,
        std::shared_ptr<Module>
    >();

}//function