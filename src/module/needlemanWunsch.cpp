/** 
 * @file needlemanWunsch.cpp
 * @author Markus Schmidt
 */
#include "module/needlemanWunsch.h"
#include <bitset>

using namespace libMA;

int iMatch = 10;
int iMissMatch = 4;
int iGap = 6;
int iExtend = 1;


std::string NeedlemanWunsch::getFullDesc() const
{
    return std::string("NeedlemanWunsch(") + 
        std::to_string(iMatch) + "," + 
        std::to_string(iMissMatch) + "," + 
        std::to_string(iGap) + "," + 
        std::to_string(iExtend) + "," + 
        std::to_string(bLocal) + "," + 
        std::to_string(fRelativePadding) + ")"
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

/*
 * the NW dynamic programming algorithm
 * 
 * if bNoGapAtBeginning || bNoGapAtEnd
 *  returns the gap at the beginning or end
 * returns 0 otherwise
 */
nucSeqIndex needlemanWunsch(
        std::shared_ptr<NucSeq> pQuery, 
        std::shared_ptr<NucSeq> pRef,
        nucSeqIndex fromQuery,
        nucSeqIndex toQuery,
        nucSeqIndex fromRef,
        nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment,
        bool bNoGapAtBeginning = false,
        bool bNoGapAtEnd = false
        DEBUG_PARAM(bool bPrintMatrix = false)
    )
{
    assert(!(bNoGapAtBeginning && bNoGapAtEnd));
    /*
     * break conditions for actually empty areas
     */
    assert(toQuery <= pQuery->length());
    assert(toRef <= pRef->length());
    if(toRef <= fromRef)
        if(toQuery <= fromQuery)
            return 0;
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
        return 0;
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
        return 0;
    }//if
#if 0//DEPRECATED
    /**
     * give up for too large areas
     * @todo this should split everything into two local alignments instead
     */
    if(uiGiveUpAfter != 0)
    {
        nucSeqIndex uiArea = (toQuery - fromQuery)*(toRef - fromRef);
        if(uiArea > uiGiveUpAfter)
        {
            int diagLen = std::min(toRef - fromRef, toQuery - fromQuery);
            int gapLen = std::max(toRef - fromRef, toQuery - fromQuery) - diagLen;

            int scoreMissmatch = -1* (iGap + iExtend*gapLen + 
                                 iMissMatch*diagLen);

            int scoreIndelOnly = -1* (iGap*2 + iExtend*(toRef - fromRef) + 
                                 iExtend*(toQuery - fromQuery));

            if(scoreIndelOnly < scoreMissmatch)
            {
                if(toRef - fromRef > toQuery - fromQuery)
                {
                    pAlignment->append(MatchType::deletion, gapLen);
                    pAlignment->append(MatchType::missmatch, diagLen);
                }//if
                else
                {
                    pAlignment->append(MatchType::insertion, gapLen);
                    pAlignment->append(MatchType::missmatch, diagLen);
                }//else
            }//if
            else
            {
                pAlignment->append(MatchType::deletion, toRef - fromRef);
                pAlignment->append(MatchType::insertion, toQuery - fromQuery);
            }//else
            return;
        }//if
    }//if
#endif

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
    if(bNoGapAtEnd)//see note above
        s[2][0][1] = 0;
    else
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
            if(bNoGapAtEnd)//see note above
                s[x][0][uiI] = 0;
            else
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
                std::cout << s[uiI][uiJ] << "\t";
            std::cout << std::endl;
        }//for
    )//DEBUG

    nucSeqIndex uiRet = 0;

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
        uiRet = (toRef-fromRef) - iY;
        DEBUG_2(
            std::cout << (toRef-fromRef) - iY << "D";
        )//DEBUG
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
    }//else
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
            return iY;
    }//while

    DEBUG_2(
        std::cout << std::endl;
    )//DEBUG

    //print the entire matrix if necessary
    DEBUG(
        if(bPrintMatrix)
        {
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
        }//if
    )//DEBUG
    return uiRet;
}//function

DEBUG(
    void debugNW(std::shared_ptr<NucSeq> q, std::shared_ptr<NucSeq> r)
    {
        auto pAlignment = std::make_shared<Alignment>();
        needlemanWunsch(
            q,
            r,
            0,
            q->length(),
            0,
            r->length(),
            pAlignment,
            false,
            false,
            true
        );
    }//function
)//DEBUG

std::shared_ptr<Container> NeedlemanWunsch::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
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
    )


    std::shared_ptr<Alignment> pRet;
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
        assert(xSeed.start() <= xSeed.end());
    }//for
    if(!bLocal)
    {
        beginRef -= (nucSeqIndex)( beginQuery * fRelativePadding );
        if(beginRef > endRef)//check for underflow
            beginRef = 0;
        endRef += (nucSeqIndex)( (pQuery->length() - endQuery) * fRelativePadding );
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
            true
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
                    pRet
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
    }//else

    return pRet;

}//function

void exportNeedlemanWunsch()
{
    DEBUG(
        boost::python::def("debugNW", &debugNW);
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
        .def_readwrite("penalty_gap_open", &iGap)
        .def_readwrite("penalty_gap_extend", &iExtend)
        .def_readwrite("score_match", &iMatch)
        .def_readwrite("penalty_missmatch", &iMissMatch)
        .def_readwrite("local", &NeedlemanWunsch::bLocal)
        .def_readwrite("relative_padding", &NeedlemanWunsch::fRelativePadding)
    ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<NeedlemanWunsch>,
        std::shared_ptr<Module> 
    >();

}//function