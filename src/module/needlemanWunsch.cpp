#include "module/needlemanWunsch.h"
using namespace libMA;


int iGap = 6;
int iExtend = 1;
int iMatch = 1;
int iMissMatch = 4;

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
    /*
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
    std::vector<std::vector<int>> s(toQuery-fromQuery+1, std::vector<int>(toRef-fromRef+1));
    std::vector<std::vector<char>> dir(toQuery-fromQuery+1, std::vector<char>(toRef-fromRef+1));

    /*
     * initialization:
     *      this part sets the scores for the last row and column (reverse order)
     * 
     * Note:
     *      if we do not want a gap at the end since the alignment ends there we need to 
     *      set the initial values along the reference to 0.
     *      we do not want a complete global alignment,
     *      merely a global alignment with respect to the query
     * 
     * @todo this can happen:

CATTACTTTATAGATTGGGAACAATCCCATTCAAAGT-------------------------------------------        reference
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII|~~
CATTACTTTATAGATTGGGAACAATCCCATTCAAAAAAGAGCGCTTCATCTTAACTTAGGGGTAGGTCCATTAGATAGCC        query

800-880
---------------------------------------------------------AAAGAAGGACTTGGCATCTGCCA        reference
                                                       ~~IIIIIIIIIIIIIIIIIIIIIII
CAATCGGACCTATACATGGGGAGCTATATTTTATATACTCGCCCACCAATGGAGTGTAAAGAAGGACTTGGCATCTGCCA        query

     */
    s[0][0] = 0;
    dir[0][0] = 1;
    s[1][0] = - (iGap + iExtend);
    dir[1][0] = 2;
    if(bNoGapAtEnd)//see note above
        s[0][1] = 0;
    else
        s[0][1] = - (iGap + iExtend);
    dir[0][1] = 3;
    for(nucSeqIndex uiI = 2; uiI < toQuery-fromQuery+1; uiI++)
    {
        s[uiI][0] = s[uiI - 1][0] - iExtend;
        dir[uiI][0] = 2;
    }//for
    for(nucSeqIndex uiI = 2; uiI < toRef-fromRef+1; uiI++)
    {
        if(bNoGapAtEnd)//see note above
            s[0][uiI] = 0;
        else
            s[0][uiI] = s[0][uiI - 1] - iExtend;
        dir[0][uiI] = 3;
    }//for
    /*
     * dynamic programming loop
     * Note:
     *      we iterate in the reverse order on reference and query
     *      so that the backtracking can be done in forward order
     *
     * This works as follows:
     *      for each cell compute the scores if resuling from an insertion deletion match/missmatch
     *      in this order. Store the score from the insertion and overwrite the score with the del
     *      match of missmatch score if any of them is higher. Also keep track of which direction
     *      we came from in the dir matrix.
     */
    for(nucSeqIndex uiI = 1; uiI < (toQuery-fromQuery)+1; uiI++)
    {
        for(nucSeqIndex uiJ = 1; uiJ < (toRef-fromRef)+1; uiJ++)
        {
            int newScore;
            //insertion
            if(dir[uiI - 1][uiJ] == 2)
                newScore = s[uiI - 1][uiJ] - iExtend;
            else
                newScore = s[uiI - 1][uiJ] - (iGap + iExtend);
            s[uiI][uiJ] = newScore;
            dir[uiI][uiJ] = 2;

            //deletion
            if(dir[uiI][uiJ - 1] == 3)
                newScore = s[uiI][uiJ - 1] - iExtend;
            else
                newScore = s[uiI][uiJ - 1] - (iGap + iExtend);
            if(newScore > s[uiI][uiJ])
            {
                s[uiI][uiJ] = newScore;
                dir[uiI][uiJ] = 3;
            }//if
            //match / missmatch
            newScore = s[uiI - 1][uiJ - 1];
            //@todo try -1s here and see what happens
            if( (*pQuery)[toQuery - uiI] == (*pRef)[toRef - uiJ] )
                newScore += iMatch;
            else
                newScore -= iMissMatch;
            if(newScore >= s[uiI][uiJ])
            {
                s[uiI][uiJ] = newScore;
                dir[uiI][uiJ] = 1;
            }//if
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
    
    /*
    * if there is no gap cost for the beginning 
    * we should start backtracking where the score is maximal
    * along the reference
    */
    if(bNoGapAtBeginning)
    {
        for(nucSeqIndex uiJ = 1; uiJ < toRef-fromRef; uiJ++)
            if(s[iX][uiJ] > s[iX][iY])
                iY = uiJ;
        uiRet = (toRef-fromRef) - iY;
        DEBUG_2(
            std::cout << (toRef-fromRef) - iY << "D";
        )//DEBUG
    }//if

    while(iX > 0 || iY > 0)
    {
        if(dir[iX][iY] == 1)
        {
            if( (*pQuery)[toQuery - iX] == (*pRef)[toRef - iY] )
            {
                pAlignment->append(MatchType::match);
                DEBUG_2(
                    std::cout << "M";
                )//DEBUG
            }
            else
            {
                pAlignment->append(MatchType::missmatch);
                DEBUG_2(
                    std::cout << "W";
                )//DEBUG
            }
            iX--;
            iY--;
        }//if
        else if(dir[iX][iY] == 2)
        {
            pAlignment->append(MatchType::insertion);
            iX--;
            DEBUG_2(
                std::cout << "I";
            )//DEBUG
        }//if
        else if(dir[iX][iY] == 3)
        {
            pAlignment->append(MatchType::deletion);
            iY--;
            DEBUG_2(
                std::cout << "D";
            )//DEBUG        
        }//if
        else{
            std::cerr << "WARNING: no direction set dynamic programming" << std::endl;
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
    return uiRet;
}//function

std::shared_ptr<Container> NeedlemanWunsch::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<Seeds> pSeeds = std::static_pointer_cast<Seeds>((*vpInput)[0]);
    std::shared_ptr<NucSeq> pQuery 
        = std::static_pointer_cast<NucSeq>((*vpInput)[1]);
    std::shared_ptr<Pack> pRefPack = 
        std::static_pointer_cast<Pack>((*vpInput)[2]);

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
    }//for
    if(!bLocal)
    {
        beginRef -= (nucSeqIndex)( beginQuery * fRelativePadding );
        if(beginRef > endRef)//check for underflow
            beginRef = 0;
        assert(pQuery->length() >= endQuery);
        endRef += (nucSeqIndex)( (pQuery->length() - endQuery) * fRelativePadding );
        if(beginRef > endRef)//check for overflow
            endRef = pRefPack->uiUnpackedSizeForwardPlusReverse()-1;
        endQuery = pQuery->length();
        beginQuery = 0;
    }//if
    pRet = std::shared_ptr<Alignment>(
        new Alignment(beginRef, endRef, beginQuery, endQuery)
    );

    //save the strip of consideration stats in the alignment
    pRet->xStats = pSeeds->xStats;
    pRet->xStats.sName = pQuery->sName;

    DEBUG(
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
        DEBUG(
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
            DEBUG(
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