#include "module/needlemanWunsch.h"
using namespace libMABS;


int iGap = 50;//20
int iExtend = 1;
int iMatch = 20;//2
int iMissMatch = 20;//2

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



void needlemanWunsch(
        std::shared_ptr<NucSeq> pQuery, 
        std::shared_ptr<NucSeq> pRef,
        nucSeqIndex fromQuery,
        nucSeqIndex toQuery,
        nucSeqIndex fromRef,
        nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment,
        nucSeqIndex uiGiveUpAfter
    )
{
    assert(toQuery <= pQuery->length());
    assert(toRef <= pRef->length());
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
//switch banded (not working yet) (1) non-banded (0)
#if 0
    nucSeqIndex uiBandSize = std::min(200, toQuery-fromQuery+1, toRef-fromRef+1);
    nucSeqIndex uiBandCenter = uiBandSize/2;
    nucSeqIndex uiBandLength = std::max(toQuery-fromQuery+1, toRef-fromRef+1);

    std::vector<std::vector<int>> s(uiBandLength, std::vector<int>(uiBandSize));
    std::vector<std::vector<char>> dir(uiBandLength, std::vector<char>(uiBandSize));

    s[0][0] = 0;
    dir[0][0] = 1;
    s[1][0] = -iGap;
    dir[1][0] = 2;
    s[0][1] = -iGap;
    dir[0][1] = 3;

#else
    std::vector<std::vector<int>> s(toQuery-fromQuery+1, std::vector<int>(toRef-fromRef+1));
    std::vector<std::vector<char>> dir(toQuery-fromQuery+1, std::vector<char>(toRef-fromRef+1));

    s[0][0] = 0;
    dir[0][0] = 1;
    s[1][0] = -iGap;
    dir[1][0] = 2;
    s[0][1] = -iGap;
    dir[0][1] = 3;
    for(nucSeqIndex uiI = 2; uiI < toQuery-fromQuery+1; uiI++)
    {
        s[uiI][0] = s[uiI - 1][0] - iExtend;
        dir[uiI][0] = 2;
    }//for
    for(nucSeqIndex uiI = 2; uiI < toRef-fromRef+1; uiI++)
    {
        s[0][uiI] = s[0][uiI - 1] - iExtend;
        dir[0][uiI] = 3;
    }//for
    for(nucSeqIndex uiI = 1; uiI < toQuery-fromQuery+1; uiI++)
    {
        for(nucSeqIndex uiJ = 1; uiJ < toRef-fromRef+1; uiJ++)
        {
            int newScore;
            //insertion
            if(dir[uiI - 1][uiJ] == 2)
                newScore = s[uiI - 1][uiJ] - iExtend;
            else
                newScore = s[uiI - 1][uiJ] - iGap;
            s[uiI][uiJ] = newScore;
            dir[uiI][uiJ] = 2;

            //deletion
            if(dir[uiI][uiJ - 1] == 3)
                newScore = s[uiI][uiJ - 1] - iExtend;
            else
                newScore = s[uiI][uiJ - 1] - iGap;
            // for the first alignment we dont want to have a malus for an
            // deletion of the reference at the beginning
            if(fromQuery == 0 && uiI == toQuery-fromQuery)
                newScore = s[uiI][uiJ - 1];
            if(newScore > s[uiI][uiJ])
            {
                s[uiI][uiJ] = newScore;
                dir[uiI][uiJ] = 3;
            }//if
            //match / missmatch
            newScore = s[uiI - 1][uiJ - 1];
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

    nucSeqIndex iX = toQuery-fromQuery;
    nucSeqIndex iY = toRef-fromRef;
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
        else if(dir[iX][iY] == 3)
        {
            pAlignment->append(MatchType::deletion);
            iY--;
            DEBUG_2(
                std::cout << "D";
            )//DEBUG        
        }//if
        else
        {
            pAlignment->append(MatchType::insertion);
            iX--;
            DEBUG_2(
                std::cout << "I";
            )//DEBUG
        }//if
    }//while
    DEBUG_2(
        std::cout << std::endl;
    )//DEBUG

#endif


}//function

std::shared_ptr<Container> NeedlemanWunsch::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    //switch local or global alignment
    bool bLocal = true;

    std::shared_ptr<Seeds> pSeeds = std::static_pointer_cast<Seeds>((*vpInput)[0]);
    std::shared_ptr<NucSeq> pQuery 
        = std::static_pointer_cast<NucSeq>((*vpInput)[1]);
    std::shared_ptr<Pack> pRefPack = 
        std::static_pointer_cast<Pack>((*vpInput)[2]);

    //no seeds => no spot found at all...
    if(pSeeds->empty())
    {
        std::shared_ptr<Alignment> pRet(new Alignment());
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
    nucSeqIndex beginRef;
    nucSeqIndex endRef;
    if(bLocal)
    {
        beginRef = pSeeds->front().start_ref();
        endRef = pSeeds->back().end_ref();
        pRet = std::shared_ptr<Alignment>(
            new Alignment(beginRef, endRef, pSeeds->front().start(), pSeeds->back().end())
        );
    }//if
    else
    {
        float fDistFac = 1.0;
        nucSeqIndex beginQuery = pSeeds->front().start();
        beginRef = 0;
        if( pSeeds->front().start_ref() > (nucSeqIndex)(beginQuery*fDistFac))
            beginRef = pSeeds->front().start_ref() - (nucSeqIndex)(beginQuery*fDistFac);
        nucSeqIndex endQuery = pSeeds->back().end();
        assert(endQuery <= pQuery->length());
        endRef = pRefPack->uiUnpackedSizeForwardPlusReverse();
        if( 
                pSeeds->back().end_ref() + (nucSeqIndex)((pQuery->length()-endQuery)*fDistFac) <
                pRefPack->uiUnpackedSizeForwardPlusReverse()
            )
            endRef = pSeeds->back().end_ref() + (nucSeqIndex)((pQuery->length()-endQuery)*fDistFac);
        pRet = std::shared_ptr<Alignment>(
            new Alignment(beginRef, endRef, 0, 0)
        );
    }//else

    //save the strip of consideration stats in the alignment
    pRet->xStats = pSeeds->xStats;

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
        return pRet;
    }

    //create the actual alignment
    nucSeqIndex endOfLastSeedQuery = 0;
    nucSeqIndex endOfLastSeedReference = 0;

    if(bLocal)
    {
        endOfLastSeedQuery = pSeeds->front().start();
        //                          pSeeds->front().start_ref() - beginRef == 0
        //endOfLastSeedReference = pSeeds->front().start_ref() - beginRef;
    }//if

    for(Seed& rSeed : *pSeeds)
    {
        needlemanWunsch(
                pQuery,
                pRef,
                endOfLastSeedQuery,
                rSeed.start(),
                endOfLastSeedReference,
                rSeed.start_ref() - beginRef,
                pRet,
                uiGiveUpAfter
            );
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
        if(len > overlap)
        {
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
        }//if
        if(rSeed.end() > endOfLastSeedQuery)
            endOfLastSeedQuery = rSeed.end();
        if(rSeed.end_ref() > endOfLastSeedReference + beginRef)
            endOfLastSeedReference = rSeed.end_ref() - beginRef;
    }//for

    if(!bLocal)
    {
        needlemanWunsch(
            pQuery,
            pRef,
            endOfLastSeedQuery,
            pQuery->length(),
            endOfLastSeedReference,
            endRef - beginRef,
            pRet,
            uiGiveUpAfter
        );
    }//if

    //cleanup the alignment
    pRet->removeDangelingDeletions();
    // @todo  due to NMW changes i might have a lot of dangeling insertions here...
    // maybe best to keep them around?

    DEBUG_2(
        std::cout << std::endl;
    )
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
        boost::python::init<nucSeqIndex>()
    )
        .def_readwrite("give_up_after", &NeedlemanWunsch::uiGiveUpAfter)
    ;
    boost::python::implicitly_convertible< 
        std::shared_ptr<NeedlemanWunsch>,
        std::shared_ptr<Module> 
    >();

    //boost::python::scope().attr("score_gap_open") = &iGap;
    //boost::python::scope().attr("score_gap_extend") = &iExtend;
    //boost::python::scope().attr("score_match") = &iMatch;
    //boost::python::scope().attr("score_missmatch") = &iMissMatch;

}//function