#include "needlemanWunsch.h"
using namespace libLAuS;



ContainerVector NeedlemanWunsch::getInputType() const
{
    return ContainerVector{
        //the sound strip of consideration
        std::shared_ptr<Container>(new Seeds()),
        //the query sequence
        std::shared_ptr<Container>(new NucSeq()),
        //the reference sequence
        std::shared_ptr<Container>(new Pack()),
    };
}//function

std::shared_ptr<Container> NeedlemanWunsch::getOutputType() const
{
    return std::shared_ptr<Container>(new Alignment());
}//function

int iDeletion = -16;//-50
int iInsertion = -16;//-50
int iDeletionContinued = -1;
int iInsertionContinued = -1;
int iMatch = 8;//20
int iMissMatch = -2;//-20

void needlemanWunsch(
        std::shared_ptr<NucSeq> pQuery, 
        std::shared_ptr<NucSeq> pRef,
        nucSeqIndex fromQuery,
        nucSeqIndex toQuery,
        nucSeqIndex fromRef,
        nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment
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
    std::vector<std::vector<int>> s(toQuery-fromQuery+1, std::vector<int>(toRef-fromRef+1));
    std::vector<std::vector<char>> dir(toQuery-fromQuery+1, std::vector<char>(toRef-fromRef+1));

    s[0][0] = 0;
    dir[0][0] = 1;
    s[1][0] = iInsertion;
    dir[1][0] = 2;
    s[0][1] = iDeletion;
    dir[0][1] = 3;
    for(nucSeqIndex uiI = 2; uiI < toQuery-fromQuery+1; uiI++)
    {
        s[uiI][0] = s[uiI - 1][0] + iInsertionContinued;
        dir[uiI][0] = 2;
    }//for
    for(nucSeqIndex uiI = 2; uiI < toRef-fromRef+1; uiI++)
    {
        s[0][uiI] = s[0][uiI - 1] + iDeletionContinued;
        dir[0][uiI] = 3;
    }//for
    for(nucSeqIndex uiI = 1; uiI < toQuery-fromQuery+1; uiI++)
    {
        for(nucSeqIndex uiJ = 1; uiJ < toRef-fromRef+1; uiJ++)
        {
            int newScore;
            //insertion
            if(dir[uiI - 1][uiJ] == 2)
                newScore = s[uiI - 1][uiJ] + iInsertionContinued;
            else
                newScore = s[uiI - 1][uiJ] + iInsertion;
            s[uiI][uiJ] = newScore;
            dir[uiI][uiJ] = 2;

            //deletion
            if(dir[uiI][uiJ - 1] == 3)
                newScore = s[uiI][uiJ - 1] + iDeletionContinued;
            else
                newScore = s[uiI][uiJ - 1] + iDeletion;
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
                newScore += iMissMatch;
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
}//function

std::shared_ptr<Container> NeedlemanWunsch::execute(
        ContainerVector vpInput
    )
{
    std::shared_ptr<Seeds> pSeeds = std::static_pointer_cast<Seeds>(vpInput[0]);
    std::shared_ptr<NucSeq> pQuery 
        = std::static_pointer_cast<NucSeq>(vpInput[1]);
    std::shared_ptr<Pack> pRefPack = 
        std::static_pointer_cast<Pack>(vpInput[2]);

    //no seeds => no spot found at all...
    if(pSeeds->empty())
    {
        std::shared_ptr<Alignment> pRet(new Alignment());
        return pRet;
    }//if

    //SOTING IS DONE IN THE LINESWEEP STEP
    //sort shadows (increasingly) by start coordinate of the match
    /*pSeeds->sort(
            [](Seed xA, Seed xB)
            {
                if(xA.start_ref() == xB.start_ref())
                    return xA.start() < xB.start();
                return xA.start_ref() < xB.start_ref();
            }//lambda
        );//sort function call*/

//DEPRECATED
#if 0
    //remove dangeling fronts and backs
    if(pSeeds->size() >= 2)
    {
        nucSeqIndex iCenter = 0;
        nucSeqIndex iSize = 0;
        nucSeqIndex iMaxDist = 10000;
        for(Seed& rSeed : *pSeeds)
        {
            if(rSeed.size() > iSize)
            {
                iSize = rSeed.size();
                iCenter = rSeed.start_ref() + rSeed.size()/2;
            }
        }
        while(pSeeds->size() > 1 && pSeeds->front().end_ref() + iMaxDist < iCenter)
        {
            DEBUG(
                std::cout << "WARNING: removed dangeling front" << std::endl;
            )
            pSeeds->pop_front();
        }//if
        while(pSeeds->size() > 1 && pSeeds->back().start_ref() > iCenter + iMaxDist)
        {
            DEBUG(
                std::cout << "WARNING: removed dangeling back" << std::endl;
            )
            pSeeds->pop_back();
        }//while
    }//if
    assert(pSeeds->size() > 0);
#endif

    DEBUG_2(
        std::cout << "seedlist: (start_ref, end_ref; start_query, end_query)" << std::endl;
        for(Seed& rSeed : *pSeeds)
        {
            std::cout << rSeed.start_ref() << ", " << rSeed.end_ref() << "; "
                << rSeed.start() << ", " << rSeed.end() << std::endl;
        }//for
    )

    nucSeqIndex beginQuery = pSeeds->front().start();
    nucSeqIndex beginRef = 0;
    if( pSeeds->front().start_ref() > beginQuery*2)
        beginRef = pSeeds->front().start_ref() - beginQuery*2;
    nucSeqIndex endQuery = pSeeds->back().end();
    assert(endQuery <= pQuery->length());
    nucSeqIndex endRef = pRefPack->uiUnpackedSizeForwardPlusReverse();
    if( 
            pSeeds->back().end_ref() + (pQuery->length()-endQuery)*2 <
            pRefPack->uiUnpackedSizeForwardPlusReverse()
        )
        endRef = pSeeds->back().end_ref() + (pQuery->length()-endQuery)*2;


    std::shared_ptr<Alignment> pRet(
            new Alignment(beginRef, endRef, iMatch, iMissMatch, iDeletion, iDeletionContinued)
        );

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

    for(Seed& rSeed : *pSeeds)
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

    needlemanWunsch(
        pQuery,
        pRef,
        endOfLastSeedQuery,
        pQuery->length(),
        endOfLastSeedReference,
        endRef - beginRef,
        pRet
    );

    //cleanup the alignment
    pRet->removeDangelingDeletions();

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
        "Picks a set of anchors for the strips of consideration.\n"
        "\n"
        "Execution:\n"
        "   Expects seg_list, query, ref\n"
        "       seg_list: the list of segments to pick the anchors from\n"
        "       query: the query as NucSeq\n"
        "       ref: the reference as Pack\n"
        "   returns alignment.\n"
        "       alignment: the final alignment\n"
    );
	boost::python::implicitly_convertible< 
		std::shared_ptr<NeedlemanWunsch>,
		std::shared_ptr<Module> 
	>();

}//function