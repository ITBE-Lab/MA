#include "needlemanWunsch.h"




std::vector<ContainerType> NeedlemanWunsch::getInputType()
{
    return std::vector<ContainerType>{
        //the sound strip of consideration
        ContainerType::stripOfConsideration,
        //the query sequence
        ContainerType:nucSeq,
        //the reference sequence
        ContainerType:packedNucSeq,
    };
}//function

std::vector<ContainerType> NeedlemanWunsch::getOutputType()
{
    return std::vector<ContainerType>{ContainerType::alignment};
}//function

int iDeletion = -15;
int iInsertion = -15;
int iDeletionContinued = -1;
int iInsertionContinued = -1;
int iMatch = 10;
int iMissMatch = -2;

void needlemanWunsch(
        std::shared_ptr<NucleotideSequence> pQuery, 
        std::shared_ptr<NucleotideSequence> pRef,
        nucSeqIndex fromQuery,
        nucSeqIndex toQuery,
        nucSeqIndex fromRef,
        nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment
    )
{
    if(toRef <= fromRef)
        if(toQuery <= fromQuery)
            return;
    DEBUG(
        std::cout << toQuery-fromQuery << std::endl;
        for(unsigned int i = fromQuery; i < toQuery; i++)
            std::cout << pQuery->charAt(i);
        std::cout << std::endl;
        std::cout << toRef-fromRef << std::endl;
        for(unsigned int i = fromRef; i < toRef; i++)
            std::cout << pRef->charAt(i);
        std::cout << std::endl;
    )//DEBUG
    if(toQuery <= fromQuery)
    {
        int iY = toRef-fromRef;
        while(iY > 0)
        {
            pAlignment->append(Alignment::MatchType::deletion);
            DEBUG(
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
            pAlignment->append(Alignment::MatchType::insertion);
            DEBUG(  
                std::cout << "I";
            )//DEBUG
            iX--;
        }//while
        return;
    }//if
    int s[toQuery-fromQuery+1][toRef-fromRef+1];
    char dir[toQuery-fromQuery+1][toRef-fromRef+1];//1=match; 2=ins; 3=del
    s[0][0] = 0;
    dir[0][0] = 1;
    s[1][0] = iInsertion;
    dir[1][0] = 2;
    s[0][1] = iDeletion;
    dir[0][1] = 3;
    for(unsigned int uiI = 2; uiI < toQuery-fromQuery+1; uiI++)
    {
        s[uiI][0] = s[uiI - 1][0] + iInsertionContinued;
        dir[uiI][0] = 2;
    }//for
    for(unsigned int uiI = 2; uiI < toRef-fromRef+1; uiI++)
    {
        s[0][uiI] = s[0][uiI - 1] + iDeletionContinued;
        dir[0][uiI] = 3;
    }//for
    for(unsigned int uiI = 1; uiI < toQuery-fromQuery+1; uiI++)
    {
        for(unsigned int uiJ = 1; uiJ < toRef-fromRef+1; uiJ++)
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

    DEBUG(
        for(unsigned int uiI = 0; uiI < toRef-fromRef+1; uiI++)
        {
            if(uiI == 0)
                std::cout << " \t \t";
            else
                std::cout << pRef->charAt(toRef - uiI) << "\t";
        }//for
        std::cout << std::endl;
        for(unsigned int uiI = 0; uiI < toQuery-fromQuery+1; uiI++)
        {
            if(uiI == 0)
                std::cout << " \t";
            else
                std::cout << pQuery->charAt(toQuery - uiI) << "\t";
            for(unsigned int uiJ = 0; uiJ < toRef-fromRef+1; uiJ++)
                std::cout << s[uiI][uiJ] << "\t";
            std::cout << std::endl;
        }//for
    )//DEBUG

    int iX = toQuery-fromQuery;
    int iY = toRef-fromRef;
    while(iX > 0 || iY > 0)
    {
        if(dir[iX][iY] == 1)
        {
            if( (*pQuery)[toQuery - iX] == (*pRef)[toRef - iY] )
            {
                pAlignment->append(Alignment::MatchType::match);
                DEBUG(
                    std::cout << "M";
                )//DEBUG
            }
            else
            {
                pAlignment->append(Alignment::MatchType::missmatch);
                DEBUG(
                    std::cout << "W";
                )//DEBUG
            }
            iX--;
            iY--;
        }//if
        else if(dir[iX][iY] == 3)
        {
            pAlignment->append(Alignment::MatchType::deletion);
            iY--;
            DEBUG(
                std::cout << "D";
            )//DEBUG        
        }//if
        else
        {
            pAlignment->append(Alignment::MatchType::insertion);
            iX--;
            DEBUG(
                std::cout << "I";
            )//DEBUG
        }//if
    }//while
    DEBUG(
        std::cout << std::endl;
    )//DEBUG
}//function

std::shared_ptr<Container> NeedlemanWunsch::execute(
        std::vector<std::shared_ptr<Container>> vpInput
    )
{

    std::shared_ptr<StripOfConsideration> pStrip = std::static_pointer_cast<StripOfConsideration>(vpInput[0]);
    std::shared_ptr<NucleotideSequence> pQuery 
        = std::static_pointer_cast<NucleotideSequence>(vpInput[1]);
    std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefPack = 
        std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[2]);

    unsigned int uiBegin = 0;
    unsigned int uiEnd = pStrip->numMatches()-1;
    while(!pStrip->matchEnabled(uiBegin))
        uiBegin++;
    while(!pStrip->matchEnabled(uiEnd))
        uiEnd--;

    nucSeqIndex beginQuery = pStrip->getMatch(uiBegin)->getPosOnQuery();
    nucSeqIndex beginRef = pStrip->getMatch(uiBegin)->getPosOnReference() - beginQuery*2;
    nucSeqIndex len = pStrip->getMatch(uiEnd)->getLength();
    nucSeqIndex endQuery = pStrip->getMatch(uiEnd)->getPosOnQuery() + len;
    nucSeqIndex endRef = pStrip->getMatch(uiEnd)->getPosOnReference() + len + (pQuery->length()-endQuery)*2;

    std::shared_ptr<NucleotideSequence> pRef = pRefPack->vExtract(beginRef, endRef);

    std::shared_ptr<Alignment> pRet(new Alignment(beginRef, endRef));

    nucSeqIndex endOfLastSeedQuery = 0;
    nucSeqIndex endOfLastSeedReference = 0;

    for(unsigned int uiI = uiBegin; uiI <= uiEnd; uiI++)
    {
        if(pStrip->matchEnabled(uiI))
        {
            needlemanWunsch(
                    pQuery,
                    pRef,
                    endOfLastSeedQuery,
                    pStrip->getMatch(uiI)->getPosOnQuery(),
                    endOfLastSeedReference,
                    pStrip->getMatch(uiI)->getPosOnReference() - beginRef,
                    pRet
                );
            unsigned int ovQ = endOfLastSeedQuery - pStrip->getMatch(uiI)->getPosOnQuery();
            if(pStrip->getMatch(uiI)->getPosOnQuery() > endOfLastSeedQuery)
                ovQ = 0;
            unsigned int ovR = endOfLastSeedReference - (pStrip->getMatch(uiI)->getPosOnReference() - beginRef);
            if(pStrip->getMatch(uiI)->getPosOnReference() - beginRef > endOfLastSeedReference)
                ovR = 0;
            unsigned int len = pStrip->getMatch(uiI)->getLength() - 1;
            nucSeqIndex overlap = std::max(ovQ, ovR);
            if(len > overlap)
            {
                pRet->append(Alignment::MatchType::match, len - overlap);
                DEBUG(
                    std::cout << len - overlap << std::endl;
                    for(unsigned int i = 0; i < len - overlap; i++)
                        std::cout << pQuery->charAt(i + pStrip->getMatch(uiI)->getPosOnQuery());
                    std::cout << std::endl;
                    for(unsigned int i = 0; i < len - overlap; i++)
                        std::cout << pRef->charAt(i + pStrip->getMatch(uiI)->getPosOnReference() - beginRef);
                    std::cout << std::endl;
                    for(unsigned int i = 0; i < len - overlap; i++)
                        std::cout << "m";
                )//DEBUG
            }//if
            if(ovQ > ovR)
                pRet->append(Alignment::MatchType::deletion, ovQ - ovR);
            DEBUG(
                for(unsigned int i = ovR; i < ovQ; i++)
                    std::cout << "d";
            )
            if(ovR > ovQ)
                pRet->append(Alignment::MatchType::insertion, ovR - ovQ);
            DEBUG(
                for(unsigned int i = ovQ; i < ovR; i++)
                    std::cout << "i";
                std::cout << std::endl;
            )//DEBUG
            endOfLastSeedQuery = pStrip->getMatch(uiI)->getPosOnQuery() + len;
            endOfLastSeedReference = pStrip->getMatch(uiI)->getPosOnReference() + len - beginRef;
        }//if
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

    std::cout << std::endl;
    return pRet;

}//function

void exportNeedlemanWunsch()
{
     //export the segmentation class
    boost::python::class_<
        NeedlemanWunsch, 
        boost::python::bases<Module>
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

}//function