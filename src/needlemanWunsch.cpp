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

int iDeletion = -1;
int iInsertion = -1;
int iMatch = 1;
int iMissMatch = -1;

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
    std::cout << toQuery-fromQuery << std::endl;
    for(unsigned int i = fromQuery; i < toQuery; i++)
        std::cout << pQuery->charAt(i);
    std::cout << std::endl;
    std::cout << toRef-fromRef << std::endl;
    for(unsigned int i = fromRef; i < toRef; i++)
        std::cout << pRef->charAt(i);
    std::cout << std::endl;
    if(toRef <= fromRef)
        if(toQuery <= fromQuery)
            return;
    if(toQuery <= fromQuery)
    {
        int iY = toRef-fromRef;
        while(iY > 0)
        {
            pAlignment->append(Alignment::MatchType::deletion);
            std::cout << "D";
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
            std::cout << "I";
            iX--;
        }//while
        return;
    }//if
    int s[toQuery-fromQuery+1][toRef-fromRef+1];
    for(unsigned int uiI = 0; uiI < toQuery-fromQuery+1; uiI++)
        s[uiI][0] = -uiI;
    for(unsigned int uiI = 0; uiI < toRef-fromRef+1; uiI++)
        s[0][uiI] = -uiI;
    for(unsigned int uiI = 1; uiI < toQuery-fromQuery+1; uiI++)
    {
        for(unsigned int uiJ = 1; uiJ < toRef-fromRef+1; uiJ++)
        {
            //insertion
            s[uiI][uiJ] = s[uiI - 1][uiJ - 1] + iInsertion;

            //deletion
            int newScore = s[uiI][uiJ - 1] + iDeletion;
            if(newScore > s[uiI][uiJ])
                s[uiI][uiJ] = newScore;
            //match / missmatch
            newScore = s[uiI - 1][uiJ - 1];
            if( (*pQuery)[toQuery - uiI] == (*pRef)[toRef - uiJ] )
                newScore += iMatch;
            else
                newScore += iMissMatch;
            if(newScore > s[uiI][uiJ])
                s[uiI][uiJ] = newScore;
        }//for
    }//for

    int iX = toQuery-fromQuery;
    int iY = toRef-fromRef;
    while(iX > 0 && iY > 0)
    {
        //std::cout << pQuery->charAt(toQuery - iX) << ":" << pRef->charAt(toRef - iY) << std::endl;
        if(s[iX-1][iY-1] >= s[iX-1][iY] && s[iX-1][iY-1] >= s[iX][iY-1])
        {
            if( (*pQuery)[toQuery - iX] == (*pRef)[toRef - iY] )
            {
                pAlignment->append(Alignment::MatchType::match);
                std::cout << "M";
            }
            else
            {
                pAlignment->append(Alignment::MatchType::missmatch);
                std::cout << "W";
            }
            iX--;
            iY--;
        }//if
        else if(s[iX][iY-1] >= s[iX-1][iY])
        {
            pAlignment->append(Alignment::MatchType::deletion);
            iY--;
            std::cout << "D";
        }//if
        else
        {
            pAlignment->append(Alignment::MatchType::insertion);
            iX--;
            std::cout << "I";
        }//if
        //std::cout << std::endl;
    }//while
    while(iX > 0)
    {
        pAlignment->append(Alignment::MatchType::insertion);
        iX--;
        std::cout << "I";
    }//while
    while(iY > 0)
    {
        pAlignment->append(Alignment::MatchType::deletion);
        iY--;
        std::cout << "D";
    }//while
    std::cout << std::endl;
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

    nucSeqIndex beginRef = pStrip->getMatch(uiBegin)->getPosOnReference() - 10;
    nucSeqIndex len = pStrip->getMatch(uiEnd)->getLength();
    nucSeqIndex endQuery = pStrip->getMatch(uiEnd)->getPosOnQuery() + len;
    nucSeqIndex endRef = pStrip->getMatch(uiEnd)->getPosOnReference() + len + 10;

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
            pRet->append(Alignment::MatchType::deletion, ovQ);
            //for(unsigned int i = 0; i < ovQ; i++)
            //    std::cout << "d";
            pRet->append(Alignment::MatchType::insertion, ovR);
            //for(unsigned int i = 0; i < ovR; i++)
            //    std::cout << "i";
            if(len > ovQ + ovR)
            {
                pRet->append(Alignment::MatchType::match, len - ovQ - ovR);
                std::cout << len - ovQ - ovR << std::endl;
                for(unsigned int i = 0; i < len - ovQ - ovR; i++)
                    std::cout << pQuery->charAt(i + pStrip->getMatch(uiI)->getPosOnQuery());
                std::cout << std::endl;
                for(unsigned int i = 0; i < len - ovQ - ovR; i++)
                    std::cout << pRef->charAt(i + pStrip->getMatch(uiI)->getPosOnReference() - beginRef);
                std::cout << std::endl;
                for(unsigned int i = 0; i < len - ovQ - ovR; i++)
                    std::cout << "m";
                std::cout << std::endl;
            }//if
            std::cout << std::endl;
            endOfLastSeedQuery = pStrip->getMatch(uiI)->getPosOnQuery() + len;
            endOfLastSeedReference = pStrip->getMatch(uiI)->getPosOnReference() + len - beginRef;
        }//if
    }//for

    needlemanWunsch(
        pQuery,
        pRef,
        endOfLastSeedQuery,
        endQuery,
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