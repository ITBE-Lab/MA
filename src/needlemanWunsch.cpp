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
    if(toQuery <= fromQuery-2)
    {
        int iY = toRef-fromRef-2;
        while(iY > 0)
        {
            pAlignment->append(Alignment::MatchType::deletion);
            iY--;
        }//while
        return;
    }//if
    if(toRef <= fromRef-2)
    {
        int iX = toQuery-fromQuery-2;
        while(iX > 0)
        {
            pAlignment->append(Alignment::MatchType::insertion);
            iX--;
        }//while
        return;
    }//if
    int s[toQuery-fromQuery-1][toRef-fromRef-1];
    for(unsigned int uiI = 0; uiI < toQuery-fromQuery-1; uiI++)
        s[uiI][0] = -uiI;
    for(unsigned int uiI = 0; uiI < toRef-fromRef-1; uiI++)
        s[0][uiI] = -uiI;
    for(unsigned int uiI = 1; uiI < toQuery-fromQuery-1; uiI++)
    {
        for(unsigned int uiJ = 1; uiJ < toRef-fromRef-1; uiJ++)
        {
            //insertion
            s[uiI][uiJ] = s[uiI - 1][uiJ - 1] + iInsertion;

            //deletion
            int newScore = s[uiI][uiJ - 1] + iDeletion;
            if(newScore > s[uiI][uiJ])
                s[uiI][uiJ] = newScore;

            //match / missmatch
            newScore = s[uiI - 1][uiJ - 1];
            if( (*pQuery)[uiI] == (*pRef)[uiJ] )
                newScore += iMatch;
            else
                newScore += iMissMatch;
            if(newScore > s[uiI][uiJ])
                s[uiI][uiJ] = newScore;

        }//for
    }//for

    int iX = toQuery-fromQuery-2;
    int iY = toRef-fromRef-2;
    while(iX != 0 || iY != 0)
    {
        if(s[iX-1][iY-1] <= s[iX-1][iY] && s[iX-1][iY-1] <= s[iX][iY-1])
        {
            if( (*pQuery)[iX] == (*pRef)[iY] )
                pAlignment->append(Alignment::MatchType::match);
            else
                pAlignment->append(Alignment::MatchType::missmatch);
            iX--;
            iY--;
        }//if
        else if(s[iX][iY-1] <= s[iX-1][iY])
        {
            pAlignment->append(Alignment::MatchType::deletion);
            iY--;
        }//if
        else
        {
            pAlignment->append(Alignment::MatchType::insertion);
            iX--;
        }//if
    }//while
    while(iX != 0)
    {
        pAlignment->append(Alignment::MatchType::insertion);
        iX--;
    }//while
    while(iY != 0)
    {
        pAlignment->append(Alignment::MatchType::deletion);
        iY--;
    }//while
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
    nucSeqIndex beginRef = pStrip->getMatch(uiBegin)->getPosOnReference();
    nucSeqIndex len = pStrip->getMatch(uiEnd)->getLength();
    nucSeqIndex endQuery = pStrip->getMatch(uiEnd)->getPosOnQuery() + len;
    nucSeqIndex endRef = pStrip->getMatch(uiEnd)->getPosOnReference() + len;

    std::shared_ptr<NucleotideSequence> pRef = pRefPack->vExtract(beginRef, endRef);

    nucSeqIndex totalSize = std::max(endQuery - beginQuery, endRef - beginRef);

    std::shared_ptr<Alignment> pRet(new Alignment(totalSize, beginRef));

    nucSeqIndex endOfLastSeedQuery = 0;
    nucSeqIndex endOfLastSeedReference = 0;

    for(unsigned int uiI = uiBegin; uiI < uiEnd; uiI++)
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
            unsigned int len = pStrip->getMatch(uiI)->getLength();
            pRet->append(Alignment::MatchType::match, len);
            endOfLastSeedQuery = pStrip->getMatch(uiI)->getPosOnQuery() + len;
            endOfLastSeedReference = pStrip->getMatch(uiI)->getPosOnReference() + len;
        }//if
    }//for

    //TODO: missig last needleman wunsch

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