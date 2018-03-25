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


NeedlemanWunsch::NeedlemanWunsch(bool bLocal)
        :
    bLocal(bLocal)
{
    matrix = parasail_matrix_create("ACGT", iMatch, -iMissMatch);
}//constructor

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
 * 
 * @TODO: at the moment ugly C code is used here (free functions).. find a way to replace that?
 */
nucSeqIndex NeedlemanWunsch::needlemanWunsch(
        std::shared_ptr<NucSeq> pQuery, 
        std::shared_ptr<NucSeq> pRef,
        nucSeqIndex fromQuery,
        nucSeqIndex toQuery,
        nucSeqIndex fromRef,
        nucSeqIndex toRef,
        std::shared_ptr<Alignment> pAlignment,
        bool bNoGapAtBeginning, //disabled at the moment
        bool bNoGapAtEnd //disabled at the moment
        DEBUG_PARAM(bool bPrintMatrix)
    )
{
    const char* seqA = pQuery->fromTo(fromQuery, toQuery).c_str();
    const char* seqB = pRef->fromTo(fromRef, toRef).c_str();
    //do the alignment
    parasail_result_t* pResult = parasail_nw_trace_scan_32(
            seqA, toQuery - fromQuery,
            seqB, toRef - fromRef,
            iGap, iExtend, matrix
        );

    //get the cigar
    parasail_cigar_t* pCigar = parasail_result_get_cigar(pResult, seqA, toQuery - fromQuery, 
        seqB, toRef - fromRef, matrix);

    //decode the cigar
    for(int i = 0; i < pCigar->len; i++)
    {
        char c = parasail_cigar_decode_op(pCigar->seq[i]);
        uint32_t uiLen = parasail_cigar_decode_len(pCigar->seq[i]);
        switch (c)
        {
            case 'M':
                pAlignment->append(MatchType::match, uiLen);
                break;
            case 'I':
                pAlignment->append(MatchType::insertion, uiLen);
                break;
            case 'D':
                pAlignment->append(MatchType::deletion, uiLen);
                break;
            case 'X':
                pAlignment->append(MatchType::missmatch, uiLen);
                break;
        }//switch
    }//for

    parasail_cigar_free(pCigar);
    parasail_result_free(pResult);
    return 0;
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