/** 
 * @file mappinQuality.cpp
 * @author Markus Schmidt
 */
#include "module/mappingQuality.h"

using namespace libMA;

using namespace libMA::defaults;
extern int iMatch;

ContainerVector MappingQuality::getInputType() const
{
    return ContainerVector{
        std::shared_ptr<Container>(new NucSeq()),
        std::shared_ptr<Container>(new ContainerVector(std::shared_ptr<Alignment>(new Alignment())))
    };
}//function

std::shared_ptr<Container> MappingQuality::getOutputType() const
{
    return std::shared_ptr<Container>(new ContainerVector(std::shared_ptr<Alignment>(new Alignment())));
}//function

std::shared_ptr<Container> MappingQuality::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<NucSeq> pQuery = std::static_pointer_cast<NucSeq>((*vpInput)[0]);
    std::shared_ptr<ContainerVector> pAlignments = std::static_pointer_cast<ContainerVector>((*vpInput)[1]);
    auto pSupplementaries = std::make_shared<ContainerVector>(std::make_shared<Alignment>());

    //if no alignment was found we cannot set any quality...
    if(pAlignments->size() == 0)
        return std::shared_ptr<ContainerVector>(
            new ContainerVector(std::shared_ptr<Alignment>(new Alignment())));

    //compute the mapping quality for the best alignment
    std::shared_ptr<Alignment> pFirst = std::static_pointer_cast<Alignment>((*pAlignments)[0]);
    pFirst->bSecondary = false;

    //set the mapping quality for all alignments to zero
    //mapping qual of best one will be overridden
    size_t uiSupplementaries = 0;
    bool bFirst = true;
    for(auto pAlign : *pAlignments)
    {
        if(bFirst)
        {
            bFirst = false;
            continue;
        }// if
        std::shared_ptr<Alignment> pCasted = std::static_pointer_cast<Alignment>(pAlign);
        pCasted->fMappingQuality = 0.0;
        // @todo this needs to be redone properly
        // at the moment we allow merely one supplementary alignment
        if (uiSupplementaries == 0 && pCasted->noOverlap(*pFirst))
        {
            pCasted->bSupplementary = true;
            pCasted->bSecondary = false;
            uiSupplementaries++;
        }// if
        else
        {
            pCasted->bSupplementary = false;
            pCasted->bSecondary = true;
        }// else
    }//for


    // mapping quality based on scores
    if(pAlignments->size() - uiSupplementaries >= 2)
    {
        size_t uiI = 1;
        std::shared_ptr<Alignment> pSecond;
        while(pSecond == nullptr || pSecond->bSupplementary)
        {
            pSecond = std::static_pointer_cast<Alignment>(
                    (*pAlignments)[uiI]
                );
            uiI++;
        }// while
        // this formula is given in the paper and is very similar to Heng li's approach in BWA-SW
        pFirst->fMappingQuality =
                static_cast<double>( pFirst->score() - pSecond->score() )
                    /
                static_cast<double>( pFirst->score() )
            ;
    }//if
    else
        // the score of the second best alignment is 0 if we do not even find one...
        pFirst->fMappingQuality = 1;// pFirst->score() / (double)(iMatch * pQuery->length());

    if(uiSupplementaries > 0)
    {
        // set supp mapping quality
        for(size_t uiI = 1; uiI < pAlignments->size(); uiI++)
        {
            std::shared_ptr<Alignment> pCasted = std::static_pointer_cast<Alignment>(
                    (*pAlignments)[uiI]
                );
            if(pCasted->bSupplementary)
                pCasted->fMappingQuality = pFirst->fMappingQuality;
        }// for
        // move supplementary alignments forward
        std::sort(
            pAlignments->begin(), pAlignments->end(),
            []
            (std::shared_ptr<Container> a, std::shared_ptr<Container> b)
            {
                return a->larger(b);
            }//lambda
        );//sort function call
    }// if

    // factors
    // penalty for too little seeds
    // (this improves the mapping quality estimation quite significantly)
    //double dA = std::max(std::min(5* pFirst->numBySeeds() / (double)pQuery->length(), 1.0), 0.01);
    //pFirst->fMappingQuality *= dA;

    /// maybe this should be moved into it's own module but whatever...
    auto pRet = std::shared_ptr<ContainerVector>(new ContainerVector(pAlignments));
    if(uiReportNBest != 0 && pRet->size() > uiReportNBest+uiSupplementaries)
    {
        //remove the smallest elements
        pRet->erase(pRet->begin()+uiReportNBest+uiSupplementaries, pRet->end());
        assert(pRet->size() == uiReportNBest+uiSupplementaries);
    }//if

    // remove secondary with too small scores
    while(pRet->size() > 1)
    {
        std::shared_ptr<Alignment> pCasted = std::static_pointer_cast<Alignment>(pRet->back());
        if(!pCasted->bSecondary)
            break;
        if(pCasted->fMappingQuality > pFirst->fMappingQuality * fMinSecScoreRatio)
            break;
        pRet->pop_back();
    }// while

    return pRet;
}//function

#ifdef WITH_PYTHON
void exportMappingQuality()
{
    //export the MappingQuality class
    boost::python::class_<
            MappingQuality, 
            boost::python::bases<Module>, 
            std::shared_ptr<MappingQuality>
        >("MappingQuality")
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<MappingQuality>,
        std::shared_ptr<Module> 
    >();
}//function
#endif
