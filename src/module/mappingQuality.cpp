/** 
 * @file mappinQuality.cpp
 * @author Markus Schmidt
 */
#include "module/mappingQuality.h"

using namespace libMA;

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

    //if no alignment was found we cannot set any quality...
    if(pAlignments->size() == 0)
        return std::shared_ptr<ContainerVector>(
            new ContainerVector(std::shared_ptr<Alignment>(new Alignment())));

    //set the mapping quality for all alignments to zero
    //mapping qual of best one will be overridden
    for(auto pAlign : *pAlignments)
    {
        std::shared_ptr<Alignment> pCasted = std::static_pointer_cast<Alignment>(pAlign);
        pCasted->fMappingQuality = 0.0;
        pCasted->bSecondary = true;
    }//for

    //compute the mapping quality for the best alignment
    std::shared_ptr<Alignment> pFirst = std::static_pointer_cast<Alignment>((*pAlignments)[0]);
    pFirst->bSecondary = false;

    // mapping quality based on scores
    if(pAlignments->size() >= 2)
    {
        std::shared_ptr<Alignment> pSecond = std::static_pointer_cast<Alignment>((*pAlignments)[1]);
        // this formula is given in the paper and is very similar to Heng li's approach in BWA-SW
        pFirst->fMappingQuality =
                ( pFirst->score() - pSecond->score() )
                    /
                (double)( pFirst->score() )
            ;
    }//if
    else
        // the score of the second best alignment is 0 if we do not even find one...
        pFirst->fMappingQuality = 1;// pFirst->score() / (double)(iMatch * pQuery->length());

    // factors
    // penalty for too little seeds
    // (this improves the mapping quality estimation quite significantly)
    //double dA = std::max(std::min(5* pFirst->numBySeeds() / (double)pQuery->length(), 1.0), 0.01);
    //pFirst->fMappingQuality *= dA;

    /// maybe this should be moved into it's own module but whatever...
    auto pRet = std::shared_ptr<ContainerVector>(new ContainerVector(pAlignments));
    if(uiReportNBest != 0 && pRet->size() > uiReportNBest)
    {
        //remove the smallest elements
        pRet->erase(pRet->begin()+uiReportNBest, pRet->end());
        assert(pRet->size() == uiReportNBest);
    }//if
    return pRet;
}//function

void exportMappingQuality()
{
    //export the MappingQuality class
    boost::python::class_<
            MappingQuality, 
            boost::python::bases<Module>, 
            std::shared_ptr<MappingQuality>
        >("MappingQuality", boost::python::init<unsigned int>())
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<MappingQuality>,
        std::shared_ptr<Module> 
    >();
}//function
