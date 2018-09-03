/** 
 * @file fileWriter.cpp
 * @author Markus Schmidt
 */
#include "module/dbWriter.h"

#ifdef WITH_POSTGRES
using namespace libMA;

ContainerVector DbWriter::getInputType() const
{
    return ContainerVector{
            std::shared_ptr<NucSeq>(new NucSeq()),
            std::shared_ptr<ContainerVector>(
                new ContainerVector(std::shared_ptr<Alignment>(new Alignment()))),
            std::shared_ptr<Pack>(new Pack())
        };
}//function

std::shared_ptr<Container> DbWriter::getOutputType() const
{
    return std::shared_ptr<Container>(new Nil());
}//function

std::shared_ptr<Container> DbWriter::execute(std::shared_ptr<ContainerVector> vpInput)
{
    std::shared_ptr<NucSeq> pQuery =
        std::dynamic_pointer_cast<NucSeq>((*vpInput)[0]); // dc
    std::shared_ptr<ContainerVector> pAlignments =
        std::dynamic_pointer_cast<ContainerVector>((*vpInput)[1]); // dc
    std::shared_ptr<Pack> pPack =
        std::dynamic_pointer_cast<Pack>((*vpInput)[2]); // dc

    for(std::shared_ptr<Container> pA : *pAlignments)
    {
        std::shared_ptr<Alignment> pAlignment = std::dynamic_pointer_cast<Alignment>(pA); // dc
        if(pAlignment->length() == 0)
            continue;
        std::string sCigar = pAlignment->cigarString(*pPack);

        uint32_t flag = pAlignment->getSamFlag(*pPack);
        
        std::string sContigOther = "*";
        std::string sPosOther = "0";
        // paired
        if(!pAlignment->xStats.pOther.expired())
        {
            flag |= pAlignment->xStats.bFirst ? FIRST_IN_TEMPLATE : LAST_IN_TEMPLATE;
            flag |= MULTIPLE_SEGMENTS_IN_TEMPLATE | SEGMENT_PROPERLY_ALIGNED;

            sContigOther = pAlignment->xStats.pOther.lock()->getContig(*pPack);
            sPosOther = std::to_string(pAlignment->xStats.pOther.lock()->getSamPosition(*pPack));
        }// if

        std::string sRefName = pAlignment->getContig(*pPack);
        // sam file format has 1-based indices bam 0-based...
        auto uiRefPos = pAlignment->getSamPosition(*pPack);

        double fMappingQual = pAlignment->fMappingQuality;

        std::string sSegment = pAlignment->getQuerySequence(*pQuery, *pPack);

        std::string sSQL = "";

        sSQL += "INSERT INTO alignment (cigar, position, mapping_quality, query_id, run_id, sam_flags, contig, query_sequence, contig_other, position_other, length) VALUES ( \'";

        sSQL += sCigar + "\', ";
        sSQL += std::to_string(uiRefPos) + ", ";
        sSQL += std::to_string(fMappingQual) + ", \'";
        sSQL += pAlignment->xStats.sName + "\', ";
        sSQL += std::to_string(iRunId) + ", ";
        sSQL += std::to_string(flag) + ", \'";
        sSQL += sRefName + "\', \'";
        sSQL += sSegment + "\', ";
        sSQL += sContigOther + "\', ";
        sSQL += sPosOther + "\', ";
        sSQL += std::to_string(pAlignment->uiEndOnQuery - pAlignment->uiBeginOnQuery) + " )";

        //std::cerr << sSQL << std::endl;

        xConnection.exec(sSQL);
    }//for

    return std::shared_ptr<Container>(new Nil());
}//function

#ifdef WITH_PYTHON
void exportDBWriter()
{
    //export the FileWriter class
    boost::python::class_<
            DbWriter, 
            boost::python::bases<Module>, 
            std::shared_ptr<DbWriter>
        >("DbWriter", boost::python::init<std::string, uint32_t>())
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<DbWriter>,
        std::shared_ptr<Module> 
    >();

}//function

#endif // WITH_PYTHON

#endif // WITH_POSTGRES