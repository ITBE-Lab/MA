#include "module/fileWriter.h"

using namespace libMABS;

ContainerVector FileWriter::getInputType() const
{
    return ContainerVector{
            std::shared_ptr<NucSeq>(new NucSeq()),
            std::shared_ptr<Alignment>(new Alignment()),
            std::shared_ptr<Pack>(new Pack())
        };
}//function

std::shared_ptr<Container> FileWriter::getOutputType() const
{
    return std::shared_ptr<Container>(new Nil());
}//function

std::shared_ptr<Container> FileWriter::execute(std::shared_ptr<ContainerVector> vpInput)
{
    std::shared_ptr<NucSeq> pQuery =
        std::static_pointer_cast<NucSeq>((*vpInput)[0]);
    std::shared_ptr<Alignment> pAlignment =
        std::static_pointer_cast<Alignment>((*vpInput)[1]);
    std::shared_ptr<Pack> pPack =
        std::static_pointer_cast<Pack>((*vpInput)[2]);

    //print alignment
    //query name
    *pFile << pQuery->sName;
    *pFile << "\t";
    //alignment flag
    *pFile << "0";
    *pFile << "\t";
    //reference name
    *pFile << pPack->nameOfSequenceWithId(pPack->uiSequenceIdForPosition(pAlignment->uiBeginOnRef));
    *pFile << "\t";
    //pos
    *pFile << pAlignment->uiBeginOnRef - pPack->startOfSequenceWithId(pPack->uiSequenceIdForPosition(pAlignment->uiBeginOnRef));
    *pFile << "\t";
    //mapping quality
    if(std::isnan(pAlignment->fMappingQuality))
        *pFile << "255";
    else
        *pFile << std::to_string( (int)(pAlignment->fMappingQuality * 254));
    *pFile << "\t";
    //cigar
    std::string sCigar = "";
    unsigned int iCont = 0;
    for(std::tuple<MatchType, nucSeqIndex> section : pAlignment->data)
    {
        if(iCont++ < pAlignment->uiDataStart)
            continue;
        sCigar += std::to_string(std::get<0>(section));
        switch (std::get<1>(section))
        {
            case MatchType::seed:
            case MatchType::match:
                sCigar += "=";
                break;
            case MatchType::missmatch:
                sCigar += "X";
                break;
            case MatchType::insertion:
                sCigar += "I";
                break;
            case MatchType::deletion:
                sCigar += "D";
                break;
        }//switch
    }//for
    *pFile << sCigar;
    *pFile << "\t";
    //Ref. name of the mate/next read ? wut? @todo
    *pFile << "*";
    *pFile << "\t";
    //Position of the mate/next read ? wut? @todo
    *pFile << "0";
    *pFile << "\t";
    //observed Template length
    *pFile << std::to_string(pAlignment->length());
    *pFile << "\t";
    //segment sequence
    *pFile << pQuery->fromTo(pAlignment->uiBeginOnQuery, pAlignment->uiEndOnQuery);
    *pFile << "\t";
    //ASCII of Phred-scaled base QUALity+33 ? wut? @todo
    *pFile << pQuery->fromToQual(pAlignment->uiBeginOnQuery, pAlignment->uiEndOnQuery);
    *pFile << std::endl;

    return std::shared_ptr<Container>(new Nil());
}//function

void exportFileWriter()
{
    //export the FileWriter class
    boost::python::class_<
            FileWriter, 
            boost::python::bases<Module>, 
            std::shared_ptr<FileWriter>
        >("FileWriter", boost::python::init<std::string>())
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<FileWriter>,
        std::shared_ptr<Module> 
    >();

}//function