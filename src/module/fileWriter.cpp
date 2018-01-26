#include "module/fileWriter.h"

using namespace libMA;

ContainerVector FileWriter::getInputType() const
{
    return ContainerVector{
            std::shared_ptr<NucSeq>(new NucSeq()),
            std::shared_ptr<ContainerVector>(
                new ContainerVector(std::shared_ptr<Alignment>(new Alignment()))),
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
    std::shared_ptr<ContainerVector> pAlignments =
        std::static_pointer_cast<ContainerVector>((*vpInput)[1]);
    std::shared_ptr<Pack> pPack =
        std::static_pointer_cast<Pack>((*vpInput)[2]);

    for(std::shared_ptr<Container> pA : *pAlignments)
    {
        std::shared_ptr<Alignment> pAlignment = std::static_pointer_cast<Alignment>(pA);
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

        std::string sRefName = pPack->nameOfSequenceWithId(pPack->uiSequenceIdForPosition(pAlignment->uiBeginOnRef));
        std::string sRefPos = std::to_string( pAlignment->uiBeginOnRef - pPack->startOfSequenceWithId(pPack->uiSequenceIdForPosition(pAlignment->uiBeginOnRef)));
        std::string sSegment = pQuery->fromTo(pAlignment->uiBeginOnQuery, pAlignment->uiEndOnQuery);
        std::string sQual = pQuery->fromToQual(pAlignment->uiBeginOnQuery, pAlignment->uiEndOnQuery);
        std::string sMapQual;
        if(std::isnan(pAlignment->fMappingQuality))
            sMapQual = "255";
        else
            sMapQual = std::to_string( (int)(pAlignment->fMappingQuality * 254));

        {
            //synchronize file output
            std::lock_guard<std::mutex> xGuard(*pLock);

            //print alignment
            //query name
            *pFile << pQuery->sName;
            *pFile << "\t";
            //alignment flag
            *pFile << "0";
            *pFile << "\t";
            //reference name
            *pFile << sRefName;
            *pFile << "\t";
            //pos
            *pFile << sRefPos;
            *pFile << "\t";
            //mapping quality
            *pFile << sMapQual;
            *pFile << "\t";
            //cigar
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
            *pFile << sSegment;
            *pFile << "\t";
            //ASCII of Phred-scaled base Quality+33
            *pFile << sQual;
            *pFile << std::endl;
        }//score xGuard
    }//for

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