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
        if(pAlignment->length() == 0)
            continue;
        std::string sCigar = "";
        for(std::tuple<MatchType, nucSeqIndex> section : pAlignment->data)
        {
            sCigar.append(std::to_string(std::get<0>(section)));
            switch (std::get<1>(section))
            {
                case MatchType::seed:
                case MatchType::match:
                    sCigar.append("=");
                    break;
                case MatchType::missmatch:
                    sCigar.append("X");
                    break;
                case MatchType::insertion:
                    sCigar.append("I");
                    break;
                case MatchType::deletion:
                    sCigar.append("D");
                    break;
            }//switch
        }//for

        char flag = 0;

        if(pPack->bPositionIsOnReversStrand(pAlignment->uiBeginOnRef))
            flag |= REVERSE_COMPLEMENTED;

        //paired
        if(pAlignment->xStats.bPaired)
        {
            flag |= pAlignment->xStats.bFirst ? FIRST_IN_TEMPLATE : LAST_IN_TEMPLATE;
            flag |= MULTIPLE_SEGMENTS_IN_TEMPLATE | SEGMENT_PROPERLY_ALIGNED;
        }//if

        std::string sRefName = pPack->nameOfSequenceWithId(pPack->uiSequenceIdForPosition(pAlignment->uiBeginOnRef));
        //1 based index... 
        std::string sRefPos = std::to_string( 1 + pAlignment->uiBeginOnRef - pPack->startOfSequenceWithId(pPack->uiSequenceIdForPosition(pAlignment->uiBeginOnRef)));
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
            *pOut << pQuery->sName << "\t";
            //alignment flag
            *pOut << std::to_string(flag) << "\t";
            //reference name
            *pOut << sRefName << "\t";
            //pos
            *pOut << sRefPos << "\t";
            //mapping quality
            *pOut << sMapQual << "\t";
            //cigar
            *pOut << sCigar  << "\t";
            //Ref. name of the mate/next read ? wut? @todo
            *pOut << "*" << "\t";
            //Position of the mate/next read ? wut? @todo
            *pOut << "0" << "\t";
            //observed Template length
            *pOut << std::to_string(pAlignment->length()) << "\t";
            //segment sequence
            *pOut << sSegment << "\t";
            //ASCII of Phred-scaled base Quality+33
            *pOut << sQual << "\n";
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