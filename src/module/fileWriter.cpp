/** 
 * @file fileWriter.cpp
 * @author Markus Schmidt
 */
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


ContainerVector RadableFileWriter::getInputType() const
{
    return ContainerVector{
            std::shared_ptr<NucSeq>(new NucSeq()),
            std::shared_ptr<ContainerVector>(
                new ContainerVector(std::shared_ptr<Alignment>(new Alignment()))),
            std::shared_ptr<Pack>(new Pack())
        };
}//function

std::shared_ptr<Container> RadableFileWriter::getOutputType() const
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
            sCigar.append(std::to_string(std::get<1>(section)));
            switch (std::get<0>(section))
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

        if(pAlignment->bSecondary)
            flag |= SECONDARY_ALIGNMENT;

        std::string sRefName = pPack->nameOfSequenceForPosition(pAlignment->uiBeginOnRef);
        auto uiRefPos = pPack->posInSequence(pAlignment->uiBeginOnRef, pAlignment->uiEndOnRef);

        DEBUG(// check if the position that is saved to the file is correct
            bool bWrong = false;
            if(pPack->bPositionIsOnReversStrand(pAlignment->uiBeginOnRef))
            {
                //@todo frill in this self check...
            }//if
            else
            {
                if( pAlignment->uiBeginOnRef != pPack->startOfSequenceWithName(sRefName) + uiRefPos)
                    bWrong = true;
            }//else

            if(bWrong)
            {
                std::cerr << "Error: Tried to write wrong index to file" << std::endl;
                std::cout << "Have: " << sRefName
                          << " (= " << pPack->startOfSequenceWithName(sRefName) << ") "
                          << uiRefPos << std::endl;
                std::cout << "Wanted: " << pAlignment->uiBeginOnRef << " " << uiRefPos << std::endl;
                if(pPack->bPositionIsOnReversStrand(pAlignment->uiBeginOnRef))
                    std::cout << "Is reverse: True" << std::endl;
                else
                    std::cout << "Is reverse: False" << std::endl;
                exit(0);
            }//if
        )//DEBUG

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
            *pOut << std::to_string(uiRefPos) << "\t";
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


std::shared_ptr<Container> RadableFileWriter::execute(std::shared_ptr<ContainerVector> vpInput)
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
        
        std::string sPaired = "";
        //paired
        if(pAlignment->xStats.bPaired)
            sPaired = pAlignment->xStats.bFirst ? 
                "first mate of read pair" : "second mate of read pair";

        std::string sRefName = pPack->nameOfSequenceWithId(pPack->uiSequenceIdForPosition(pAlignment->uiBeginOnRef));
        //1 based index... 
        std::string sRefPos = std::to_string( 1 + pAlignment->uiBeginOnRef - pPack->startOfSequenceWithId(pPack->uiSequenceIdForPosition(pAlignment->uiBeginOnRef)));
        std::string sSegmentQuery = pQuery->fromTo(pAlignment->uiBeginOnQuery, pAlignment->uiEndOnQuery);
        std::string sQueryPos = std::to_string(pAlignment->uiBeginOnQuery);
        std::string sSegmentRef = pPack->vExtract(pAlignment->uiBeginOnRef, pAlignment->uiEndOnRef)->toString();
        std::string sQual = pQuery->fromToQual(pAlignment->uiBeginOnQuery, pAlignment->uiEndOnQuery);
        std::string sMapQual;
        if(std::isnan(pAlignment->fMappingQuality))
            sMapQual = "255";
        else
            sMapQual = std::to_string( (int)(pAlignment->fMappingQuality * 254));

        std::string sQueryLine;
        std::string sMatchLine;
        std::string sRefLine;
        unsigned int uiAlignmentCounter = 0;
        unsigned int uiQueryCounter = 0;
        unsigned int uiRefCounter = 0;

        {
            //synchronize output
            std::lock_guard<std::mutex> xGuard(*pLock);
            *pOut << "Score: " << std::to_string(pAlignment->score()) << "\nBegin on reference sequence: " << sRefName << " at position: " << sRefPos << "\nBegin on Query: " << sQueryPos << (pAlignment->bSecondary ? " Secondary\n" : "\n");
            for(std::tuple<MatchType, nucSeqIndex> section : pAlignment->data)
            {
                for(unsigned int i=0; i< std::get<1>(section);i++)
                {
                    if (uiAlignmentCounter % uiNucsPerLine == 0)
                        *pOut << std::to_string(uiAlignmentCounter)
                            << "-" << std::to_string(uiAlignmentCounter+uiNucsPerLine) << "\n";
                    switch(std::get<0>(section))
                    {
                        case MatchType::match:
                            sQueryLine += sSegmentQuery[uiQueryCounter++];
                            sRefLine += sSegmentRef[uiRefCounter++];
                            sMatchLine += "|";
                            break;
                        case MatchType::seed:
                            sQueryLine += sSegmentQuery[uiQueryCounter++];
                            sRefLine += sSegmentRef[uiRefCounter++];
                            sMatchLine += "I";
                            break;
                        case MatchType::missmatch:
                            sQueryLine += sSegmentQuery[uiQueryCounter++];
                            sRefLine += sSegmentRef[uiRefCounter++];
                            sMatchLine += " ";
                            break;
                        case MatchType::insertion:
                            sQueryLine += sSegmentQuery[uiQueryCounter++];
                            sRefLine += "-";
                            sMatchLine += " ";
                            break;
                        case MatchType::deletion:
                            sQueryLine += "-";
                            sRefLine += sSegmentRef[uiRefCounter++];
                            sMatchLine += " ";
                            break;
                    }//switch
                    if (++uiAlignmentCounter % uiNucsPerLine == 0)
                    {
                        *pOut << sQueryLine << "\tQuery\n";
                        *pOut << sMatchLine << "\n";
                        *pOut << sRefLine << "\tReference\n\n";
                        sQueryLine = "";
                        sRefLine = "";
                        sMatchLine = "";
                    }//if
                }//for
            }//for
            if(uiAlignmentCounter % uiNucsPerLine == 0)
                return std::shared_ptr<Container>(new Nil());
            while (uiAlignmentCounter++ % uiNucsPerLine != 0)
            {
                sQueryLine += "-";
                sRefLine += "-";
                sMatchLine += " ";
            }
            *pOut << sQueryLine << "\tQuery\n";
            *pOut << sMatchLine << "\n";
            *pOut << sRefLine << "\tReference\n\n";
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