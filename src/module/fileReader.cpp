#include "module/fileReader.h"

using namespace libMABS;

ContainerVector FileReader::getInputType() const
{
    return ContainerVector{std::shared_ptr<Container>(new Nil())};
}//function

std::shared_ptr<Container> FileReader::getOutputType() const
{
    return std::shared_ptr<Container>(new ContainerVector(std::shared_ptr<NucSeq>(new NucSeq())));
}//function


void FileReader::getLines(std::function<void(std::string&)> fDo)
{
    std::ifstream xFile (sFileName);
    if (xFile.is_open())
    {
        std::string sLine;
        while ( std::getline (xFile, sLine) )
            fDo(sLine);
        xFile.close();
    }//if
    else std::cout << "Unable to open file"; 
}//function

std::shared_ptr<Container> FileReader::execute(std::shared_ptr<ContainerVector> vpInput)
{
    std::shared_ptr<ContainerVector> sRet(
            new ContainerVector(std::shared_ptr<NucSeq>(new NucSeq()))
        );
    bool bFASTA = false;
    bool bFASTAQ = false;
    bool bQualNext = false;
    bool bKnowFormat = false;
    std::shared_ptr<NucSeq> pLast;
    getLines(
        [&]
        (std::string& sLine)
        {
            if(sLine.length() == 0)
                return;

            //figure out the format if that has not already happened
            if(!bKnowFormat)
            {
                if(sLine[0] == '>')
                    bFASTA = true;
                if(sLine[0] == '@')
                    bFASTAQ = true;
                bKnowFormat = true;
            }//if

            //actually parse the file
            if( (bFASTAQ && sLine[0] == '@') || (bFASTA && sLine[0] == '>') )
            {
                if(pLast != nullptr)
                {
                    pLast->vTranslateToNumericFormUsingTable(pLast->xNucleotideTranslationTable, 0);
                    sRet->push_back(pLast);
                }//if
                pLast = std::shared_ptr<NucSeq>(new NucSeq());
                pLast->sName = sLine.substr(1);
                bQualNext = false;
            }//if
            else if(bFASTAQ && sLine[0] == '+')
                bQualNext = true;
            else if(bFASTAQ && bQualNext)
            {
                for(size_t i=0; i < sLine.length(); i++)
                    pLast->quality(i) = (uint8_t)sLine[i];
            }//else
            else
            {
                size_t uiLineSize = sLine.length();
                std::vector<uint8_t> xQuality(uiLineSize, 126);//uiLineSize uint8_t's with value 127
                pLast->vAppend((const uint8_t*)sLine.c_str(), xQuality.data(), uiLineSize);
            }//else
        }//lambda
    );
    if(pLast != nullptr)
    {
        pLast->vTranslateToNumericFormUsingTable(pLast->xNucleotideTranslationTable, 0);
        sRet->push_back(pLast);
    }//if
    return sRet;
}//function

void exportFileReader()
{
    //export the FileReader class
    boost::python::class_<
            FileReader, 
            boost::python::bases<Module>, 
            std::shared_ptr<FileReader>
        >("FileReader", boost::python::init<std::string>())
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<FileReader>,
        std::shared_ptr<Module> 
    >();

}//function