#include "module/fileReader.h"

using namespace libMA;

ContainerVector FileReader::getInputType() const
{
    return ContainerVector{std::shared_ptr<Container>(new Nil())};
}//function

std::shared_ptr<Container> FileReader::getOutputType() const
{
    return std::shared_ptr<Container>(new NucSeq());
}//function

std::shared_ptr<Container> FileReader::execute(std::shared_ptr<ContainerVector> vpInput)
{
    std::shared_ptr<NucSeq> pRet(new NucSeq());
    //FASTA format
    if(pFile->good() && pFile->peek() == '>')
    {
        std::string sLine;
        std::getline (*pFile, sLine);
        pRet->sName = sLine.substr(1);
        while(pFile->good() && pFile->peek() != '>')
        {
            std::getline (*pFile, sLine);
            size_t uiLineSize = sLine.length();
            std::vector<uint8_t> xQuality(uiLineSize, 126);//uiLineSize uint8_t's with value 127
            pRet->vAppend((const uint8_t*)sLine.c_str(), xQuality.data(), uiLineSize);
        }//while
        pRet->vTranslateToNumericFormUsingTable(pRet->xNucleotideTranslationTable, 0);
        return pRet;
    }//if
    //FASTAQ format
    if(pFile->good() && pFile->peek() == '@')
    {
        std::string sLine;
        std::getline (*pFile, sLine);
        pRet->sName = sLine.substr(1);
        while(pFile->good() && pFile->peek() != '+')
        {
            std::getline (*pFile, sLine);
            size_t uiLineSize = sLine.length();
            std::vector<uint8_t> xQuality(uiLineSize, 126);//uiLineSize uint8_t's with value 127
            pRet->vAppend((const uint8_t*)sLine.c_str(), xQuality.data(), uiLineSize);
        }//while
        pRet->vTranslateToNumericFormUsingTable(pRet->xNucleotideTranslationTable, 0);
        //quality
        unsigned int uiPos = 0;
        while(pFile->good() && pFile->peek() != '@')
        {
            std::getline (*pFile, sLine);
            for(size_t i=0; i < sLine.length(); i++)
                pRet->quality(i + uiPos) = (uint8_t)sLine[i];
            uiPos += sLine.length();
        }//while
        return pRet;
    }//if
    //unknown format @todo print appropriate error message
    std::cout << "could not read next query" << std::endl;
    return std::shared_ptr<Nil>(new Nil());
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