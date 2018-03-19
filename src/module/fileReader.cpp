/** 
 * @file fileReader.cpp
 * @author Markus Schmidt
 */
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

size_t len(std::string& sLine)
{
    size_t uiLineSize = sLine.length();
    while( 
            uiLineSize > 0 &&
            sLine[uiLineSize-1] != 'A' &&
            sLine[uiLineSize-1] != 'C' &&
            sLine[uiLineSize-1] != 'T' &&
            sLine[uiLineSize-1] != 'G' &&
            sLine[uiLineSize-1] != 'a' &&
            sLine[uiLineSize-1] != 'c' &&
            sLine[uiLineSize-1] != 't' &&
            sLine[uiLineSize-1] != 'g'
        )
        uiLineSize--;
    return uiLineSize;
}

std::shared_ptr<Container> FileReader::execute(std::shared_ptr<ContainerVector> vpInput)
{
    std::shared_ptr<NucSeq> pRet(new NucSeq());
    //FASTA format
    if(pFile->good() && !pFile->eof() && pFile->peek() == '>')
    {
        std::string sLine;
        std::getline (*pFile, sLine);
        //make sure that the name contains no spaces
        //in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr(1, sLine.find(' '));
        while(pFile->good() && !pFile->eof() && pFile->peek() != '>')
        {
            std::getline (*pFile, sLine);
            DEBUG(
                for(auto character : sLine)
                {
                    bool bOkay = false;
                    for(char c : std::vector<char>{ 'A', 'C', 'T', 'G', 'a', 'c', 't', 'g'})
                        if(c == character)
                            bOkay = true;
                    if(!bOkay)
                    {
                        std::cerr << "Invalid symbol in fasta: " << sLine << std::endl;
                        throw AlignerException("Invalid symbol in fasta");
                    }//if
                }//for
            )//DEBUG
            size_t uiLineSize = len(sLine);
            std::vector<uint8_t> xQuality(uiLineSize, 126);//uiLineSize uint8_t's with value 127
            pRet->vAppend((const uint8_t*)sLine.c_str(), xQuality.data(), uiLineSize);
        }//while
        pRet->vTranslateToNumericFormUsingTable(pRet->xNucleotideTranslationTable, 0);

        //run self tests for the nucSeq
        DEBUG_2(
            std::cout << pRet->fastaq() << std::endl;
        )//DEBUG_@
        DEBUG(
            pRet->check();
        )//DEBUG

        return pRet;
    }//if
    //FASTAQ format
    if(pFile->good() && !pFile->eof() && pFile->peek() == '@')
    {
        std::string sLine;
        std::getline (*pFile, sLine);
        //make sure that the name contains no spaces
        //in fact everythin past the first space is considered description rather than name
        pRet->sName = sLine.substr(1, sLine.find(' '));
        while(pFile->good() && !pFile->eof() && pFile->peek() != '+')
        {
            std::getline (*pFile, sLine);
            size_t uiLineSize = len(sLine);
            std::vector<uint8_t> xQuality(uiLineSize, 126);//uiLineSize uint8_t's with value 127
            pRet->vAppend((const uint8_t*)sLine.c_str(), xQuality.data(), uiLineSize);
        }//while
        pRet->vTranslateToNumericFormUsingTable(pRet->xNucleotideTranslationTable, 0);
        //quality
        unsigned int uiPos = 0;
        while(pFile->good() && !pFile->eof() && pFile->peek() != '@')
        {
            std::getline (*pFile, sLine);
            size_t uiLineSize = len(sLine);
            for(size_t i=0; i < uiLineSize; i++)
                pRet->quality(i + uiPos) = (uint8_t)sLine[i];
            uiPos += uiLineSize;
        }//while
        return pRet;
    }//if
    //if we reach this point we have read all content of the file
    throw ModuleDryException();
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