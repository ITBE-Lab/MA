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
    /*
     * Has next and next require synchronized access.
     * This is done by the module synchronization.
     */
    if(pFile->hasNext())
        return pFile->next();

    //if we reach this point we have read all content of the file
    return Nil::pEoFContainer;
}//function

#ifdef WITH_PYTHON
void exportFileReader()
{
    //export the FileReader class
    boost::python::class_<
            FileReader, 
            boost::python::bases<Module>, 
            std::shared_ptr<FileReader>
        >("FileReader", boost::python::init<std::string>())
    DEBUG(
        .def("testBufReader", &FileReader::testBufReader)
        .staticmethod("testBufReader")
    )
    ;

    boost::python::implicitly_convertible< 
        std::shared_ptr<FileReader>,
        std::shared_ptr<Module> 
    >();

}//function
#endif