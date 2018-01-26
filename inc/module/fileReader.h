#ifndef FILE_READER_H
#define FILE_READER_H

#include "module/module.h"
#include "container/nucSeq.h"
#include <fstream>

namespace libMA
{

    class FileReader: public Module
    {
    public:
        std::shared_ptr<std::ifstream> pFile;

        FileReader(std::string sFileName)
                :
            pFile(new std::ifstream(sFileName))
        {
            if (!pFile->is_open())
            {
                std::cout << "Unable to open file" << std::endl;
                //@todo exception
            }//if
        }//constructor

        ~FileReader()
        {
            pFile->close();
        }//deconstructor

        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Nil
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(NucSeq)
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "FileReader";
        }//function

        bool outputsVolatile() const
        {
            return true;
        }//function

    };//class

}//namespace

void exportFileReader();

#endif