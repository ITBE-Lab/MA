#ifndef FILE_READER_H
#define FILE_READER_H

#include "module/module.h"
#include "container/nucSeq.h"
#include <fstream>

namespace libMABS
{

    class FileReader: public Module
    {
    public:
        std::string sFileName;

        FileReader(std::string sFileName)
                :
            sFileName(sFileName)
        {}//constructor

        void getLines(std::function<void(std::string&)> fDo);

        std::shared_ptr<Container> execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - Nil
         */
        ContainerVector getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - ContainerVector(NucSeq)
         */
        std::shared_ptr<Container> getOutputType() const;

        std::string getName() const
        {
            return "FileReader";
        }//function

    };//class

}//namespace

void exportFileReader();

#endif