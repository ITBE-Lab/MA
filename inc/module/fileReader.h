/** 
 * @file fileReader.h
 * @brief Reads queries from a file.
 * @author Markus Schmidt
 */
#ifndef FILE_READER_H
#define FILE_READER_H

#include "module/module.h"
#include "container/nucSeq.h"

/// @cond DOXYGEN_SHOW_SYSTEM_INCLUDES
#include <fstream>
/// @endcond

namespace libMA
{

    /**
     * @brief Reads Queries from a file.
     * @details
     * Reads (multi-)fasta or fastaq format.
     */
    class FileReader: public Module
    {
    private:
        class BufferedReader
        {
            private:
                const size_t uiBufferSize = 1048576; // == 2^20 ~= 0.1 GB buffer size
                std::vector<char> vBuffer;
                std::vector<>
            public:

        }; // class
    public:
        std::shared_ptr<std::ifstream> pFile;

        /**
         * @brief creates a new FileReader.
         */
        FileReader(std::string sFileName)
                :
            pFile(new std::ifstream(sFileName))
        {
            if (!pFile->is_open())
            {
                throw AlignerException("Unable to open file" + sFileName);
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

        // @override
        std::string getName() const
        {
            return "FileReader";
        }//function

        // @override
        std::string getFullDesc() const
        {
            return std::string("FileReader");
        }//function

        bool outputsVolatile() const
        {
            return true;
        }//function

    };//class

}//namespace

#ifdef WITH_PYTHON
void exportFileReader();
#endif

#endif