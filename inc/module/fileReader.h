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
                const size_t uiFileBufferSize = 1048576; // == 2^20 ~= 0.1 GB buffer size
                const size_t uiQueryBufferSize = 100;
                unsigned int uiNucSeqBufPosRead = 0;
                unsigned int uiNucSeqBufPosWrite = 0;
                unsigned int uiCharBufPosRead = uiFileBufferSize;
                std::vector<char> vBuffer;
                std::vector<std::shared_ptr<NucSeq>> vpNucSeqBuffer;
                std::mutex xMutex;
                std::condition_variable cv;
                std::ifstream xFile;
                std::thread xThread;
                inline void reFillBuffer()
                {
                    xFile.read(vBuffer.data(), uiFileBufferSize);
                    uiCharBufPosRead = 0;
                }// function

                inline unsigned int searchEndline()
                {
                    unsigned int uiI = 0;
                    while(
                            uiCharBufPosRead + uiI < uiFileBufferSize &&
                            vBuffer[uiCharBufPosRead + uiI] != '\n'
                          )
                        uiI++;
                    return uiI;
                }// function
            public:

                BufferedReader(std::string sFileName)
                        :
                    vpNucSeqBuffer(),
                    xFile(sFileName)
                {
                    if (!xFile.is_open())
                    {
                        throw AlignerException("Unable to open file" + sFileName);
                    }//if
                    vpNucSeqBuffer.reserve(uiQueryBufferSize);
                    vBuffer.reserve(uiFileBufferSize);
                    // the async file reader
                    xThread = std::thread(
                        [&]
                        ()
                        {
                            std::unique_lock<std::mutex> xLock(xMutex);
                            while(pFile->good() && !pFile->eof())
                            {
                                std::shared_ptr<NucSeq> pCurr;
                                if(uiCharBufPosRead >= uiFileBufferSize)
                                    reFillBuffer();
                                // find end of name
                                assert(vBuffer[uiCharBufPosRead] == '>');
                                do
                                {
                                    unsigned int uiCharBufPosReadLen = searchEndline();
                                    for(
                                            unsigned int uiI = uiCharBufPosRead;
                                            uiI < uiCharBufPosRead + uiCharBufPosReadLen;
                                            uiI++
                                        )
                                        pCurr->sName += vBuffer[uiI];
                                    uiCharBufPosRead += uiCharBufPosReadLen;
                                    if(uiCharBufPosRead >= uiFileBufferSize)
                                        reFillBuffer();
                                }// do
                                while(uiCharBufPosRead >= uiFileBufferSize);

                                // remove the description from the query name
                                pCurr->sName = pCurr->sName.substr(1, pCurr->sName.find(' '));

                                // find end of nuc section
                                do
                                {
                                    unsigned int uiCharBufPosReadLen = searchEndline();
                                    //memcpy the data over
                                    pCurr->vAppend(
                                            &vBuffer[uiCharBufPosRead],
                                            uiCharBufPosReadLen
                                        );
                                    uiCharBufPosRead += uiCharBufPosReadLen;
                                    if(uiCharBufPosRead >= uiFileBufferSize)
                                        reFillBuffer();
                                }// do
                                while(
                                        uiCharBufPosRead >= uiFileBufferSize ||
                                        vBuffer[uiCharBufPosRead] != '>'
                                    );
                                // @todo write into buffer and wait on mutex

                                pCurr->vTranslateToNumericFormUsingTable(
                                        pCurr->xNucleotideTranslationTable,
                                        0
                                    );
                                vpNucSeqBuffer[uiNucSeqBufPosWrite] = pCurr;
                                while(
                                        ( (uiNucSeqBufPosWrite + 1) % uiQueryBufferSize) 
                                            !=
                                        uiNucSeqBufPosRead
                                    )
                                    cv.wait(xLock);
                                uiNucSeqBufPosWrite = (uiNucSeqBufPosWrite + 1) % uiQueryBufferSize;
                            }// while
                        }// lambda
                    );// std::thread
                }// constructor

                ~BufferedReader()
                {
                    xThread.join();
                }// deconstructor

                bool hasNext()
                {
                    if(uiNucSeqBufPosRead != uiNucSeqBufPosWrite)
                        return true;

                    // no data immediately available -> we have to wait
                    std::lock_guard<std::mutex> xLock(xMutex);
                    return uiNucSeqBufPosRead != uiNucSeqBufPosWrite;
                }// function

                std::shared_ptr<NucSeq> next()
                {
                    auto pQuery = vpNucSeqBuffer[uiNucSeqBufPosRead];
                    vpNucSeqBuffer[uiNucSeqBufPosRead] = nullptr;
                    uiNucSeqBufPosRead = (uiNucSeqBufPosRead + 1) % uiQueryBufferSize;
                    cv.notify_one();
                    return pQuery;
                }// function
        };// class
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

        void testBufReader()
        {

        }// function

    };//class

}//namespace

#ifdef WITH_PYTHON
void exportFileReader();
#endif

#endif