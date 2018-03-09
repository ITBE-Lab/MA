/** 
 * @file fileWriter.h
 * @brief Writes alignments to a file.
 * @author Markus Schmidt
 */
#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include "module/module.h"
#include "container/nucSeq.h"
#include "container/alignment.h"
#include "container/pack.h"
#include <fstream>

namespace libMA
{
    /**
     * @brief wrapper for various out streams.
     * @details
     * this exists, so that alignments can be written to a GUI element if so desired.
     */
    class OutStream
    {
    public:
        virtual OutStream& operator<<(std::string) {return *this;};
    };//class

    /**
     * @brief wraps the std outstream.
     */
    class StdOutStream : public OutStream
    {
    public:
        StdOutStream& operator<<(std::string s)
        {
            std::cout << s;
            return *this;
        }//function
    };//class

    /**
     * @brief wraps a file outstream.
     * @details
     * truncates the file if it already exists
     */
    class FileOutStream : public OutStream
    {
    public:
        std::ofstream file;

        FileOutStream(std::string sFileName)
            :
            file(sFileName, std::ofstream::out | std::ofstream::trunc)
        {
            if (!file.good())
            {
                std::cout << "Unable to open file\n";
                ///@todo exception
            }//if
        }//constructor

        ~FileOutStream()
        {
            file.close();
        }//deconstructor

        FileOutStream& operator<<(std::string s)
        {
            file << s;
            return *this;
        }//function
    };//class

    /**
     * @brief Writes SAM output.
     */
    class FileWriter: public Module
    {
    public:
        #define MULTIPLE_SEGMENTS_IN_TEMPLATE        0x001
        #define SEGMENT_PROPERLY_ALIGNED             0x002
        #define SEGMENT_UNMAPPED                     0x004
        #define NEXT_SEGMENT_UNMAPPED                0x008
        #define REVERSE_COMPLEMENTED                 0x010
        #define NEXT_REVERSE_COMPLEMENTED            0x020
        #define FIRST_IN_TEMPLATE                    0x040
        #define LAST_IN_TEMPLATE                     0x080
        #define SECONDARY_ALIGNMENT                  0x100
        #define NOT_PASSING_FILTERS                  0x200
        #define PCR_OR_DUPLICATE                     0x400
        #define SUPPLEMENTARY_ALIGNMENT              0x800

        //holds a file ourstream if necessary
        std::shared_ptr<OutStream> pOut;
        std::shared_ptr<std::mutex> pLock;

        /**
         * @brief creates a new FileWriter.
         * @details 
         * if sFileName is "stdout" the writer will output to stdout instead of the file.
         * Otherwise sFileName is used as the filename to write to.
         * The file will be truncated is it already exists.
         */
		FileWriter(std::string sFileName)
			:
			pLock(new std::mutex)
		{
			if (sFileName != "stdout")
                pOut = std::shared_ptr<OutStream>(new FileOutStream(sFileName));
            else
                pOut = std::shared_ptr<OutStream>(new StdOutStream());
			*pOut << "@HD VN:1.5 SO:unknown\n";
		}//constructor

        /**
         * @brief creates a new FileWriter.
         * @details 
         * Allows more control of the output by using a given OutStream.
         */
		FileWriter(std::shared_ptr<OutStream> pOut)
			:
            pOut(pOut),
			pLock(new std::mutex)
		{
			*pOut << "@HD VN:1.5 SO:unknown\n";
		}//constructor

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
            return "FileWriter";
        }//function

    };//class

    /**
     * @brief Writes human readable alignment output.
     * @details this format is loosely based on blast.
     */
    class RadableFileWriter: public Module
    {
    public:
        //holds a file ourstream if necessary
        std::shared_ptr<OutStream> pOut;
        std::shared_ptr<std::mutex> pLock;
        unsigned int uiNucsPerLine = 80;

        /**
         * @brief creates a new RadableFileWriter.
         */
		RadableFileWriter(std::shared_ptr<OutStream> pOut)
			:
            pOut(pOut),
			pLock(new std::mutex)
		{}//constructor

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
            return "ReadableFileWriter";
        }//function

    };//class

}//namespace

void exportFileWriter();

#endif