#ifndef PAIRED_READS_H
#define PAIRED_READS_H

namespace libMABS
{

    /**
     * @brief Picks paired alignments for paired reads
     * @ingroup module
     * @details
     * Given two vector of alignments this module picks one alignment of each vector, 
     * so that their combined score and penalty for the insert size is maximal.
     */
    class PairedReads: public Module
    {
    public:
        /**
         * @brief Penalty for unpaired reads.
         * @details
         * the score penalty that is substracted from the SW scores of the alignment of two *
         * unpaired reads.
         * If the penalty is smaller than the penalty added by the insert size 
         * between the paired seed we make an unpaired alignment.
         * The default value is taken from BWA-MEM as we use the same formula for pairing.
         */
        double u = 17;
        ///@brief use normal distribution for the insert size
        bool bNormalDist = true;
        ///@brief use uniform distribution for the insert size
        bool bUniformDist = false;
        ///@brief the mean of the insert size
        double mean = 0;
        ///@brief the standard deviation of the insert size
        double std = 1;

        /**
         * @brief The probability for a insert size >= d.
         * @details
         * Can be set to use normal or uniform distribution.
         * Note: both distributions are set using mean and standard deviation.
         */
        double p(nucSeqIndex d) const;

        PairedReads()
        {}//constructor

        std::shared_ptr<Container> execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - ContainerVector(Alignment)
         * - ContainerVector(Alignment)
         */
        ContainerVector getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * ContainerVector(Alignment)
         * with two elements
         */
        std::shared_ptr<Container> getOutputType() const;

        std::string getName() const
        {
            return "PairedReads";
        }
    };//class
}//namspace libMABS

/**
 * @brief export the PairedReads and PairedAlignment @ref Module "module" to python.
 * @ingroup export
 */
void exportPairedReads();



#endif