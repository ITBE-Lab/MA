/** 
 * @file pairedReads.cpp
 * @author Markus Schmidt
 */
#include "module/pairedReads.h"


using namespace libMA;

extern int iMatch;

#define SQ_12 std::sqrt(12)

double PairedReads::p(nucSeqIndex d) const
{
    if (bNormalDist)
        return 1 - 0.5 * (  1 + std::erf( (d - mean) / (std * std::sqrt(2)) )  );
    if (bUniformDist)
        /*
         * uniform distribution has following CDF:
         * (d - a) / (b - a)
         * where a and b are min,max.
         * the mean is thus: (b+a)/2
         * the std is: (b-a)*sqrt(12)
         */
        return 1 - std::max(  0.0, std::min(1.0, ( d - mean - std/(SQ_12*2) ) / ( std/SQ_12 ) )  );

    return 0.0;
}//function

ContainerVector PairedReads::getInputType() const
{
    return ContainerVector{
        std::shared_ptr<ContainerVector>(
            new ContainerVector(std::shared_ptr<Container>(new Alignment()))),
        std::shared_ptr<ContainerVector>(
            new ContainerVector(std::shared_ptr<Container>(new Alignment())))
    };
}//function

std::shared_ptr<Container> PairedReads::getOutputType() const
{
    return std::shared_ptr<Container>(new ContainerVector(std::shared_ptr<Alignment>(new Alignment())));
}//function

std::shared_ptr<Container> PairedReads::execute(
        std::shared_ptr<ContainerVector> vpInput
    )
{
    std::shared_ptr<ContainerVector> pAlignments1 = std::static_pointer_cast<ContainerVector>((*vpInput)[0]);
    std::shared_ptr<ContainerVector> pAlignments2 = std::static_pointer_cast<ContainerVector>((*vpInput)[1]);

    //the indices of the best alignment pair
    unsigned int uiI1 = 0, uiI2 = 0;
    double maxScore = 0;
    bool bPaired = false;
    // try out all combinations
    // this assumes that not more than three or four alignments are reported
    // otherwise this here might be a serious bottleneck
    for(unsigned int i = 0; i < pAlignments1->size(); i++)
    {
        std::shared_ptr<Alignment> pAlignment1 = std::static_pointer_cast<Alignment>((*pAlignments1)[i]);
        for(unsigned int j = 0; j < pAlignments2->size(); j++)
        {
            std::shared_ptr<Alignment> pAlignment2 = std::static_pointer_cast<Alignment>((*pAlignments2)[j]);

            // get the distance of the alignments on the reference
            nucSeqIndex d = pAlignment1->beginOnRef() - pAlignment2->beginOnRef();
            // make sure that we do not have underflows here
            if (pAlignment2->beginOnRef() > pAlignment1->beginOnRef())
                d = pAlignment2->beginOnRef() - pAlignment1->beginOnRef();

            // compute the score by the formula given in:
            // "Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM" 
            // (Heng Li)
            double score = pAlignment1->score() + pAlignment2->score() - std::min(-iMatch * std::log(p(d))/std::log(4), u);

            // check if we have a new best pair
            if(score > maxScore)
            {
                // if so update the score and indices
                maxScore = score;
                uiI1 = i;
                uiI2 = j;
                // it's possible that the best pair is two individual alignments
                bPaired = -iMatch * std::log(p(d))/std::log(4) <= u;
            }//if
        }//for
    }//for

    // set the paired property in the respective alignment stats
    std::static_pointer_cast<Alignment>((*pAlignments1)[uiI1])->xStats.bPaired = bPaired;
    std::static_pointer_cast<Alignment>((*pAlignments2)[uiI2])->xStats.bPaired = bPaired;
    //set which read was first...
    std::static_pointer_cast<Alignment>((*pAlignments1)[uiI1])->xStats.bFirst = true;
    std::static_pointer_cast<Alignment>((*pAlignments2)[uiI2])->xStats.bFirst = false;

    // return the best pair
    return std::shared_ptr<ContainerVector>(new ContainerVector
            {
                (*pAlignments1)[uiI1],
                (*pAlignments2)[uiI2]
            }//container vector
        );
}//function

void exportPairedReads()
{
    //export the PairedReads class
    boost::python::class_<
            PairedReads, 
            boost::python::bases<Module>, 
            std::shared_ptr<PairedReads>
        >("PairedReads")
    ;
    //tell boost that these two classes are convertible
    boost::python::implicitly_convertible< 
        std::shared_ptr<PairedReads>,
        std::shared_ptr<Module> 
    >();
}//function
