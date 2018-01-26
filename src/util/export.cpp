#include "util/export.h"

using namespace libMA;

/**
 * @brief The boost-python main method.
 *
 * A main function that exposes all Containers and Modules to python
 * by calling the respective export*() functions.
 * The main method is generated using boost-python.
 */
BOOST_PYTHON_MODULE(libMA)
{
        DEBUG_3(
                std::cout.setf(std::ios::unitbuf);
        )
        exportContainer();
        exportModule();
        exportFM_index();
        exportSequence();
        exportBinarySeeding();
        exportPack();
        exportIntervalTree();
        exportExceptions();
        exportSeed();
        exportAlignment();
        exportLinesweep();
        exportNeedlemanWunsch();
        exportStripOfConsideration();
        exportChaining();
        exportSMW();
        exportExtractAllSeeds();
        exportExecOnVector();
        exportReSeed();
        exportFileReader();
        exportFileWriter();
        exportMappingQuality();
        exportPairedReads();
        exportSplitter();
}//function

std::vector<std::shared_ptr<Pledge>> setUpCompGraph(
    std::shared_ptr<Pledge> pPack,
    std::shared_ptr<Pledge> pFMDIndex,
    std::vector<std::shared_ptr<Pledge>> aQueries,
    std::shared_ptr<Module> pOut,
    unsigned int uiThreads,
    unsigned int uiMaxAmbiguity,
    unsigned int uiNumSOC,
    bool bPariedNormal,
    bool bPariedUniform,
    unsigned int uiPairedMean,
    double fPairedStd,
    double dPairedU,
    bool bSeedSetPairs,
    unsigned int uiReportNBest
)
{
    if(uiNumSOC < uiReportNBest)
        throw new AlignerException("cannot report more alignments than computed (increase strip of consideration amount)");

    //setup all modules

    //modules required for any alignment
    std::shared_ptr<Module> pLockQuery(new Lock(std::shared_ptr<Container>(new NucSeq())));
    std::shared_ptr<Module> pSeeding(new BinarySeeding(bSeedSetPairs));
    std::shared_ptr<Module> pSOC(new StripOfConsideration(uiMaxAmbiguity, 0, uiNumSOC, 0.0f, 0.0f));
    std::shared_ptr<Module> pCouple(new ExecOnVec(std::shared_ptr<Module>(new LinearLineSweep())));
    //we only want to report the best alignment
    std::shared_ptr<Module> pOptimal(new ExecOnVec(
        std::shared_ptr<Module>(new NeedlemanWunsch(0)), true, uiReportNBest));

    //modules for the paired alignment
    bool bPaired = bPariedNormal || bPariedUniform;
    std::shared_ptr<Module> pPaired(new PairedReads(
        dPairedU, 
        bPariedNormal, 
        bPariedUniform, 
        uiPairedMean, 
        fPairedStd
    ));

    //setup the computational graph
    std::vector<std::shared_ptr<Pledge>> aRet;
    for(unsigned int i=0; i < uiThreads; i++)
    {
        //lock the query in for this subgraph
        std::shared_ptr<Pledge> pQuery = Module::promiseMe(
                pLockQuery, 
                std::vector<std::shared_ptr<Pledge>>
                {
                    aQueries[0]
                }
            );
        //one more module to unlock the locked-in query
        //this requires the pledge from the lock
        //therefore we have to create one module for each subgraph
        std::shared_ptr<Module> pUnLock(new UnLock(pQuery));
        //the seeding stage
        std::shared_ptr<Pledge> pSeeds = Module::promiseMe(
                pSeeding, 
                std::vector<std::shared_ptr<Pledge>>
                {
                    pFMDIndex,
                    pQuery
                }
            );
        //the filtering stage
        std::shared_ptr<Pledge> pSOCs = Module::promiseMe(
                pSOC, 
                std::vector<std::shared_ptr<Pledge>>
                {
                    pSeeds,
                    pQuery,
                    pPack,
                    pFMDIndex
                }
            );
        //the coupling stage
        std::shared_ptr<Pledge> pCoupled = Module::promiseMe(
                pCouple, 
                std::vector<std::shared_ptr<Pledge>>
                {
                    pSOCs
                }
            );
        std::shared_ptr<Pledge> pAlignments;
        if(bPaired)
        {
            if(aQueries.size() != 2)
                throw new AlignerException("two input files are required for paired alignments");
            //lock the query in for this subgraph
            std::shared_ptr<Pledge> pQuery2 = Module::promiseMe(
                    pLockQuery, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        aQueries[1]
                    }
                );
            //synchronize both locks
            //this automatically makes it so that unlock unlocks both
            Pledge::synchronize(pQuery, pQuery2);
            //the seeding stage
            std::shared_ptr<Pledge> pSeeds2 = Module::promiseMe(
                    pSeeding, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pFMDIndex,
                        pQuery2
                    }
                );
            //the filtering stage
            std::shared_ptr<Pledge> pSOCs2 = Module::promiseMe(
                    pSOC, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pSeeds2,
                        pQuery2,
                        pPack,
                        pFMDIndex
                    }
                );
            //the coupling stage
            std::shared_ptr<Pledge> pCoupled2 = Module::promiseMe(
                    pCouple, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pSOCs2
                    }
                );
            //the optimal matching stage
            std::shared_ptr<Pledge> pAlignments1 = Module::promiseMe(
                    pOptimal, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pCoupled,
                        pQuery,
                        pPack
                    }
                );
            //the optimal matching stage
            std::shared_ptr<Pledge> pAlignments2 = Module::promiseMe(
                    pOptimal, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pCoupled2,
                        pQuery2,
                        pPack
                    }
                );
            //pick the best alignments
            pAlignments = Module::promiseMe(
                    pPaired, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pAlignments1,
                        pAlignments2
                    }
                );
            
        }//if
        else
        {
            if(aQueries.size() != 1)
                throw new AlignerException("one input files is required for unpaired alignments");
            //the optimal matching stage
            pAlignments = Module::promiseMe(
                    pOptimal, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pCoupled,
                        pQuery,
                        pPack
                    }
                );
        }//else

        //write the output to a file
        std::shared_ptr<Pledge> pNil = Module::promiseMe(
                pOut, 
                std::vector<std::shared_ptr<Pledge>>
                {
                    pQuery,
                    pAlignments,
                    pPack
                }
            );
        //unlock the query so that this subgraph can be executed multiple times
        std::shared_ptr<Pledge> pRet = Module::promiseMe(
                pUnLock, 
                std::vector<std::shared_ptr<Pledge>>
                {
                    pNil
                }
            );
        //save the 
        aRet.push_back(pRet);
    }//for

    return aRet;
}//function