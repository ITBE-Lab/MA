/** 
 * @file export.cpp
 * @author Markus Schmidt
 */
#include "util/export.h"

using namespace libMA;


extern int iGap;
extern int iExtend;
extern int iMatch;
extern int iMissMatch;

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
        exportMinimizers();
        exportMinimizersHash();
        exportOtherSeeding();
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
    unsigned int uiReportNBest,
    bool bLocal,
    int iMatch_,
    int iMisMatch_,
    int iGap_,
    int iExtend_,
    unsigned int uiMaxGapArea
)
{
    iMatch = iMatch_;
    iExtend = iExtend_;
    iGap = iGap_;
    iMissMatch = iMisMatch_;

    if(uiNumSOC < uiReportNBest)
        throw AlignerException("cannot report more alignments than computed (increase strip of consideration amount)");

    //setup all modules

    //modules required for any alignment
    std::shared_ptr<Module> pLockQuery(new Lock(std::shared_ptr<Container>(new NucSeq())));
    std::shared_ptr<Module> pSeeding(new BinarySeeding(bSeedSetPairs));
    std::shared_ptr<Module> pSOC(new StripOfConsideration(uiMaxAmbiguity, uiNumSOC));
    std::shared_ptr<Module> pCouple(
            new ExecOnVec(std::shared_ptr<Module>(new LinearLineSweep(uiMaxGapArea)))
        );
    //we only want to report the best alignment
    std::shared_ptr<Module> pDoOptimal(new ExecOnVec(
        // false here is counter intuitive by experiments did show that soting by soc score is better...
        std::shared_ptr<Module>(new NeedlemanWunsch(bLocal)), false, uiReportNBest));
    std::shared_ptr<Module> pMapping(new MappingQuality());

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
        //the optimal matching stage
        std::shared_ptr<Pledge> pOptimal = Module::promiseMe(
                    pDoOptimal, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pCoupled,
                        pQuery,
                        pPack
                    }
                );
        //assign a mapping quality
        std::shared_ptr<Pledge> pAlignments = Module::promiseMe(
                    pMapping, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pQuery,
                        pOptimal
                    }
                );
        if(bPaired)
        {
            if(aQueries.size() != 2)
                throw AlignerException("two input files are required for paired alignments");
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
            std::shared_ptr<Pledge> pOptimal = Module::promiseMe(
                    pDoOptimal, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pCoupled2,
                        pQuery2,
                        pPack
                    }
                );
            //the optimal matching stage
            std::shared_ptr<Pledge> pAlignments2 = Module::promiseMe(
                    pMapping, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pQuery2,
                        pOptimal
                    }
                );
            //pick the best alignments
            pAlignments = Module::promiseMe(
                    pPaired, 
                    std::vector<std::shared_ptr<Pledge>>
                    {
                        pAlignments,
                        pAlignments2
                    }
                );
            
        }//if
        else
        {
            if(aQueries.size() != 1)
                throw AlignerException("one input files is required for unpaired alignments");
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