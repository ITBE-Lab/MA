/** 
 * @file binarySeeding.h
 * @brief Implements a segmentation algorithm.
 * @author Markus Schmidt
 */
#ifndef BINARY_SEEDING_H
#define BINARY_SEEDING_H

#include "util/system.h"
#include "module/module.h"
#include "container/segment.h"

namespace libMA
{
    class PerfectMatch;
    
    #define complement(x) (uint8_t)NucSeq::nucleotideComplement(x)

    /**
     * @brief Computes a maximally covering set of seeds.
     * @details
     * Can use either the extension scheme by Li et Al. or ours.
     * @ingroup module
     */
    class BinarySeeding : public Module{
    public:
        bool bLrExtension;
        const unsigned int uiMinAmbiguity;

        /**
         * @brief The simplified extension scheme presented in our Paper.
         * @details
         * Computes Two segments for each index as follows:
         * - extend backwards first then forwards
         * - extend forwards first then backwards
         * Starts both extensions on center.
         * Returns an interval spanning the entire covered area.
         * Segments are saved in pSegmentVector.
         */
        inline Interval<nucSeqIndex> maximallySpanningExtension(
                nucSeqIndex center,
                std::shared_ptr<FMIndex> pFM_index,
                std::shared_ptr<NucSeq> pQuerySeq,
                std::shared_ptr<SegmentVector> pSegmentVector
            )
        {
            // query sequence itself 
            const uint8_t *q = pQuerySeq->pGetSequenceRef(); 
            
            /* Initialize ik on the foundation of the single base q[x].
            * In order to understand this initialization you should have a look 
            * to the corresponding PowerPoint slide.
            */
            // start I(q[x]) in T (start in BWT used for backward search) + 1, 
            // because very first string in SA-array starts with $
            // size in T and T' is equal due to symmetry
            SAInterval ik(
                                pFM_index->L2[complement(q[center])] + 1, 
                                pFM_index->L2[(int)q[center]] + 1, 
                                pFM_index->L2[(int)q[center] + 1] - pFM_index->L2[(int)q[center]]
                            );

            /*
            * extend ik right, until there are no more matches
            */
            nucSeqIndex end = center;
            for(nucSeqIndex i = center+1; i < pQuerySeq->length(); i++)
            {
                DEBUG_3(
                    std::cout << i-1 << " -> " << ik.start() << " " << ik.end() << std::endl;
                    std::cout << i-1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
                )
                assert(ik.size() > 0);
                SAInterval ok = pFM_index->extend_backward(ik, complement(q[i]));

                DEBUG_3(
                    std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
                    std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
                )
                /*
                * In fact, if ok.getSize is zero, then there are no matches any more.
                */
                if (ok.size() <= uiMinAmbiguity)
                    break; // the SA-index interval size is too small to be extended further
                end = i;
                ik = ok;
            }//for
            DEBUG_3(
                std::cout << "swap" << std::endl;
            )
            //this is required in order to extend the other way
            ik = ik.revComp();
            nucSeqIndex start = center;
            /*
            * extend ik left, until there are no more matches
            */
            if(center > 0)
            {
                for(nucSeqIndex i = center-1; i >= 0; i--)
                {
                    DEBUG_3(
                        std::cout << i+1 << " -> " << ik.start() << " " << ik.end() << std::endl;
                        std::cout << i+1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
                    )
                    assert(ik.size() > 0);
                    SAInterval ok = pFM_index->extend_backward(ik, q[i]);
                    DEBUG_3(
                        std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
                        std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
                    )

                    /*
                    * In fact, if ok.getSize is zero, then there are no matches any more.
                    */
                    if (ok.size() <= uiMinAmbiguity)
                        break; // the SA-index interval size is too small to be extended further
                    start = i;
                    ik = ok;
                    //cause nuxSeqIndex is unsigned
                    if(i == 0)
                        break;
                }//for
            }//if
            std::shared_ptr<Segment> pRightLeft(new Segment(start,end-start,ik));
            assert(start >= 0);
            assert(end < pQuerySeq->length());
            assert(pRightLeft->end() < pQuerySeq->length());
            pSegmentVector->push_back(pRightLeft);
            DEBUG_3(
                std::cout << "--other way--" << std::endl;
            )
            /* Initialize ik on the foundation of the single base q[x].
            * In order to understand this initialization you should have a look 
            *to the corresponding PowerPoint slide.
            */
            // start I(q[x]) in T (start in BWT used for backward search) + 1, 
            // because very first string in SA-array starts with $
            // size in T and T' is equal due to symmetry
            ik = SAInterval(
                                pFM_index->L2[q[center]] + 1, 
                                pFM_index->L2[(int)complement(q[center])] + 1, 
                                pFM_index->L2[(int)q[center] + 1] - pFM_index->L2[(int)q[center]]
                            );
            start = center;
            /*
            * extend ik left, until there are no more matches
            */
            if(center > 0)
            {
                for(nucSeqIndex i = center-1; i >= 0; i--)
                {
                    DEBUG_3(
                        std::cout << i+1 << " -> " << ik.start() << " " << ik.end() << std::endl;
                        std::cout << i+1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
                    )
                    assert(ik.size() > 0);
                    SAInterval ok = pFM_index->extend_backward(ik, q[i]);
                    DEBUG_3(
                        std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
                        std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
                    )

                    /*
                    * In fact, if ok.getSize is zero, then there are no matches any more.
                    */
                    if (ok.size() <= uiMinAmbiguity)
                        break; // the SA-index interval size is too small to be extended further
                    start = i;
                    ik = ok;
                    //cause nuxSeqIndex is unsigned
                    if(i == 0)
                        break;
                }//for
            }//if
            DEBUG_3(
                std::cout << "swap" << std::endl;
            )
            //this is required in order to extend the other way
            ik = ik.revComp();
            end = center;
            /*
            * extend ik right, until there are no more matches
            */
            for(nucSeqIndex i = center+1; i < pQuerySeq->length(); i++)
            {
                DEBUG_3(
                    std::cout << i-1 << " -> " << ik.start() << " " << ik.end() << std::endl;
                    std::cout << i-1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
                )
                assert(ik.size() > 0);
                SAInterval ok = pFM_index->extend_backward(ik, complement(q[i]));

                DEBUG_3(
                    std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
                    std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
                )

                /*
                * In fact, if ok.getSize is zero, then there are no matches any more.
                */
                if (ok.size() <= uiMinAmbiguity)
                    break; // the SA-index interval size is too small to be extended further
                end = i;
                ik = ok;
            }//for
            std::shared_ptr<Segment> pLeftRight(new Segment(start,end-start,ik.revComp()));
            assert(start >= 0);
            assert(end < pQuerySeq->length());
            assert(pLeftRight->end() < pQuerySeq->length());
            pSegmentVector->push_back(pLeftRight);

            //to return the covered area
            Interval<nucSeqIndex> ret(center,0);
            if(pLeftRight->start() < pRightLeft->start())
                ret.start(pLeftRight->start());
            else
                ret.start(pRightLeft->start());

            if(pLeftRight->end() > pRightLeft->end())
                ret.end(pLeftRight->end());
            else
                ret.end(pRightLeft->end());

            return ret;
        }//function

        /**
         * @brief The extension scheme from Li et Al.
         * @details
         * Computes all non-enclosed segments overlapping center.
         * Returns an interval spanning the entire covered area.
         * Segments are saved in pSegmentVector.
         */
        Interval<nucSeqIndex> smemExtension(
                nucSeqIndex center,
                std::shared_ptr<FMIndex> pFM_index,
                std::shared_ptr<NucSeq> pQuerySeq,
                std::shared_ptr<SegmentVector> pSegmentVector
            )
            {
                //to remember the covered area
                Interval<nucSeqIndex> ret(center,0);

                // query sequence itself
                const uint8_t *q = pQuerySeq->pGetSequenceRef(); 

                assert(center < pQuerySeq->length());
                
                /* Initialize ik on the foundation of the single base q[x].
                * In order to understand this initialization you should have a look 
                *to the corresponding PowerPoint slide.
                */
                // start I(q[x]) in T (start in BWT used for backward search) + 1, 
                // because very first string in SA-array starts with $
                // size in T and T' is equal due to symmetry
                SAInterval ik(
                                    pFM_index->L2[complement(q[center])] + 1, 
                                    pFM_index->L2[(int)q[center]] + 1, 
                                    pFM_index->L2[(int)q[center] + 1] - pFM_index->L2[(int)q[center]]
                                );

                /*
                * forward extension first
                * this way we need to swap only once (forward to backwards) instead of swapping
                * (backwards to forwards to backwards)
                * 
                * curr is used to remember the Suffix array interval each time we loose some hits by extending
                */
                std::vector<Segment> curr;
                // extend until the end of the query
                for(nucSeqIndex i = center+1; i < pQuerySeq->length(); i++)
                {
                    DEBUG_3(
                        std::cout << i-1 << " -> " << ik.start() << " " << ik.end() << std::endl;
                        std::cout << i-1 << " ~> " << ik.revComp().start() << " " << ik.revComp().end() << std::endl;
                    )
                    assert(ik.size() > 0);
                    //this is the extension
                    SAInterval ok = pFM_index->extend_backward(ik, complement(q[i]));

                    // checking weather we lost some intervals
                    // if so -> remember the interval just before we lost the hits
                    if(ok.size() != ik.size()) 
                        // save the reverse complement cause when extending the saved interval we will extend
                        // in the other direction
                        curr.push_back(Segment(center, i-center-1, ik.revComp()));
                    // if were at the end of the query and we still have hits we need to make sure to record them
                    if(i == pQuerySeq->length()-1 && ok.size() != 0)
                        // save the reverse complement cause when extending the saved interval we will extend
                        // in the other direction
                        curr.push_back(Segment(center, i-center, ok.revComp()));

                    DEBUG_3(
                        std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
                        std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
                    )
                    /*
                    * In fact, if ok.getSize is zero, then there are no matches any more.
                    * thus we can stop extending forwards
                    */
                    if (ok.size() <= uiMinAmbiguity)
                        break; // the SA-index interval size is too small to be extended further
                    // if we get here we can forget the old interval and save the current interval.
                    ik = ok;
                    // remember that we covered this area
                    ret.end(i);
                }//for
                DEBUG_3(
                    std::cout << "swap" << std::endl;
                )
                /*
                * This is the backwards extension part
                * Here we need to extend the intervals in reverse order with respect to how we discovered them.
                * (reversing is done by push_front insted of push_back)
                *
                * we will use prev and curr in this way:
                *         each iteration we will extend all intervals in prev
                *         and save the intervals that need to be extended further in curr
                *         at the end of the iteration we will swap prev and curr
                *         then clear curr
                */
                std::reverse(curr.begin(), curr.end());
                std::vector<Segment> prev;
                //pointers for easy swapping of the lists

                // FIXME: for some reason valgrind does NOT like this
                // maybe it cant deal with the pointers?
                std::vector<Segment> *pPrev, *pCurr, *pTemp;
                pPrev = &curr;
                pCurr = &prev;
                // quick check that we can extend backwards at all (center is unsigned thus this is necessary)
                if(center != 0)
                {
                    // extend until we reach the start of the query
                    for(nucSeqIndex i = center-1; i >= 0; i--)
                    {
                        assert(pCurr->empty());

                        /*
                        * we need to remember weather finished extending some interval in this step.
                        * because:
                        *         if we already have found one with this length
                        *             then all following intervals that we find have to be enclosed
                        *             (this is due to the fact that they we know they start further right but 
                        *              end at the same point)
                        */
                        bool bHaveOne = false;

                        /*
                        * for all remembered intervals 
                        * (ordered by the start on the query)
                        */
                        for(Segment& ik : *pPrev)
                        {
                            DEBUG_3(
                                std::cout << i+1 << " -> " << ik.saInterval().start() << " " << ik.saInterval().end() << std::endl;
                                std::cout << i+1 << " ~> " << ik.saInterval().revComp().start() << " " << ik.saInterval().revComp().end() << std::endl;
                            )
                            // actually extend the current interval
                            SAInterval ok = pFM_index->extend_backward(ik.saInterval(), q[i]);
                            DEBUG_3(
                                std::cout << i << " -> " << ok.start() << " " << ok.end() << std::endl;
                                std::cout << i << " ~> " << ok.revComp().start() << " " << ok.revComp().end() << std::endl;
                            )
                            DEBUG_3(
                                std::cout << ik.start() << ", " << ik.end() << ": " << ik.saInterval().size() << " -> " << ok.size() << std::endl;
                            )
                            // check if the extension resulted in a non enclosed interval
                            if(ok.size() <= uiMinAmbiguity && !bHaveOne)
                            {
                                // save the interval
                                pSegmentVector->push_back(std::shared_ptr<Segment>(new Segment(ik)));
                                assert(ik.start() <= ik.end());
                                assert(ik.end() <= pQuerySeq->length());
                                // we need to remember that we already found a interval this iteration
                                bHaveOne = true;
                            }// if
                            // check if we can extend this interval further
                            else if(ok.size() > uiMinAmbiguity)
                            {
                                // if so add the intervals to the list
                                Segment xSeg = Segment(i, ik.size()+1, ok);
                                // FIXME: memory leak here according to valgrind ?!?
                                pCurr->push_back(xSeg);
                                assert(xSeg.end() <= pQuerySeq->length());
                            }// if
                        }// for


                        // swap out the lists and clear the things we just worked on
                        pTemp = pPrev;
                        pPrev = pCurr;
                        pCurr = pTemp;
                        // FIXME: memory leak here according to valgrind ?!?
                        pCurr->clear();
                        pCurr->shrink_to_fit();

                        // if there are no more intervals to extend
                        if(pPrev->empty())
                            break;

                        // remember that we covered this area
                        ret.start(i);

                        // cause nuxSeqIndex is unsigned we have to avoid underflow
                        if(i == 0)
                            break;
                    }// for
                }// if

                //if we reach the beginning of the query it is possible that there are still intervals that contain matches.
                //we need to save the longest of those, which is conveniently (due to our sorting) the first one in the list
                if(!pPrev->empty())
                {
                    assert(pPrev->front().size() >= pPrev->back().size());

                    pSegmentVector->push_back(std::shared_ptr<Segment>(new Segment(pPrev->front())));
                    assert(pPrev->front().start() <= pPrev->front().end());
                    assert(pPrev->front().end() <= pQuerySeq->length());
                    
                    DEBUG_2(
                        std::cout << pPrev->front().start() << ":" << pPrev->front().end() << std::endl;
                    )
                }//if

                //return the area that we covered
                return ret;
            }//function

        /*
        *    does nothing if the given interval can be found entirely on the genome.
        *    if the interval cannot be found this method splits the interval in half and repeats the step with the first half,
        *    while queuing the second half as a task in the thread pool.
        */
        void procesInterval(
                Interval<nucSeqIndex> xAreaToCover,
                std::shared_ptr<SegmentVector> pSegmentVector,
                std::shared_ptr<FMIndex> pFM_index,
                std::shared_ptr<NucSeq> pQuerySeq
            );

    public:
        /**
         * @brief Initialize a BinarySeeding Module
         * @details
         * if bLrExtension is True our extension scheme is used,
         * otherwise the extension scheme by Li et Al. is used.
         * Our approach is faster and computes seeds of higher quality.
         * However Li et Al.s approach will increase the overall accuracy of the alignment.
         */
        BinarySeeding(bool bLrExtension, unsigned int uiMinAmbiguity)
                :
            bLrExtension(bLrExtension),
            uiMinAmbiguity(uiMinAmbiguity)
        {}//constructor
        
        std::shared_ptr<Container> EXPORTED execute(std::shared_ptr<ContainerVector> vpInput);

        /**
         * @brief Used to check the input of execute.
         * @details
         * Returns:
         * - FMIndex
         * - NucSeq
         */
        ContainerVector EXPORTED getInputType() const;

        /**
         * @brief Used to check the output of execute.
         * @details
         * Returns:
         * - SegmentVector
         */
        std::shared_ptr<Container> EXPORTED getOutputType() const;

        std::string getName() const
        {
            return "BinarySeeding";
        }

        std::string getFullDesc() const
        {
            return std::string("BinarySeeding(") + 
                std::to_string(bLrExtension) + "," +
                std::to_string(uiMinAmbiguity) + ")";
        }//function
    };//class

}//namespace


/**
 * @brief exports the Segmentation @ref Module "module" to python.
 * @ingroup export
 */
void exportBinarySeeding();

#endif