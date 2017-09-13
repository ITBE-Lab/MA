
#define NUM_THREADS_ALIGNER 6

#ifndef GRAPHICAL_METHOD_H
#define GRAPHICAL_METHOD_H

#define NUM_THREADS_ALIGNER 6

#define DEBUG_ENABLED

#include "intervalTree.h"
#include <algorithm>
#include "segmentation.h"
#include "meta_programming.h"
#include "threadPool.h"
#include "module.h"
#include <boost/python.hpp>
#include <memory>


#if 0
class PerfectMatch;//TODO: replace perfect match by seed

/* container for the perfect matches.
* extends the match data structure by an id and a disabled flag.
* the container is necessary since every match will be stored in 2 buckets each.
* the buckets might get processed simultaneously and the disabled flag might be different for each combination of match and bucket.
*/
struct PerfectMatchContainer{
	//the actual perfect match from the segmentation step
	std::shared_ptr<const PerfectMatch> pxPerfectMatch;
	//the bucket wide unique id for this match can be used to access the match in the match array of the according bucket
	//once the matches are sorted though this id gets invalidated
	size_t uiId;
	//is the match enabled in this bucket?
	bool bDisabled;
};

/* each perfect match "casts a shadow" at the left and right border of the bucket
 * each shadow is stored in one of these data structures.
*/
class ShadowInterval: public Interval<nucSeqIndex>{
	//is this shadow cast at the left or at the right border of the bucket?
	bool bLeft;
	//the bucket wide unique id for the corresponding match. it can be used to access the match in the match array of the according bucket
	//once the matches are sorted though this id gets invalidated
	uint16_t uiCorrespondingPerfectMatchId;
	//while swiping the interfering shadows will get stored in this list
	std::list<ShadowInterval*> lpxInterferingIntervals;
};

/* a perfect match calculated in the segmentation process
*/
class PerfectMatch{
private:
	//length of the match
	const nucSeqIndex uiLength;
	//the beginning of the match on the reference
	const nucSeqIndex uiPosOnReference;
	//the beginning of the match on the query
	const nucSeqIndex uiPosOnQuery;

public:
	PerfectMatch(const nucSeqIndex uiLenght, const nucSeqIndex uiPosOnReference, const nucSeqIndex uiPosOnQuery)
		:
		uiLength(uiLenght),
		uiPosOnReference(uiPosOnReference),
		uiPosOnQuery(uiPosOnQuery)
	{}//constructor

	//add the length of this match to some number
	nucSeqIndex operator+(const nucSeqIndex uiX) const { return this->uiLength + uiX; }//function
	//in order to determine the bucket this match should be positioned in we need to subtract the position of this match on the query from the position on the reference
	//since we work with unsigned int's we will add the length of the query in order to ensure positive numbers.
	nucSeqIndex getPositionForBucketing(nucSeqIndex uiQueryLength) const { return uiPosOnReference + (uiQueryLength - uiPosOnQuery); }//function

	/*print the position of this match relative to the beginning of the bucket.
	*/
	void printAsChainView(std::ofstream &xOut, const nucSeqIndex uiBucketBegin, const nucSeqIndex uiQueryLenght, const nucSeqIndex uiBucketSize) const
	{

		/* we want to print the distance between the position of the match on the reference and the beginning of the bucket on the reference
		* in order to make it possible to work with unsigned int's the buckets are offset by the query length 
		* (in order to find the corresponding bucket we subtract the position of the match on the query from the position on the reference.
		* this could become negative so we add the length of the query).
		*/
		assert(uiPosOnReference + uiQueryLenght + uiBucketSize >= uiBucketBegin);
		xOut << "(" << uiPosOnReference + uiQueryLenght + uiBucketSize - uiBucketBegin << "," << uiPosOnQuery << "," << uiLength << ")";
	}//function

	/*determine the start and end positions this match casts on the left border of the given bucket
	  pxMatch is the container this match is stored in.
	*/
	ShadowInterval getLeftShadow(nucSeqIndex uiBucketStart, PerfectMatchContainer pxMatch, nucSeqIndex uiBucketSize, nucSeqIndex uiQueryLength) const
	{
		ShadowInterval xRet;
		xRet.uiStart = uiPosOnQuery + uiBucketStart;
		xRet.uiEnd = uiPosOnReference + uiLength + uiQueryLength + uiBucketSize;
		xRet.bLeft = true;
		xRet.uiCorrespondingPerfectMatchId = pxMatch.uiId;
		return xRet;
	}

	/*determine the start and end positions this match casts on the right border of the given bucket
	pxMatch is the container this match is stored in.
	*/
	ShadowInterval getRightShadow(nucSeqIndex uiBucketStart, PerfectMatchContainer pxMatch, nucSeqIndex uiBucketSize, nucSeqIndex uiQueryLength) const
	{
		ShadowInterval xRet;
		xRet.uiStart = uiPosOnReference;
		xRet.uiEnd = uiPosOnQuery + uiLength + uiBucketStart + uiQueryLength + uiBucketSize * 2;
		xRet.bLeft = false;
		xRet.uiCorrespondingPerfectMatchId = pxMatch.uiId;
		return xRet;
	}

	//getters
	nucSeqIndex getLength() const { return uiLength; }//function
	nucSeqIndex getPosOnQuery() const { return uiPosOnQuery; }//function
	nucSeqIndex getPosOnReference() const { return uiPosOnReference; }//function
};
#endif

class SeedBucket
{
private:
	nucSeqIndex uiTotalScore;

	std::list<std::shared_ptr<Seed>> lContent;

	std::mutex xMutex;

	//bool bUsed;

public:
	SeedBucket()
		:
		uiTotalScore(0),
		lContent(),
		xMutex()
		//,bUsed(false)
	{}//constructor

	////disable copying of buckets
	SeedBucket(const SeedBucket&) = delete;

	void addSeed(const std::shared_ptr<Seed> pNew)
	{
		std::lock_guard<std::mutex> xGuard(xMutex);
		lContent.push_back(pNew);
		uiTotalScore += pNew->size();
		//end of scope xGuard
	}//function

	nucSeqIndex getValue()
	{
		return uiTotalScore;
	}//function

	void forall(std::function<void(const std::shared_ptr<Seed>)> fDo)
	{
		for (auto pSeed : lContent)
		{
			fDo(pSeed);
		}//for
	}//function

	//bool isUsed() const { return bUsed; }//function
	//void setUsed() { bUsed = true; }//function
};//class

class StripOfConsideration : public Container, public Interval<nucSeqIndex>
{
private:

	/*the summed up value of the content.
	* this value gets modified while processing the bucket.
	* for example we might detect that two matches contradict each other, then we will disable one of them and subtract it's value
	* but initially this is just the plain sum of all matches lying in this bucket
	*/
	nucSeqIndex uiValueOfContent;

	//contains all matches that belong into this bucket
	std::list<std::shared_ptr<Seed>> lSeeds;

#if 0
	/*"cast" the "shadows" of all matches against the left and right border of the bucket, will store the outcome in axShadows
	*/
	void calculateShadows(nucSeqIndex uiQueryLength)
	{
		axShadows.clear();
		//each match will "cast" 2 "shadows"
		axShadows.resize(apxPerfectMatches.size() * 2);

		for (unsigned int uiI = 0; uiI < apxPerfectMatches.size(); uiI++)
		{
			axShadows[uiI*2] = apxPerfectMatches[uiI].pxPerfectMatch->getLeftShadow(uiBegin, apxPerfectMatches[uiI], uiEnd - uiBegin, uiQueryLength);
			axShadows[uiI * 2 + 1] = apxPerfectMatches[uiI].pxPerfectMatch->getRightShadow(uiBegin, apxPerfectMatches[uiI], uiEnd - uiBegin, uiQueryLength);
		}//for
	}//function

	/*sorts the shadows by their beginnings (increasingly). if 2 matches start at the exact same position they will be sorted by their end (decreasingly)
	*/
	void sortShadows()
	{

		if (axShadows.size() <= 1)
			return;

		//sort (increasingly) by start coordinate of the match
		std::sort(axShadows.begin(), axShadows.end(),
			[](const ShadowInterval xA, const ShadowInterval xB)
			{
				assert(xA.bLeft != xB.bLeft || xA.uiCorrespondingPerfectMatchId != xB.uiCorrespondingPerfectMatchId);

				if (xA.uiStart == xB.uiStart)
					return xA.uiEnd > xB.uiEnd;
				return xA.uiStart < xB.uiStart;
			}//lambda
		);//sort function call

		assert(axShadows.begin()->uiStart < (--axShadows.end())->uiStart ||
			(axShadows.begin()->uiStart == (--axShadows.end())->uiStart && axShadows.begin()->uiEnd >= (--axShadows.end())->uiEnd)
			);
	}//function

	/* calculate the exact score of the content in this bucket. 
	 * (considering disabled matches and overlapping matches)
	 * the array of perfect matches will be sorted therefore all perfect match ids will get invalidated
	 *
	 * !!!will invalidate all perfect match ids!!!
	 */
	void calculateScore()
	{
		if (apxPerfectMatches.size() <= 1)
			return;

		//sort (increasingly) by start coordinate of the match
		std::sort(apxPerfectMatches.begin(), apxPerfectMatches.end(),
			[](const PerfectMatchContainer xA, const PerfectMatchContainer xB)
			{
				return xA.pxPerfectMatch->getPosOnQuery() < xB.pxPerfectMatch->getPosOnQuery();
			}//lambda
		);//sort function call

		/*what we do:
		*	walk over the matches which are sorted by their beginnings.
		*	always remember the end of the last match.
		*	check if the current match overlaps with the last one (either on the query or on the reference)
		*	if so subtract the overlapped area from the overall score
		*/

		/*
		* nucSeqIndex getCombinedArea(ShadowInterval *pxCurr) contains very similar code to this one, but for the shadows
		*/
		nucSeqIndex uiLastPosOnQuery = 0;
		nucSeqIndex uiLastPosOnRef = 0;
		nucSeqIndex uiLastLength = 0;

		uiValueOfContent = 0;
		for (unsigned int uiI = 0; uiI < apxPerfectMatches.size(); uiI++)
		{

			if (!apxPerfectMatches[uiI].bDisabled)
			{

				nucSeqIndex uiThisPosOnQuery = apxPerfectMatches[uiI].pxPerfectMatch->getPosOnQuery();
				nucSeqIndex uiThisPosOnRef = apxPerfectMatches[uiI].pxPerfectMatch->getPosOnReference();
				nucSeqIndex uiThisLength = apxPerfectMatches[uiI].pxPerfectMatch->getLength();

				uiValueOfContent += apxPerfectMatches[uiI].pxPerfectMatch->getLength();
				if (uiLastLength != 0)
				{
					nucSeqIndex uiQueryOverlap = (uiLastPosOnQuery + uiLastLength) - uiThisPosOnQuery;
					if (uiThisPosOnQuery > uiLastLength + uiLastPosOnQuery)
						uiQueryOverlap = 0;

					nucSeqIndex uiRefOverlap = (uiLastPosOnRef + uiLastLength) - uiThisPosOnRef;
					if (uiThisPosOnRef > uiLastLength + uiLastPosOnRef)
						uiRefOverlap = 0;


					if (uiQueryOverlap > uiRefOverlap)
						uiValueOfContent -= uiQueryOverlap;
					else
						uiValueOfContent -= uiRefOverlap;

#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
					if (uiQueryOverlap > uiRefOverlap)
						uiAmountOfOverlapsInAligment += uiQueryOverlap;
					else
						uiAmountOfOverlapsInAligment += uiRefOverlap;
					if (uiThisPosOnQuery > uiLastLength + uiLastPosOnQuery && uiThisPosOnRef > uiLastLength + uiLastPosOnRef)
					{
						uiAmountOfGapsInAligment++;
						if (uiThisPosOnQuery - (uiLastLength + uiLastPosOnQuery) > uiLongestGapInALignment)
							uiLongestGapInALignment = uiThisPosOnQuery - (uiLastLength + uiLastPosOnQuery);
						if (uiThisPosOnRef - (uiLastLength + uiLastPosOnRef) > uiLongestGapInALignment)
							uiLongestGapInALignment = uiThisPosOnRef - (uiLastLength + uiLastPosOnRef);
					}
#endif

				}//if

				uiLastLength = uiThisLength;
				uiLastPosOnQuery = uiThisPosOnQuery;
				uiLastPosOnRef = uiThisPosOnRef;
			}//if
		}//for

	}//function


	nucSeqIndex getCombinedArea(ShadowInterval *pxCurr)
	{
		nucSeqIndex uiRet = 0;

		/*what we do:
		*	walk over the shadows which are sorted by their beginnings.
		*	add the length of the match (corresponding to the current shadow) to the score
		*	always remember the end of the last match.
		*	check if the current match overlaps with the last one (either on the query or on the reference)
		*	if so subtract the overlapped area from the score
		*/

		/*
		* calculateScore() contains very similar code to this one, but for the actual matches
		*/
		nucSeqIndex uiLastPosOnQuery = 0;
		nucSeqIndex uiLastPosOnRef = 0;
		nucSeqIndex uiLastLength = 0;
		for (auto pxInterval : pxCurr->lpxInterferingIntervals)
		{
			if (!apxPerfectMatches[pxInterval->uiCorrespondingPerfectMatchId].bDisabled)
			{
				nucSeqIndex uiThisPosOnQuery = apxPerfectMatches[pxInterval->uiCorrespondingPerfectMatchId].pxPerfectMatch->getPosOnQuery();
				nucSeqIndex uiThisPosOnRef = apxPerfectMatches[pxInterval->uiCorrespondingPerfectMatchId].pxPerfectMatch->getPosOnReference();
				nucSeqIndex uiThisLength = apxPerfectMatches[pxInterval->uiCorrespondingPerfectMatchId].pxPerfectMatch->getLength();

				uiRet += apxPerfectMatches[pxInterval->uiCorrespondingPerfectMatchId].pxPerfectMatch->getLength();
				if (uiLastLength != 0)
				{
					nucSeqIndex uiQueryOverlap = (uiLastPosOnQuery + uiLastLength) - uiThisPosOnQuery;
					if (uiThisPosOnQuery > uiLastLength + uiLastPosOnQuery)
						uiQueryOverlap = 0;

					nucSeqIndex uiRefOverlap = (uiLastPosOnRef + uiLastLength) - uiThisPosOnRef;
					if (uiThisPosOnRef > uiLastLength + uiLastPosOnRef)
						uiRefOverlap = 0;


					if (uiQueryOverlap > uiRefOverlap)
						uiRet -= uiQueryOverlap;
					else
						uiRet -= uiRefOverlap;
				}//if

				uiLastLength = uiThisLength;
				uiLastPosOnQuery = uiThisPosOnQuery;
				uiLastPosOnRef = uiThisPosOnRef;
			}//if
			//end if
		}//for

		return uiRet;
	}//function

	/*helper function for  lineSweep() 
	* this will deal with all intervals that we have passed.
	* for all passed intervals it will:
	*	iterate over all matches that contradict the passed one and decide if the passed or all contradicting matches need to be disabled (using the total length).
	*
	* this only works since the lpxShadowIntervalsSurroundingCurrPos is sorted by the end of the shadows and because we add the shadows in the order of their starting point.
	*	therefore we ensure that we reach all shadows that are completely within some other shadow before that one.
	*	thus we ensure that the decision to keep the smaller intervals with respect to even smaller intervals within them is already made when making the decision of weather to keep the big or a the small intervals.
	*	so this is a greedy approach that makes the right decision in every step, as long as the order of the decisions is correct. (this is enshured by linesweep())
	*/
	void disableIntervalsBefore(std::list<ShadowInterval*> &lpxShadowIntervalsSurroundingCurrPos, nucSeqIndex uiWhere)
	{
		auto ppxSurroundingCurrPosItt = lpxShadowIntervalsSurroundingCurrPos.begin();
		assert(ppxSurroundingCurrPosItt == lpxShadowIntervalsSurroundingCurrPos.end() || uiWhere >= (*ppxSurroundingCurrPosItt)->uiStart);
		//remove all intervals at the beginning of the list that we passed when jumping to the beginning of this interval
		while (ppxSurroundingCurrPosItt != lpxShadowIntervalsSurroundingCurrPos.end() && (*ppxSurroundingCurrPosItt)->uiEnd < uiWhere)
		{
			//do not consider disabled matches
			if (apxPerfectMatches[(*ppxSurroundingCurrPosItt)->uiCorrespondingPerfectMatchId].bDisabled)
			{
				lpxShadowIntervalsSurroundingCurrPos.erase(ppxSurroundingCurrPosItt++);
				continue;
			}//if

			//get the summed length of all interfering intervals. (consider that some of those might overlap)
			nucSeqIndex uiSizeInterferingIntervals = getCombinedArea(*ppxSurroundingCurrPosItt);

			//decide to disable this or all interfering intervals based on the length
			if (uiSizeInterferingIntervals < apxPerfectMatches[(*ppxSurroundingCurrPosItt)->uiCorrespondingPerfectMatchId].pxPerfectMatch->getLength())
			{
				for (auto pxInterval : (*ppxSurroundingCurrPosItt)->lpxInterferingIntervals)
					if (!apxPerfectMatches[pxInterval->uiCorrespondingPerfectMatchId].bDisabled)
					{
						disableMatch(pxInterval->uiCorrespondingPerfectMatchId);
						pxInterval->lpxInterferingIntervals.clear();
					}//if
				//end for
			}//if
			else
			{
				disableMatch((*ppxSurroundingCurrPosItt)->uiCorrespondingPerfectMatchId);
				(*ppxSurroundingCurrPosItt)->lpxInterferingIntervals.clear();
			}//else

			//erase this interval from the list. we never have to visit it again
			lpxShadowIntervalsSurroundingCurrPos.erase(ppxSurroundingCurrPosItt++);
		}//while
	}//function

	/*helper function for lineSweep()
	* given an shadow and the list of shadows surrounding the current position this will add markers to all shadows that completely surround the given shadow.
	* then it will insert the current shadow at the correct position of the list (the list is sorted by the end of the shadows)
	* the function expects the list to be in a valid state considering the given interval.
	* this is ensured by lineSweep() using disableIntervalsBefore()
	*/
	void addInvalidMarkers(std::list<ShadowInterval*> &lpxShadowIntervalsSurroundingCurrPos, ShadowInterval *pxCurr)
	{
		//add edges between the current and all elements in the list that end after the current.
		auto ppxSurroundingCurrPosItt = lpxShadowIntervalsSurroundingCurrPos.end();
		while (ppxSurroundingCurrPosItt != lpxShadowIntervalsSurroundingCurrPos.begin() && (*--ppxSurroundingCurrPosItt)->uiEnd >= pxCurr->uiEnd)
		{
			assert((*ppxSurroundingCurrPosItt)->bLeft == pxCurr->bLeft);

			//ignore already disabled intervals
			if (apxPerfectMatches[(*ppxSurroundingCurrPosItt)->uiCorrespondingPerfectMatchId].bDisabled)
				continue;

			assert(pxCurr->uiCorrespondingPerfectMatchId != (*ppxSurroundingCurrPosItt)->uiCorrespondingPerfectMatchId);

			(*ppxSurroundingCurrPosItt)->lpxInterferingIntervals.emplace_back(pxCurr);

			//iterator gets decreased in the second condition of the while loop
		}//while

		//increment the iterator because insert will insert BEFORE the current element, but we want to insert after.
		//there a two special cases where we do not want to increment:
		// 1) the iterator currently points to the past-the-end position (e.g. if the list was empty from the beginning)
		// 2) the iterator points the the very first element and the element we want to insert has to go before the first one
		if (ppxSurroundingCurrPosItt != lpxShadowIntervalsSurroundingCurrPos.end()
			&&
			!(ppxSurroundingCurrPosItt == lpxShadowIntervalsSurroundingCurrPos.begin() && (*ppxSurroundingCurrPosItt)->uiEnd >= pxCurr->uiEnd)
			)
			ppxSurroundingCurrPosItt++;
		//end if
		lpxShadowIntervalsSurroundingCurrPos.insert(ppxSurroundingCurrPosItt, pxCurr);
	}//function

	/*this function is executing the line sweep algorithm in order to remove matches that contradict each other
	* in a way that the matches adding up to the maximal possible sum remain.
	* therefore the shadows cast by the matches against the left & right wall of the current bucket are used
	* two helper functions are used disableIntervalsBefore() , addInvalidMarkers()
	* the function expects calculateShadows() and sortShadows() to be called beforehand.
	* check the according ppt slides to understand the algorithm.
	*/
	void lineSweep()
	{
		//we distinguish between left and right shadows for performance reasons
		std::list<ShadowInterval*> lpxShadowIntervalsSurroundingCurrPosL;
		std::list<ShadowInterval*> lpxShadowIntervalsSurroundingCurrPosR;


		auto pxCurrShadowIntervalIterator = axShadows.begin();
		//iterate over all shadows
		while (pxCurrShadowIntervalIterator != axShadows.end())
		{
			//ignore disabled ones
			if (apxPerfectMatches[pxCurrShadowIntervalIterator->uiCorrespondingPerfectMatchId].bDisabled)
			{
				pxCurrShadowIntervalIterator++;
				continue;
			}//if

			//remove all shadows which ends have been passed and make the decision which ones to keep
			disableIntervalsBefore(
				pxCurrShadowIntervalIterator->bLeft ? lpxShadowIntervalsSurroundingCurrPosL : lpxShadowIntervalsSurroundingCurrPosR
				, pxCurrShadowIntervalIterator->uiStart
				);

			//add a marker in all shadows that surround the current one
			// &* converts from iterator to pointer
			addInvalidMarkers(
				pxCurrShadowIntervalIterator->bLeft ? lpxShadowIntervalsSurroundingCurrPosL : lpxShadowIntervalsSurroundingCurrPosR
				, &*pxCurrShadowIntervalIterator
				);

			//enable this block if you want to check lpxShadowIntervalsSurroundingCurrPos for consistency in each iteration
#if 0
			/*check the list for consistency*/
			if (lpxShadowIntervalsSurroundingCurrPos.size() > 0)
			{
				auto ppxSurroundingCurrPosItt = lpxShadowIntervalsSurroundingCurrPos.begin();
				auto pxDragpointer = *ppxSurroundingCurrPosItt;
				if (ppxSurroundingCurrPosItt != lpxShadowIntervalsSurroundingCurrPos.end())
					ppxSurroundingCurrPosItt++;
				while (ppxSurroundingCurrPosItt != lpxShadowIntervalsSurroundingCurrPos.end())
				{
					assert(pxDragpointer->uiEnd <= (*ppxSurroundingCurrPosItt)->uiEnd);
					pxDragpointer = *ppxSurroundingCurrPosItt;
					ppxSurroundingCurrPosItt++;
				}//for
			}//if
			/*end check list for consistency*/
#endif


			pxCurrShadowIntervalIterator++;
		}//while

		//remove all remaining shadows and make the decision which ones to keep
		if (lpxShadowIntervalsSurroundingCurrPosL.size() != 0)
		{
			auto endOfList = *(--lpxShadowIntervalsSurroundingCurrPosL.end());
			disableIntervalsBefore(lpxShadowIntervalsSurroundingCurrPosL, endOfList->uiEnd + 1);
		}//if

		//remove all remaining shadows and make the decision which ones to keep
		if (lpxShadowIntervalsSurroundingCurrPosR.size() != 0)
		{
			auto endOfList = *(--lpxShadowIntervalsSurroundingCurrPosR.end());
			disableIntervalsBefore(lpxShadowIntervalsSurroundingCurrPosR, endOfList->uiEnd + 1);
		}//if
	}//function

	/*this function will check weather any of the matches in the current bucket contradict each other.
	* it works much like the actual linesweep algorithm. check linesweep() in order to understand whats going on.
	*/
	bool noShadowsOverlap()
	{
		std::list<ShadowInterval*> lpxShadowIntervalsSurroundingCurrPosL;
		std::list<ShadowInterval*> lpxShadowIntervalsSurroundingCurrPosR;

		auto pxCurrShadowIntervalIterator = axShadows.begin();
		while (pxCurrShadowIntervalIterator != axShadows.end())
		{
			if (apxPerfectMatches[pxCurrShadowIntervalIterator->uiCorrespondingPerfectMatchId].bDisabled)
			{
				pxCurrShadowIntervalIterator++;
				continue;
			}//if

			auto lpxShadowIntervalsSurroundingCurrPos = (pxCurrShadowIntervalIterator->bLeft ? lpxShadowIntervalsSurroundingCurrPosL : lpxShadowIntervalsSurroundingCurrPosR);
			auto ppxSurroundingCurrPosItt = lpxShadowIntervalsSurroundingCurrPos.end();
			if (ppxSurroundingCurrPosItt != lpxShadowIntervalsSurroundingCurrPos.begin() && (*--ppxSurroundingCurrPosItt)->uiEnd >= pxCurrShadowIntervalIterator->uiEnd)
			{
				if ((*ppxSurroundingCurrPosItt)->bLeft != pxCurrShadowIntervalIterator->bLeft)
					return false;

				if (apxPerfectMatches[(*ppxSurroundingCurrPosItt)->uiCorrespondingPerfectMatchId].bDisabled)
					continue;

				if (pxCurrShadowIntervalIterator->uiCorrespondingPerfectMatchId == (*ppxSurroundingCurrPosItt)->uiCorrespondingPerfectMatchId)
					return false;

				//iterator gets decreased in the second condition of the while loop
			}//while
			if (ppxSurroundingCurrPosItt != lpxShadowIntervalsSurroundingCurrPos.end()
				&&
				!(ppxSurroundingCurrPosItt == lpxShadowIntervalsSurroundingCurrPos.begin() && (*ppxSurroundingCurrPosItt)->uiEnd >= pxCurrShadowIntervalIterator->uiEnd)
				)
				ppxSurroundingCurrPosItt++;
			//end if
			lpxShadowIntervalsSurroundingCurrPos.insert(ppxSurroundingCurrPosItt, &*pxCurrShadowIntervalIterator);

			pxCurrShadowIntervalIterator++;
		}//while

		return true;
	}//function

	#endif

public:

	StripOfConsideration()
			:
		Interval(0, 0),
		lSeeds()
	{}//constructor

	StripOfConsideration(nucSeqIndex uiStart, nucSeqIndex uiSize)
		:
		Interval(uiStart, uiSize),
		lSeeds()
	{}//constructor

	
	/*used to identify the strip of cinsideration datatype in the aligner pipeline*/
	ContainerType getType(){return ContainerType::stripOfConsideration;}
	
	void addElement(SeedBucket& rxBucket)
	{
		rxBucket.forall(
			[&](const std::shared_ptr<Seed> pSeed)
			{
				lSeeds.push_back(pSeed);
				uiValueOfContent += pSeed->getValue();
			}//lambda
		);
	}//function

	inline nucSeqIndex getValue() const
	{
		return uiValueOfContent;
	}//function

	inline void setValue(const nucSeqIndex uiNewVal)
	{
		uiValueOfContent = uiNewVal;
	}//function

	inline void subtractFromValue(const nucSeqIndex uiNewVal)
	{
		uiValueOfContent -= uiNewVal;
	}//function

	std::list<std::shared_ptr<Seed>>& seeds()
	{
		return lSeeds;
	}//function

	void forAllSeeds(std::function<void(std::list<std::shared_ptr<Seed>>::iterator pSeed)> fDo)
	{
		for(std::list<std::shared_ptr<Seed>>::iterator pSeed = lSeeds.begin(); pSeed != lSeeds.end(); pSeed++)
			fDo(pSeed);
	}//function
#if 0
	/*adds a match to the current bucket. 
	* will add the length of the match to the score of the bucket.
	* note that the sum of all matches is not a good score for the bucket, since 2 matches might contradict or overlap.
	* the sum should rather be seen as the maximal possible score for this bucket.
	* call process() to calculate the actual score.
	*/
	void addElement(std::shared_ptr<SeedBucket> pxNew)
	{
		if(pxNew == nullptr)
			return;
		uiValueOfContent += pxNew->getValue();
		lpxSeedBuckets.push_back(pxNew);
	}//function

	bool withinBucket(std::shared_ptr<const Seed> pxNew, nucSeqIndex uiQueryLength)
	{
		return uiBegin <= pxNew->getPositionForBucketing(uiQueryLength) && uiEnd >= pxNew->getPositionForBucketing(uiQueryLength);
	}//function

	/* this returns the score of the bucket.
	* before calling process() this will return the maximal possible score.
	* after calling process the actual score will be returned.
	* it might not be necessary to process each bucket (depending on the possible maximum score)
	*/
	nucSeqIndex getValueOfContet() const { return uiValueOfContent; }//function
	nucSeqIndex getBegin() { return this->uiBegin; }//function
	nucSeqIndex getEnd() { return this->uiEnd; }//function

	nucSeqIndex averagePosOnRef() 
	{
		nucSeqIndex uiRet = 0;
		nucSeqIndex uiCount = 0;
		for (unsigned int ui = 0; ui < apxPerfectMatches.size(); ui++)
		{
			if (!apxPerfectMatches[ui].bDisabled)
			{
				uiRet += apxPerfectMatches[ui].pxPerfectMatch->getPosOnReference();
				uiCount++;
			}//if
		}//for

		return uiRet / uiCount;
	}//function

	/*print the content of the bucket, relative to the location of the bucket
	*/
	void printAsChainView(std::ofstream &xOut, nucSeqIndex uiBucketSize, nucSeqIndex uiQueryLength) const {
		xOut << "(" << uiQueryLength << "," << uiBucketSize << ")";
		for (unsigned int ui = 0; ui < apxPerfectMatches.size(); ui++)
		{
			if (!apxPerfectMatches[ui].bDisabled)
				apxPerfectMatches[ui].pxPerfectMatch->printAsChainView(xOut, uiBegin, uiQueryLength, uiBucketSize);
		}//for
	}//function


	/*calculate the actual score of all matches within this bucket. 
	* this will label some matches as "disabled", meaning that they are not part of the optimal solution
	* once this is finished getValueOfContent() will return the actual score instead of the maximal possible one.
	*/
	void process(nucSeqIndex uiQueryLength)
	{

		//move all elements from the buckets into the vector
		for (auto &pxBucket : lpxPerfectMatcheBuckets)
		{
			pxBucket->forall([&](std::shared_ptr<const PerfectMatch> pxMatch)
				{
					if (withinBucket(pxMatch, uiQueryLength))
					{
						PerfectMatchContainer xCont = {
							pxMatch,
							(size_t)apxPerfectMatches.size(),
							false
						};
						apxPerfectMatches.push_back(xCont);
					}//if
					else
						uiValueOfContent -= pxMatch->getLength();
				}//lambda
			);//forall
		}//for
		lpxPerfectMatcheBuckets.clear();

		if (uiValueOfContent == 0)
			return;

		/*in order to run line sweep,
		* the shadows have to be created and sorted.
		*/
		calculateShadows(uiQueryLength);
		sortShadows();
		//actually run lineSweep()
		lineSweep();

		/*once lineSweep() is done getVaueOfContet() will actually return correct the score considering contradicting matches
		* but NOT considering overlaps
		* still it would be possible to break here and use the score as a good approximation 
		*/

		//quick check if the given solution is ok
		//NOTE: this does not check weather the given solution is the optimal one
		assert(noShadowsOverlap());

		//shadows are not required anymore
		axShadows.clear();

		//calculate the actual score, considering disabled AND overlapping matches
		calculateScore();
	}//function

	unsigned int numMatches()
	{
		return apxPerfectMatches.size();
	}//function

	std::shared_ptr<const PerfectMatch> getMatch(unsigned int at)
	{
		return apxPerfectMatches[at].pxPerfectMatch;
	}//function

	bool matchEnabled(unsigned int at)
	{
		return !apxPerfectMatches[at].bDisabled;
	}//function
#endif
};//class


class StripOfConsiderationVector: public Container
{//TODO: get rid of this wierd class
public:
	std::vector<std::shared_ptr<StripOfConsideration>> x;
	
	StripOfConsiderationVector(std::vector<std::shared_ptr<StripOfConsideration>> x)
			:
		x(x)
	{}//constructor
	
	StripOfConsiderationVector(std::shared_ptr<StripOfConsideration> x)
			:
		x(std::vector<std::shared_ptr<StripOfConsideration>>{x})
	{}//constructor

	StripOfConsiderationVector()
	{}//constructor

	/*used to identify the FM_indexWrapper datatype in the aligner pipeline*/
    ContainerType getType(){return ContainerType::stripOfConsiderationList;}
};//class

#if 0
/*given the perfect matches this can be used to find the matches that represent the optimal solution
* this is an alternative to chaining
*/
class GraphicalMethod
{
private:
	/*the matches are stored in buckets depending on their position on the reference - the position on the query
	* this speeds up the process significantly, since the lineSweep algorithm that is used requires ~O(n^2) 
	* (the sigma is rather in the direction of O(n) though, especially if we give the algorithm few contradicting matches. for no contradicting matches we have O(n))
	* the number of contradicting matches within each strip should be very small
	*/
	std::vector<std::shared_ptr<StripOfConsideration>> apxStripsOfConsideration;
	const nucSeqIndex uiReferenceLength;
	const nucSeqIndex uiQueryLength;
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
	AlignmentQuality *pxQuality;
#endif

public:

	GraphicalMethod(nucSeqIndex uiReferenceLength, nucSeqIndex uiQueryLength
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
		, AlignmentQuality *pxQuality
#endif
		)
		:
		apxStripsOfConsideration(),
		uiReferenceLength(uiReferenceLength),
		uiQueryLength(uiQueryLength)
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
		, pxQuality(pxQuality)
#endif
	{}//constructor

	void addStripOfConsideration(std::shared_ptr<StripOfConsideration> pxNew)
	{
		apxStripsOfConsideration.push_back(pxNew);
	}//function

	/*before processing, StripOfConsideration.getVaueOfContet() will return the maximal possible score
	* after being processed the bucket will return the exact score.
	* this processes each bucket
	*
	* when in doubt call process() or smartProcess() rather than this
	*
	* (processAllStrips() could be used instead of process() if the strips shall be printed afterwards)
	*/
	void processAllStrips()
	{
		ThreadPool xPool(NUM_THREADS_ALIGNER);

		unsigned int uiPorcessNbuckets = 1000;
		for (unsigned int ui = 0; ui < apxStripsOfConsideration.size(); ui += uiPorcessNbuckets)
		{
			xPool.enqueue(
				[/*WARNING: do not catch anything here: the lambda function is enqueued to into a thread pool,
				 local variables might not exist anymore during it's execution*/]
				(size_t, std::vector<std::shared_ptr<StripOfConsideration>> *papxBucket, unsigned int uiQueryLength, unsigned int uiPorcessNbuckets, unsigned int uiStartHere)
				{
					for (auto ui = uiStartHere; ui < uiPorcessNbuckets + uiStartHere && ui < papxBucket->size(); ui++)
						(*papxBucket)[ui]->process(uiQueryLength);

					//BOOST_LOG_TRIVIAL(info) << uiPorcessNbuckets + uiStartHere << " of " << paxBucket->size() << " buckets were processed";

				}, //lambda
				&apxStripsOfConsideration, uiQueryLength, uiPorcessNbuckets, ui
			);
		}//for

		//wait for the processing of the buckets to finish
		//end of scope xPool
	}//function

	/*this sorts the buckets according to their score
	* note that buckets will return their maximal possible score instead of their exact score before being processed.
	*/
	void sortStrips()
	{
		//sort (increasing) by the summed value of the content in each bucket.
		std::sort(apxStripsOfConsideration.begin(), apxStripsOfConsideration.end(),
			[](const std::shared_ptr<StripOfConsideration> pxA, const std::shared_ptr<StripOfConsideration> pxB)
			{
				return pxA->getValueOfContet() > pxB->getValueOfContet();
			}//lambda
		);
		assert(apxStripsOfConsideration[0]->getValueOfContet() >= apxStripsOfConsideration[apxStripsOfConsideration.size() - 1]->getValueOfContet());
	}//function

	/* this will process buckets until the maximal possible score of the best non processed bucket is smaller than the score of the best processed bucket.
	* this way we can be sure that we found the optimal bucket.
	* do NOT sort the buckets after calling this, the function will place the best bucket at the beginning of the bucket array.
	* this is way faster than processAllBuckets() since we can expect to find only very few positions on the reference where we can find most of the query.
	* therefore this is likely to not process any more than 10 buckets.
	*/
	void smartProcess()
	{
		if(apxStripsOfConsideration.empty())
			return;
		sortStrips();
		unsigned int uiBestBucket = 0;
		unsigned int uiProcessedBucketUntil = 0;
		//while the score of the best processed bucket is lower than the score of the best unprocessed bucket
		while (uiProcessedBucketUntil < apxStripsOfConsideration.size() && apxStripsOfConsideration[uiBestBucket]->getValueOfContet() <= apxStripsOfConsideration[uiProcessedBucketUntil]->getValueOfContet())
		{
			{
				ThreadPool xPool(NUM_THREADS_ALIGNER);
				//process #NUM_THREADS_ALIGNER buckets
				for (unsigned int ui = uiProcessedBucketUntil; ui < NUM_THREADS_ALIGNER + uiProcessedBucketUntil && ui < apxStripsOfConsideration.size(); ui++)
					xPool.enqueue(
					[/*WARNING: do not catch anything here: the lambda function is enqueued to into a thread pool,
					 local variables might not exist anymore during it's execution*/]
					(size_t, std::shared_ptr<StripOfConsideration> pxBucket, unsigned int uiQueryLength)
					{
						pxBucket->process(uiQueryLength);

					}, //lambda
					apxStripsOfConsideration[ui], uiQueryLength
					);
			}//end of scope xPool

			//check if we found a new best bucket
			for (unsigned int ui = uiProcessedBucketUntil; ui < NUM_THREADS_ALIGNER + uiProcessedBucketUntil && ui < apxStripsOfConsideration.size(); ui++)
				if (apxStripsOfConsideration[uiBestBucket]->getValueOfContet() < apxStripsOfConsideration[ui]->getValueOfContet())
					uiBestBucket = ui;

			uiProcessedBucketUntil += NUM_THREADS_ALIGNER;
		}//while

		//move the best bucket to the beginning of the array
		auto pxTemp = apxStripsOfConsideration[0];
		apxStripsOfConsideration[0] = apxStripsOfConsideration[uiBestBucket];
		apxStripsOfConsideration[uiBestBucket] = apxStripsOfConsideration[0];



#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
		pxQuality->uiLongestSeedInAlignment = apxStripsOfConsideration[0]->uiLongestSeedInALignment;
		pxQuality->uiLongestGapInAlignment = apxStripsOfConsideration[0]->uiLongestGapInALignment;
		pxQuality->uiAmountOfGapsInAligment = apxStripsOfConsideration[0]->uiAmountOfGapsInAligment;
		pxQuality->uiAmountOfSeedsInAligment = apxStripsOfConsideration[0]->uiAmountOfSeedsInAligment;
		pxQuality->uiAmountOfOverlapsInAligment = apxStripsOfConsideration[0]->uiAmountOfOverlapsInAligment;
#endif

	}//function
	/*before processing, Bucket.getVaueOfContet() will return the maximal possible score
	* after being processed the bucket will return the exact score.
	* this processes each bucket
	* and sorts the array of buckets in descending order
	*/
	void process()
	{
		processAllStrips();
		sortStrips();
	}//function

	/*will return the bucket with the n-th best score or nullptr if that bucket does not exist
	* you should process() or smartProcess() the buckets otherwise this will just return the n-th bucket (counting from the position on the reference)
	*/
	std::shared_ptr<StripOfConsideration> getNthBestBucket(unsigned int uiN)
	{
		if (uiN >= apxStripsOfConsideration.size())
			return nullptr;
		return apxStripsOfConsideration[uiN];
	}//function

	/*prints the value of the content of each bucket
	*/
	void printBuckets(std::ofstream &xOut, std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxSequencePack, const nucSeqIndex uiStripSize)
	{
		for (unsigned int ui = 0; ui < apxStripsOfConsideration.size(); ui++)
		{
			//do not print empty buckets
			if (apxStripsOfConsideration[ui]->getValueOfContet() == 0)
				continue;

			if (ui * uiStripSize >= uiQueryLength && pxSequencePack->bPositionIsOnReversStrand(ui * uiStripSize - uiQueryLength))
				xOut << "(" << uiReferenceLength - (ui + 1) * uiStripSize + 1 << ","
				<< uiReferenceLength - ui * uiStripSize
				<< "," << apxStripsOfConsideration[ui]->getValueOfContet() << ",reverse)";
			else
				xOut << "(" << ui * uiStripSize << "," << (ui + 1) * uiStripSize - 1 << "," << apxStripsOfConsideration[ui]->getValueOfContet() << ",forward)";
		}//for
		xOut << std::endl;
	}//function

};

class AnchorMatchList{
private:
	typedef std::vector<std::shared_ptr<PerfectMatchBucket>> PerfectMatchBuckets;
	typedef std::list<std::shared_ptr<PerfectMatch>> AnchorSegments;

	PerfectMatchBuckets apxPerfectMatchBuckets;
	AnchorSegments lpxAnchorMatches;
	nucSeqIndex uiStripSize; 
	nucSeqIndex uiQueryLength;
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
	AlignmentQuality *pxQuality;
#endif


	std::shared_ptr<StripOfConsideration> collectStripOfConsideration(std::shared_ptr<const PerfectMatch> pxAnchorMatch)
	{
		std::shared_ptr<StripOfConsideration> pxNew(new StripOfConsideration(pxAnchorMatch, uiQueryLength, uiStripSize));
		//BOOST_LOG_TRIVIAL(info) << "			starting to collect " << uiThreadId;
		for (unsigned int uiC = pxNew->getBegin() / uiStripSize - 1; uiC <= pxNew->getEnd() / uiStripSize; uiC++)
		{
			pxNew->addElement(apxPerfectMatchBuckets[uiC]);
		}//for

		return pxNew;
	}//function

	bool someOtherAnchorAlike(std::list<std::shared_ptr<const PerfectMatch>> &lpxUsedAnchors, std::shared_ptr<const PerfectMatch> pxAnchor)
	{
		for (auto pxCurr : lpxUsedAnchors)
		{
			auto uiAnchorPos = pxAnchor->getPositionForBucketing(uiQueryLength);
			auto uiCurrPos = pxCurr->getPositionForBucketing(uiQueryLength);
			if (uiAnchorPos <= uiCurrPos + uiStripSize && uiAnchorPos + uiStripSize >= uiCurrPos)
				return true;
		}//for
		return false;
	}//function

public: 
	AnchorMatchList(size_t uiNumThreads, nucSeqIndex uiStripSize, nucSeqIndex uiQueryLength, nucSeqIndex uiRefLength
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
		, AlignmentQuality *pxQuality
#endif
		)
			:
		apxPerfectMatchBuckets( (uiQueryLength + uiRefLength)/uiStripSize ),
		lpxAnchorMatches(),
		uiStripSize(uiStripSize),
		uiQueryLength(uiQueryLength)
#if confGENEREATE_ALIGNMENT_QUALITY_OUTPUT
		, pxQuality(pxQuality)
#endif
	{
		for (unsigned int uiC = 0; uiC < apxPerfectMatchBuckets.size(); uiC++)
		{
			apxPerfectMatchBuckets[uiC].reset(new PerfectMatchBucket());
		}//for
	}//constructor

	void addMatch(std::shared_ptr<const PerfectMatch> pxNew)
	{
		if(pxNew == nullptr)
			return;
		apxPerfectMatchBuckets[pxNew->getPositionForBucketing(uiQueryLength) / uiStripSize]->addMatch(pxNew);
	}//function

	void addAnchorSegment(std::shared_ptr<PerfectMatch> pxNew)
	{
		lpxAnchorMatches.push_back(pxNew);
	}//function

	void findAnchors(std::vector<std::shared_ptr<StripOfConsideration>> &aRet)
	{
		for (auto pxAnchorMatch : lpxAnchorMatches)
		{
			if (apxPerfectMatchBuckets[pxAnchorMatch->getPositionForBucketing(uiQueryLength) / uiStripSize]->isUsed())
				continue;

			apxPerfectMatchBuckets[pxAnchorMatch->getPositionForBucketing(uiQueryLength) / uiStripSize]->setUsed();
			aRet.push_back(collectStripOfConsideration(pxAnchorMatch));
		}//for
	}//function

};//class
#endif

class Bucketing: public Module
{
private:

public:
	unsigned int uiNumThreads = 8;
	nucSeqIndex uiStripSize = 1000;
	unsigned int uiMaxHitsPerInterval = 1000;
	bool bSkipLongBWTIntervals = true;
	
private:
	nucSeqIndex getPositionForBucketing(nucSeqIndex uiQueryLength, const std::shared_ptr<Seed> pS) const 
	{ 
		return pS->start_ref() + (uiQueryLength - pS->start()); 
	}//function

	void addSeed(
			nucSeqIndex uiQueryLength, 
			const std::shared_ptr<Seed> pNew, 
			std::vector<SeedBucket>& raxSeedBuckets
		)
	{
		raxSeedBuckets[getPositionForBucketing(uiQueryLength, pNew) / uiStripSize].addSeed(pNew);
	}//function

	void forEachNonBridgingSeed(
			std::shared_ptr<SegmentTreeInterval> pxNode,
			bool bAnchorOnly,
			std::shared_ptr<FM_Index> pxFM_index,
			std::shared_ptr<FM_Index> pxRev_FM_Index,std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
			std::shared_ptr<NucleotideSequence> pxQuerySeq,
			std::function<void(std::shared_ptr<Seed>)> fDo
		);
	
	void saveSeeds(
			std::shared_ptr<SegmentTreeInterval> pxNode,
			std::shared_ptr<FM_Index> pxFM_index,
			std::shared_ptr<FM_Index> pxRev_FM_Index,
			std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pxRefSequence,
			std::shared_ptr<NucleotideSequence> pxQuerySeq,
			std::vector<SeedBucket>& raxSeedBuckets
		);

public:

	Bucketing(){}//constructor

	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> vpInput);

    std::vector<ContainerType> getInputType();

	std::vector<ContainerType> getOutputType();
};//class
#if 0
class LineSweepContainer: public Module
{
public:

	LineSweepContainer(){}//constructor

	std::shared_ptr<Container> execute(std::vector<std::shared_ptr<Container>> pInput);

    std::vector<ContainerType> getInputType();

    std::vector<ContainerType> getOutputType();
};//class
#endif
void exportGraphicalMethod();

#endif
