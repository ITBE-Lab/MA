/** 
 * @file graphicalMethod.h
 * @brief Implements StripOfConsideration.
 * @author Markus Schmidt
 */

#ifndef GRAPHICAL_METHOD_H
#define GRAPHICAL_METHOD_H


#include "threadPool.h"
#include "intervalTree.h"

class SeedBucket
{
private:
	nucSeqIndex uiTotalScore;

	std::list<Seed> lxContent;

	std::mutex xMutex;

	//bool bUsed;

public:
	SeedBucket()
		:
		uiTotalScore(0),
		lxContent(),
		xMutex()
		//,bUsed(false)
	{}//constructor

	////disable copying of buckets
	SeedBucket(const SeedBucket&) = delete;

	void addSeed(const Seed xNew)
	{
		std::lock_guard<std::mutex> xGuard(xMutex);
		lxContent.push_back(xNew);
		uiTotalScore += xNew.size();
		//end of scope xGuard
	}//function

	nucSeqIndex getValue()
	{
		return uiTotalScore;
	}//function

	void forall(std::function<void(const Seed&)> fDo)
	{
		for (auto xSeed : lxContent)
		{
			fDo(xSeed);
		}//for
	}//function
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
	std::list<Seed> axSeeds;

public:

	StripOfConsideration()
			:
		Interval(0, 0),
		uiValueOfContent(0),
		axSeeds()
	{}//constructor

	StripOfConsideration(nucSeqIndex uiStart, nucSeqIndex uiSize)
		:
		Interval(uiStart, uiSize),
		uiValueOfContent(0),
		axSeeds()
	{}//constructor

	
	/*used to identify the strip of cinsideration datatype in the aligner pipeline*/
	ContainerType getType(){return ContainerType::stripOfConsideration;}
	
	void addElement(SeedBucket& rxBucket)
	{
		rxBucket.forall(
			[&](const Seed& rxSeed)
			{
				axSeeds.push_back(rxSeed);
				uiValueOfContent += rxSeed.getValue();
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

	std::list<Seed>& seeds()
	{
		return axSeeds;
	}//function

	void forAllSeeds(std::function<void(std::list<Seed>::iterator pSeed)> fDo)
	{
		for(std::list<Seed>::iterator pSeed = axSeeds.begin(); pSeed != axSeeds.end(); pSeed++)
			fDo(pSeed);
	}//function
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

void exportGraphicalMethod();

#endif