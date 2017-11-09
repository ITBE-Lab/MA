/**
 * @file rmq.h
 * @brief implements a n-dimensional range maximum query
 * @author Markus Schmidt
 */

#ifndef RMQ_H
#define RMQ_H

#include "seed.h"
#include <memory>
#include <vector>
#include <iostream>
#include <algorithm>
#include "exception.h"
#include "debug.h"

#define COST_INS_DEL 2//lambda in the paper
#define COST_POSS_MATCH 1
#define SCORE_MATCH 3

class Chain;

/**
 * @brief A Range Maximum Query implementation.
 * @details
 * A RMQ is a N dimensional binary search tree,
 * with canonical nodes as leafs in the N-th dimension.
 * @see Chaining algorithms for multiple genome comparison 
 * [Mohamed Ibrahim Abouelhoda, Enno Ohlebusch]
 */
template<typename coordinate>
class RMQ
{
public:
	class RMQData{
	public:
		coordinate x, y;
		std::shared_ptr<Chain> pChain;
		int64_t score;

		RMQData(coordinate x, coordinate y, std::shared_ptr<Chain> pChain, int64_t score)
				:
			x(x),
			y(y),
			pChain(pChain),
			score(score)
		{}//constructor

		const bool operator< (const RMQData& other) const
		{
			return score < other.score;
			/*if(pChain == nullptr)
				return true;
			if(other.pChain == nullptr)
				return false;
			return *pChain < *other.pChain;*/
		}//function
	};//class

private:
    class Node
    {
    public:
		//array of left and right ptrs
		std::vector<unsigned int> yArray;

		Node(unsigned int start, unsigned int end, std::vector<RMQData>& data)
				:
			yArray()
			//,xPriorityQueue()
		{
			for(unsigned int i = start; i <= end; i++)
				yArray.push_back(i);
			std::sort(
					yArray.begin(), yArray.end(),
					[&]
					(unsigned int a, unsigned int b)
					{
						return data[a].y < data[b].y;
					}//lambda
				);//sort call
		}//constructor

		Node()
				:
			yArray()
			//,xPriorityQueue()
		{}//constructor

		virtual std::shared_ptr<Node> lPtr()
		{
			throw NullPointerException("cannot access child of leaf");
		}//function

        virtual std::shared_ptr<Node> rPtr()
		{
			throw NullPointerException("cannot access child of leaf");
		}//function

		virtual RMQData& rmq(
				coordinate x1, coordinate y1, coordinate x2, 
				coordinate y2, std::vector<RMQData>& data
			)
		{
			assert(yArray.size() > 0);
			if(y2 < data[yArray.front()].y)
				return data[0];
			if(y1 > data[yArray.back()].y)
				return data[0];
			unsigned int start = 0, end = yArray.size()-1;
			DEBUG_3(
				for(unsigned int i : yArray)
					assert(data.size() > i);
				for(unsigned int i : yArray)
					std::cout << data[i].y << ", ";
				std::cout << std::endl;
			)
			while( start != end )
			{
				unsigned int center = (start + end) / 2;
				if(data[yArray[center]].y < y1)
					start = center;
				else
					end = center;
				DEBUG_3(
					std::cout << center << std::endl;
				)
			}//while
			assert(data[yArray[start]].y >= y1);
			assert(start == 0 || data[yArray[start - 1]].y < y1);
			if(data[yArray[start]].y > y2)
				return data[0];//here we assume that the first element is equal to a nullptr...
			unsigned int max = yArray[start];
			unsigned int i = start + 1;
			DEBUG_2(
				if(data[yArray[start]].pChain != nullptr)
					std::cout << "Seed (ref/query/size/score/x/y):" 
						<< data[yArray[start]].pChain->s.end_ref() << " "
						<< data[yArray[start]].pChain->s.end() << " "
						<< data[yArray[start]].pChain->s.size() << " "
						<< data[yArray[start]].score << " "
						<< data[yArray[start]].x << " "
						<< data[yArray[start]].y << std::endl;
				else
					std::cout << "Seed dummy" << std::endl; 
			)
			while(i < yArray.size() && data[yArray[i]].y <= y2)
			{
				if(data[max].score < data[yArray[i]].score)
					max = yArray[i];
				DEBUG_2(
					if(data[yArray[i]].pChain != nullptr)
						std::cout << "Seed (ref/query/size/score/x/y):" 
							<< data[yArray[i]].pChain->s.end_ref() << " "
							<< data[yArray[i]].pChain->s.end() << " "
							<< data[yArray[i]].pChain->s.size() << " "
							<< data[yArray[i]].score << " "
							<< data[yArray[i]].x << " "
							<< data[yArray[i]].y << std::endl;
					else
						std::cout << "Seed dummy" << std::endl; 
				)
				i++;
			}//while
			return data[max];
		}//function

	};//class

    class Branch : public Node
    {
	public:
		std::shared_ptr<Node> pL, pR;
		//std::vector<std::tuple<unsigned int, unsigned int>> yArrayPtrs;
		coordinate iSplitVal;
		
		Branch(
				coordinate iSplitVal, std::shared_ptr<Node> pL, 
				std::shared_ptr<Node> pR, std::vector<RMQData>& data
			)
				:
			Node(),
			pL(pL),
			pR(pR),
			//yArrayPtrs(),
			iSplitVal(iSplitVal)
		{
			//kind of like mergesort for the yArrays...
			unsigned int i = 0, j = 0;
			while(i < pL->yArray.size() && j < pR->yArray.size())
			{
				//yArrayPtrs.push_back(std::make_tuple(i, j));
				if( data[pL->yArray[i]].y < data[pR->yArray[j]].y )
					Node::yArray.push_back(pL->yArray[i++]);
				else
					Node::yArray.push_back(pR->yArray[j++]);
			}//while
			while(i < pL->yArray.size())
			{
				//yArrayPtrs.push_back(std::make_tuple(i, j));
				Node::yArray.push_back(pL->yArray[i++]);
			}//while
			while(j < pR->yArray.size())
			{
				//yArrayPtrs.push_back(std::make_tuple(i, j));
				Node::yArray.push_back(pR->yArray[j++]);
			}//while
			assert(pL != nullptr);
			assert(pR != nullptr);
		}//constructor

        std::shared_ptr<Node> lPtr()
        {
            return pL;
        }//function

        std::shared_ptr<Node> rPtr()
        {
            return pR;
		}//function

		RMQData& rmq(
				coordinate x1, coordinate y1, coordinate x2, 
				coordinate y2, std::vector<RMQData>& data
			)
		{
			assert(pL != nullptr);
			assert(pR != nullptr);
			if(x2 < iSplitVal)
				return pL->rmq(x1, y1, x2, y2, data);
			if(x1 > iSplitVal)
				return pR->rmq(x1, y1, x2, y2, data);
			RMQData& a = pL->rmq(x1, y1, x2, y2, data);
			RMQData& b = pR->rmq(x1, y1, x2, y2, data);
			if(a.score < b.score)
				return b;
			return a;
		}//function
	};//class

	
	std::shared_ptr<Node> createRecursively(
			unsigned int start, unsigned int end, std::vector<RMQData>& data
		)
	{
		if(data[start].x == data[end-1].x)
			return std::shared_ptr<Node>(new Node(start, end-1, data));
		
		unsigned int center = ( start + end ) / 2;
		return std::shared_ptr<Node>(new Branch(
				data[center].x,
				createRecursively(start, center, data),
				createRecursively(center, end, data),
				data
			));
	}//function

	std::vector<RMQData>& data;
	std::shared_ptr<Node> root;

public:

	RMQ(std::vector<RMQData>& rNewdata)
			:
		data(rNewdata)
	{
		std::sort(
			data.begin(), data.end(),
			[]
			(const RMQData& a, const RMQData& b)
			{
				return a.x < b.x;
			}//lambda
		);//function call
		root = createRecursively(0, data.size(), data);
	}//constructor

	RMQData& rmq(coordinate x1, coordinate y1, coordinate x2, coordinate y2)
	{
		if(x2 < x1)
		{
			coordinate swap = x1;
			x1 = x2;
			x2 = swap;
		}//if
		if(y2 < y1)
		{
			coordinate swap = y1;
			y1 = y2;
			y2 = swap;
		}//if
		DEBUG_2(
			std::cout << "rmq (" << x1 << "," << y1 << ") - (" << x2 << ", " 
				<< y2 << ")" << std::endl; 
			)
		return root->rmq(x1, y1, x2, y2, data);
	}//function
};//class
void exportTest();

#endif