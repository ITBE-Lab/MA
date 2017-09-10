#pragma once

#include <functional>
#include <iostream>
#include "exception.h"


template<typename T>
class Interval
{
private :
	// start position of interval
	T iStart;
	// size of interval
	T iSize;
public:
	Interval(T iStart, T iSize)
		:
		iStart(iStart),
		iSize(iSize)
	{}// constructor

	Interval(const Interval& c)
		:
		iStart(c.iStart),
		iSize(c.iSize)
	{}// copy constructor

	inline const T start() const
	{
		return iStart;
	}// function

	inline const T start_boost1() const
	{
		return iStart;
	}// function

	inline void start(const T iVal)
	{
		iStart = iVal;
	}// function

	inline void start_boost2(const T iVal)
	{
		iStart = iVal;
	}// function

	inline const T end() const
	{
		return iStart + iSize;
	}// function

	inline const T end_boost1() const
	{
		return iStart + iSize;
	}// function

	inline void end(const T iVal)
	{
		iSize = iVal - iStart;
	}// function
	
	inline void end_boost2(const T iVal)
	{
		iSize = iVal - iStart;
	}// function

	inline void size(const T iVal)
	{
		iSize = iVal;
	}// function

	inline void size_boost2(const T iVal)
	{
		iSize = iVal;
	}// function

	inline const T size() const
	{
		return iSize;
	}// function

	inline const T size_boost1() const
	{
		return iSize;
	}// function

	inline void set(T iStart, T iSize)
	{
		start(iStart);
		size(iSize);
	}//function

	inline const T operator[](const std::size_t i) const
	{
		if(i == 0)
			return start();
		if(i == 1)
			return end();
		throw new NullPointerException("can only access index 0 and 1 of interval");
	}//operator

	inline Interval& operator=(const Interval& rxOther)
	{
		iStart = rxOther.iStart;
		iSize = rxOther.iSize;
		return *this;
	}// operator

	inline bool operator==(const Interval& rxOther)
	{
		return iStart == rxOther.iStart && iSize == rxOther.iSize;
	}// operator
};//class
