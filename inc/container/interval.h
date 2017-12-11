/** 
 * @file interval.h
 * @brief Implements a generic Interval class.
 * @author Markus Schmidt
 */
#pragma once

#include <functional>
#include <iostream>
#include "exception.h"

/**
 * @brief A generic multipurpose Interval.
 */
template<typename T>
class Interval
{
public:
	/// @brief Start position of interval.
	T iStart;
	/// @brief Size of interval.
	T iSize;
	/**
	 * @brief Creates a new Interval.
	 */
	Interval(T iStart, T iSize)
		:
		iStart(iStart),
		iSize(iSize)
	{}// constructor

	/**
	 * @brief Copys from another Interval.
	 */
	Interval(const Interval& c)
		:
		iStart(c.iStart),
		iSize(c.iSize)
	{}// copy constructor

	/**
	 * @brief Default empty constructor.
	 */
	Interval()
	{}// default constructor

	inline static Interval start_end(T start, T end)
	{
		Interval xRet(start, 0);
		xRet.end(end);
		return xRet;
	}//function

	/**
	 * @returns the end of the interval.
	 * @brief Interval end.
	 * @note end = start = size
	 */
	inline const T end() const
	{
		return iStart + iSize;
	}// function

	///@brief Wrapper for boost-python.
	inline const T end_boost1() const
	{
		return iStart + iSize;
	}// function

	/**
	 * @brief Allows chaning the end of the interval.
	 * @note end = start = size
	 */
	inline void end(const T iVal)
	{
		iSize = iVal - iStart;
	}// function
	
	///@brief Wrapper for boost-python.
	inline void end_boost2(const T iVal)
	{
		end(iVal);
	}// function

	/**
	 * @returns the start of the interval.
	 * @brief Interval start.
	 */
	inline const T start() const
	{
		return iStart;
	}// function

	///@brief Wrapper for boost-python.
	inline const T start_boost1() const
	{
		return iStart;
	}// function

	/**
	 * @brief allows chaning the beginning of the interval.
	 */
	inline void start(const T iVal)
	{
		T iEnd = end();
		iStart = iVal;
		end(iEnd);
	}// function

	///@brief Wrapper for boost-python.
	inline void start_boost2(const T iVal)
	{
		start(iVal);
	}// function

	/**
	 * @brief allows chaning the size of the interval.
	 */
	inline void size(const T iVal)
	{
		iSize = iVal;
	}// function

	///@brief wrapper for boost-python
	inline void size_boost2(const T iVal)
	{
		iSize = iVal;
	}// function

	/**
	 * @returns the center of the interval.
	 * @brief The center of the interval.
	 */
	inline const T center() const
	{
		return start() + size()/2;
	}// function

	/**
	 * @returns the size of the interval.
	 * @brief The size of the interval.
	 */
	inline const T size() const
	{
		return iSize;
	}// function

	///@brief wrapper for boost-python
	inline const T size_boost1() const
	{
		return iSize;
	}// function

	/**
	 * @brief allows chaning the start and size of the interval.
	 */
	inline void set(T iStart, T iSize)
	{
		start(iStart);
		size(iSize);
	}//function

	/*
	 * @returns the interval start and end for i = 0 and 1 respectively.
	 * @brief Interval start for 0 and end for 1.
	 * @details
	 * Throws a nullpointer exception for any other i.
	 */
	inline const T operator[](const std::size_t i) const
	{
		if(i == 0)
			return start();
		if(i == 1)
			return end();
		throw new NullPointerException("can only access index 0 and 1 of interval");
	}//operator

	/*
	 * @brief copys from another Interval.
	 */
	inline Interval& operator=(const Interval& rxOther)
	{
		iStart = rxOther.iStart;
		iSize = rxOther.iSize;
		return *this;
	}// operator

	/*
	 * @brief compares two Intervals.
	 * @returns true if start and size are equal, false otherwise.
	 */
	inline bool operator==(const Interval& rxOther)
	{
		return iStart == rxOther.iStart && iSize == rxOther.iSize;
	}// operator
};//class


template<typename T>
void exportInterval()
{
    //export the Seed class
    boost::python::class_<
            Interval<T>
        >(
            "Interval",
            boost::python::init<T, T>()
        )
        .def_readwrite("start", &Interval<T>::iStart)
        .def_readwrite("size", &Interval<T>::iSize)
    ;
}//function
