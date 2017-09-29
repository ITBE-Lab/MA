/** 
 * @file exception.h
 * @brief Defines various exception classes that are used throughout the aligner.
 * @author Arne Kutzner
 * @author Markus Schmidt
 * @details
 * The exceptions are exposed to python.
 */

#pragma once

#include <string>
#include <cmath>
#include <boost/python.hpp>

/** 
 * @brief An annotated exception class on the foundation of std::exception.
 * @details
 * can be overloaded to provide different tpyes of exeptions, 
 * where each exception object can have it's unique description.
 */
template <int N>
class Annotated_exception : public std::exception
{
private :
	/* automated memory deallocation !
	 */
	std::string text;

public :
	/**
	 * @brief takes the string that shall be printed in case the exception is thrown.
	 * @details
	 * prepends information about the exception type to the string.
	 */
	Annotated_exception( const char* info )  
	{
		text = std::string( (N == 0) ? "NCBI data receiver/"		:
							(N == 1) ? "Smith Waterman Aligner"		:
							(N == 2) ? "Fasta Reader"				:
							(N == 3) ? "Null Pointer"				:
							(N == 4) ? "NCBI XML"					:
							(N == 5) ? "Boost ASIO"					:
							(N == 6) ? "FTP Download"				:
							(N == 7) ? "Module In/Out"					
									 : "Unknown Source"
						  );
		((text += " (") += info) += ")";
	} // constructor

	~Annotated_exception() throw() 
	{} // destructor
	
	/**
	 * @brief Information about the exception.
	 * @returns instance specific information about the exception.
	 */
	virtual const char* what() const throw() 
	{ 
		return text.c_str(); 
	} // method
};

/**
 * @brief Class for memory manager exceptions.
 */
class DataReceiverException : public Annotated_exception<0> 
{
public :
	DataReceiverException(const char* info) : Annotated_exception( info )
	{}
};

/**
 * @brief The aligner throws exceptions of this class.
 */
class AlignerException : public Annotated_exception<1> 
{
public :
	AlignerException(const char* info) : Annotated_exception( info )
	{}
};

/**
 * @brief Exceptions for the fasta reader.
 */
class fasta_reader_exception : public Annotated_exception<2> 
{
public :
	fasta_reader_exception(const char* info) : Annotated_exception( info )
	{}
};

/**
 * @brief Exceptions for null pointer.
 */
class NullPointerException : public Annotated_exception<3> 
{
public :
	NullPointerException(const char* info) : Annotated_exception( info )
	{}
};


/**
 * @brief Exceptions related to the NCBI XML parser.
 */
class NXBI_XML_Exception : public Annotated_exception<4> 
{
public :
	NXBI_XML_Exception(const char* info) : Annotated_exception( info )
	{}
};

/**
 * @brief Exceptions related to boost ASIO.
 */
class BOOST_ASIO_Exception : public Annotated_exception<5>
{
public:
	BOOST_ASIO_Exception( const char* info ) : Annotated_exception( info )
	{}
};

/**
 * @brief Exceptions for Download.
 */
class Download_Exeption : public Annotated_exception<6>
{
public:
	Download_Exeption(const char* info) : Annotated_exception(info){}
};

/**
 * @brief Exceptions for I/O.
 */
class ModuleIO_Exception : public Annotated_exception<7>
{
public:
		ModuleIO_Exception(const char* info) : Annotated_exception(info){}
};


/**
 * @brief A little method for null-pointer exception testing using our exception class.
 */
template<typename T>
T notNull( T pointer )
{
	if ( pointer == NULL )
	{
		throw NullPointerException( "" );
	} // if
} // generic function



/**
 * @brief Boost-python export function.
 * @details
 * Exports all exceptions to make them available in python.
 */
void exportExceptions();