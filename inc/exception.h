#pragma once

#include <string>
#include <boost/python.hpp>

/* An annotated exception class on the foundation of std::exception.
 */
template <int N>
class annotated_exception : public std::exception
{
private :
	/* automated memory deallocation !
	 */
	std::string text;

public :
	annotated_exception( const char* info )  
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

	~annotated_exception() throw() 
	{} // destructor
	
	virtual const char* what() const throw() 
	{ 
		return text.c_str(); 
	} // method
};

/* Class for memory manager exceptions.
 */
class DataReceiverException : public annotated_exception<0> 
{
public :
	DataReceiverException(const char* info) : annotated_exception( info )
	{}
};

/* The aligner throws exceptions of this class.
 */
class AlignerException : public annotated_exception<1> 
{
public :
	AlignerException(const char* info) : annotated_exception( info )
	{}
};

/* Exceptions for the fasta reader
 */
class fasta_reader_exception : public annotated_exception<2> 
{
public :
	fasta_reader_exception(const char* info) : annotated_exception( info )
	{}
};

/* Exceptions for null pointer
 */
class NullPointerException : public annotated_exception<3> 
{
public :
	NullPointerException(const char* info) : annotated_exception( info )
	{}
};


/* Exceptions related to the NCBI XML parser
 */
class NXBI_XML_Exception : public annotated_exception<4> 
{
public :
	NXBI_XML_Exception(const char* info) : annotated_exception( info )
	{}
};

/* Exceptions related to boost ASIO
 */
class BOOST_ASIO_Exception : public annotated_exception<5>
{
public:
	BOOST_ASIO_Exception( const char* info ) : annotated_exception( info )
	{}
};

//markus
/* Exceptions for Download
	*/
class Download_Exeption : public annotated_exception<6>
{
public:
	Download_Exeption(const char* info) : annotated_exception(info){}
};

class ModuleIO_Exception : public annotated_exception<7>
{
public:
		ModuleIO_Exception(const char* info) : annotated_exception(info){}
};
//end markus


/* A little method for null-pointer exception testing using our exception class
 */
template<typename T>
T notNull( T pointer )
{
	if ( pointer == NULL )
	{
		throw NullPointerException( "" );
	} // if
} // generic function




void exportExceptions();