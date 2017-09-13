#pragma once

#include <boost/log/trivial.hpp>

#ifndef DEBUG_LEVEL
	#define DEBUG_LEVEL 1
#endif

#if DEBUG_LEVEL >= 1
#define DEBUG(x) x
#else //DEBUG_LEVEL
#define DEBUG(x)
#endif //DEBUG_LEVEL

#if DEBUG_LEVEL >= 2
#define DEBUG_2(x) x
#else //DEBUG_LEVEL
#define DEBUG_2(x)
#endif //DEBUG_LEVEL

#if DEBUG_LEVEL >= 3
#define DEBUG_3(x) x
#else //DEBUG_LEVEL
#define DEBUG_3(x)
#endif //DEBUG_LEVEL

/* Overloading of the operator "<<" for arrays, so that we get automatic formated output for arrays.
 * From: http://stackoverflow.com/questions/19152178/printing-an-stdarray
 */
 template <class T, std::size_t N>
 std::ostream& operator<< ( std::ostream& rxOstream, const std::array<T, N> &raArray )
 {
	 std::copy( raArray.cbegin(), raArray.cend(), std::ostream_iterator<T>( rxOstream, " " ) );
	 return rxOstream;
 } // template function

/* +++++++++++++++++++++++++++++++++++++++++++++++++
 * 1. Application of some function f to all elements of some tuple.
 * Iteration starts with index 0.
 * Works with the empty tuple as well!
 * Taken from http://pastebin.com/h0Je0453
 */
template <int I, int TYPE_SIZE, typename Tuple>
struct iterate_over_tuple_impl : public iterate_over_tuple_impl<I + 1, TYPE_SIZE, Tuple>
{
	typedef typename std::tuple_element<I, Tuple >::type tp;
	
	template <typename Function>
	void operator () ( Function&& f, Tuple& t )
	{
		/* Application of the function to the i-th element of the tuple.
		 * Improvement: Deliver the I and TYPE_SZIE as integral template parameter.
		 */
		f( std::get<I>( t ), I, TYPE_SIZE );
		
		iterate_over_tuple_impl<I + 1, TYPE_SIZE, Tuple>::operator () ( std::forward<Function>(f), t );
	} // operator
}; // struct

template <int I, typename Tuple>
struct iterate_over_tuple_impl<I, I, Tuple>
{
	template <typename Function>
	void operator () ( Function&& f, Tuple& t ) 
	{}
}; // struct

/* Fills the element of the given tuple t by repeatedly calling the function f.
 * Starts with tuple element 0 and wok up to tuple element sizeof...(TupleTypes).
 */
template <typename Functor, typename... TupleTypes>
void iterateOverTuple( Functor&& functor, std::tuple<TupleTypes...>& tuple )
{
	iterate_over_tuple_impl<0, sizeof...(TupleTypes), std::tuple<TupleTypes...> >()
	(	std::forward<Functor>( functor ),
		tuple
	); // function call
}; // struct

/* Like iterateOverTuple but with an additional currying of 2 arguments.
 * (So, iterateOverTupleCurry2 requires a functor that expects two arguments.)
 */
template <typename Functor, typename... TupleTypes>
void iterateOverTupleCurry2( Functor&& functor, std::tuple<TupleTypes...>& tuple )
{
	iterate_over_tuple_impl<0, sizeof...(TupleTypes), std::tuple<TupleTypes...> >()
	(	std::bind( std::forward<Functor>( functor ), std::placeholders::_1, std::placeholders::_2 ), 
		tuple
	); // function call
}; // struct

/* +++++++++++++++++++++++++++++++++++++++++++++++++
 * 2. Pretty printer for tuples (uses similar iteration strategy like above iteration).
 * Iteration starts with index 0.
 * Does not work with the empty tuple!
 * Taken from http://pastebin.com/h0Je0453
 */
template<typename Type, unsigned N, unsigned Last>
struct tuple_printer
{
	static void print( std::ostream& out, const Type& value )
	{
		out << std::get<N>( value ) << ", ";
		tuple_printer<Type, N + 1, Last>::print( out, value );
	} // method
}; // struct

template<typename Type, unsigned N>
struct tuple_printer<Type, N, N>
{
	static void print( std::ostream& out, const Type& value )
	{
		out << std::get<N>( value );
	} // method
}; // struct

template<typename... Types>
std::ostream& operator<<(std::ostream& out, const std::tuple<Types...>& value)
{
	out << "(";
	tuple_printer<std::tuple<Types...>, 0, sizeof...(Types)-1>::print( out, value );
	out << ")";
	
	return out;
} // function

/* +++++++++++++++++++++++++++++++++++++++++++++++++
 * 3. Store parameter pack without expanding it.
 * Taken from http://stackoverflow.com/questions/4691657/is-it-possible-to-store-a-template-parameter-pack-without-expanding-it
 * Example of application.
 * typedef variadic_typedef<int, float> myTypes;
 * typedef convert_in_tuple<myTypes>::type int_float_tuple;
 */
template <typename... Args>
struct variadic_typedef
{
	// this single type represents a collection of types,
	// as the template arguments it took to define it
};

template <typename... Args>
struct convert_in_tuple
{
	// base case, nothing special,
	// just use the arguments directly
	// however they need to be used
	typedef std::tuple<Args...> type;
};

template <typename... Args>
struct convert_in_tuple<variadic_typedef<Args...>>
{
	// expand the variadic_typedef back into
	// its arguments, via specialization
	// (doesn't rely on functionality to be provided
	// by the variadic_typedef struct itself, generic)
	typedef typename convert_in_tuple<Args...>::type type;
};

/* +++++++++++++++++++++++++++++++++++++++++++++++++
 * 4. Unpack a tuple to call a matching function.
 * FIX ME: Better naming scheme, more documentation.
 * http://stackoverflow.com/questions/7858817/unpacking-a-tuple-to-call-a-matching-function-pointer
 */
template<int ...>
struct seq {};

template<int N, int ...S>
struct gens : gens<N - 1, N - 1, S...> {};

template<int ...S>
struct gens<0, S...> {
	typedef seq<S...> type;
};

/* For a solution without std::function look over here:
 * http://stackoverflow.com/questions/9535680/functions-functors-as-template-parameters-can-they-be-stored
 */
template<typename ...Args>
struct TupleUnpackToParameterByReference // DEPRECATED - use the below form now.
{	/* CHECK ME: Shouldn't we use a reference here? Hmm.. There is some discussion about on stackoverflow.
	 */
	const std::function<void( const Args& ... )> &func;

	template<int ...S>
	void callFunc( const std::tuple<Args...> &params, seq<S...> )
	{
		func( std::get<S>( params ) ... );
	}

	void operator()( const std::tuple<Args...> &params )
	{
		callFunc( params, typename gens<sizeof...( Args )>::type( ) );
	}

	TupleUnpackToParameterByReference( const std::function<void( const Args& ... )> &func )
		: func( func )
	{} // constructor
}; // struct

/* Unpacks a tuple for function application.
 */
template<typename ...Args>
struct TupleUnpackAndCallFunctor
{	const std::tuple<Args...> &_params;

	template<typename Functor, int ...S>
	void callFunc( Functor&& f, seq<S...> )
	{
		f( std::get<S>( _params ) ... );
	} // method

	template<typename Functor>
	void operator()( Functor&& f )
	{
		callFunc( std::forward<Functor>(f), typename gens<sizeof...(Args)>::type() );
	} // operator ()

	TupleUnpackAndCallFunctor( const std::tuple<Args...> &params )
		: _params( params )
	{ } // constructor

	/* Constructor that involves directly the functor call.
	 */
	template<typename Functor>
	TupleUnpackAndCallFunctor( Functor&& f, const std::tuple<Args...> &params )
		: _params( params )
	{ 
		callFunc( std::forward<Functor>(f), typename gens<sizeof...(Args)>::type() ); 
	} // constructor
}; // struct

/* +++++++++++++++++++++++++++++++++++++++++++++++++
 * 5. Apply some function to all arguments one by one.
 */
template <typename Function, typename Arg>
void metaApplyFunctionToAllArgumentsIndicatingFinalArgument( Function&& f, Arg&& arg )
{
	f( arg, true );
}; // struct

template <typename Function, typename Arg, typename... Args >
void metaApplyFunctionToAllArgumentsIndicatingFinalArgument( Function&& f, Arg&& arg, Args&&... args )
{
	f( arg, false );

	metaApplyFunctionToAllArgumentsIndicatingFinalArgument( f, std::forward<Args>(args)... );
}; // struct

/* Repeatedly executes the functor and catches the exception.
 * This function is reasonable. If we repeatedly want to do something, that can fail due to external reasons.
 * REMARK: If we wish smaller code size and faster compialtion, we should replace the functor in the template with a std::function argument.
 */
template <typename ReturnedValueType, typename FunctorType>
ReturnedValueType metaTryNTimesOrThrow( FunctorType&& functor, unsigned int uiNumberOfTries )
{
	std::string sFullExceptionText; // We keep the text of all exceptions.
	
	for ( decltype(uiNumberOfTries) uiCounter = 0; uiCounter < uiNumberOfTries; uiCounter++ )
	{
		try
		{
			return functor();

			/* We could successfully execute the functor without experiencing any problem.
			 * So, we return.
			 */
			
		} // try
		catch ( std::exception &rcException )
		{
			((((sFullExceptionText += "Try ") += std::to_string( uiCounter )) += " failed due to reason: ") += rcException.what()) += "\n";
		} // catch (std::exception)
		catch ( ... )
		{
			(((sFullExceptionText += "Try ") += std::to_string( uiCounter )) += " failed due unknown reason: ") += "\n";
		} // catch (eclipse)

		BOOST_LOG_TRIVIAL( trace ) << "metaTryNTimesOrThrow starts try number " << uiCounter;
	} // for

	
	throw std::runtime_error( std::string( "Gave up after ") + std::to_string( uiNumberOfTries ) + " tries with history:\n" + sFullExceptionText );
} // method

/* Try catch encapsulation for given functor.
 * In the case of some problem we report via logging mechanism.
 */
template <typename Functor>
void metaTryAndCatch( Functor&& functor )
{	try 
	{
		functor();
	} // try
	catch ( std::exception &rxException )
	{
		BOOST_LOG_TRIVIAL( warning ) << "Execution of functor failed due to reason: " << rxException.what();
	} // catch 
	catch ( ... )
	{
		BOOST_LOG_TRIVIAL( warning ) << "Execution of functor failed due to unknown reason.";
	} // catch
} // meta function

/* The following code belongs to sqlite.h but is "parked" over here.
 */
#if 0
/* From http://stackoverflow.com/questions/687490/how-do-i-expand-a-tuple-into-variadic-template-functions-arguments
 */
/**
 * Static Function Tuple Argument Unpacking
 *
 * This recursive template unpacks the tuple parameters into
 * variadic template arguments until we reach the count of 0 where the function
 * is called with the correct parameters
 *
 * @tparam N Number of tuple arguments to unroll
 */
template <unsigned int N, typename... Types>
struct apply_func
{
	/* the type of argument N - 1
	 */
	typedef typename std::tuple_element<N - 1, std::tuple<Types...> >::type ArgType;

	template < typename T, typename TA, typename... Args >
	static void applyTuple ( T& pObj,
							 TA &rxQuery,
							 Args && ...args )
	{
		// apply_func<N - 1, Types...>::applyTuple( pObj, rxQuery, dispatcher<N - 1, ArgType> ().aColumnElement<ArgType> ( rxQuery ), args... );
		apply_func<N - 1, Types...>::applyTuple( pObj, rxQuery, ArgType(), args... );
	}
};

//-----------------------------------------------------------------------------

/**
 * Static Function Tuple Argument Unpacking End Point
 *
 * This recursive template unpacks the tuple parameters into
 * variadic template arguments until we reach the count of 0 where the function
 * is called with the correct parameters
 */
template <typename... Types>
struct apply_func<0, Types...>
{
	template < typename T, typename TA, typename... Args >
	static void applyTuple ( T& pObj, 
							 TA &rxQuery,
							 Args && ... args )
	{
		pObj.emplace_back( args... );
	}
};

//-----------------------------------------------------------------------------

/**
 * Static Function Call Forwarding Using Tuple Pack Parameters
 */
// Actual apply function
template < typename T, typename TA, typename... Args >
void applyTuple( T &pObj, TA &rxQuery )
{
   apply_func<sizeof...(Args), Args ...>::applyTuple<T, TA>( pObj, rxQuery );
}


/* Helper function
 */
template <int I, typename T>
struct dispatcher
{
	template<typename A>
	inline A aColumnElement( CppSQLite3Query &rxQuery )
	{}

	template<>
	inline int aColumnElement<int>( CppSQLite3Query &rxQuery )
	{
		return rxQuery.getIntField( I );
	} // method

	template<>
	inline long long aColumnElement<long long>( CppSQLite3Query &rxQuery )
	{
		return rxQuery.getInt64Field( I );
	} // method

	template<>
	inline GeneticSequence aColumnElement<GeneticSequence>( CppSQLite3Query &rxQuery )
	{
		return GeneticSequence(); // std::string( rxQuery.getStringField( I ) ); // GeneticSequence(); // 
	} // method

	T get( CppSQLite3Query &rxQuery ) 
	{		
		return aColumnElement<T>( rxQuery );
	} // operator
}; // struct



				// iterate_over_tuple( dispatcher::aColumnElement(), xQuery, xSingleRowAsTuple );
			
				/* C++11, we move the fresh tuple into the result table.
				 * The content of xSingleRowAsTuple is void after this operation.
				 * http://stackoverflow.com/questions/18847424/c-create-custom-function-dispatcher-from-variadic-template
				 */
				//  applyTuple< std::vector< std::tuple<Types ... > >, CppSQLite3Query, Types ... >( pResultTableRef->xInternalTable, xQuery );
#if 0	
				typedef std::tuple_element<1, std::tuple<Types...> >::type myType;
				typedef std::tuple_element<0, std::tuple<Types...> >::type myType0;
				
				
				pResultTableRef->xInternalTable.emplace_back( dispatcher<0, myType0>().aColumnElement<myType0>(xQuery), dispatcher<1, myType>().aColumnElement<myType>(xQuery) );

				if ( sizeof...(Types) > 2 )
				{
					// typedef std::tuple_element<2, std::tuple<Types...> >::type myType2;
				}

				// pResultTableRef->xInternalTable.emplace_back( std::move(GeneticSequence()), dispatcher<1, myType>().aColumnElement<myType>(xQuery) );
#endif
				// pResultTableRef->xInternalTable.emplace_back( xSingleRowAsTuple );
#endif