/**
 * @brief Template based metaprogramming used in the context of the DB module.
 * @author Arne Kutzner, Markus Schmidt
 * @date Nov 2019
 */
#include <tuple>
#include <functional>

/* METAPROGRAMMING
 * Unpack a tuple to call a matching function.
 * FIX ME: Better naming scheme, more documentation.
 * http://stackoverflow.com/questions/7858817/unpacking-a-tuple-to-call-a-matching-function-pointer
 * FIXME: Use the construct from the standard
 */
template <int...> struct Seq
{};

template <int N, int... S> struct Gens : Gens<N - 1, N - 1, S...>
{};

template <int... S> struct Gens<0, S...>
{
    typedef Seq<S...> type;
};


/* For a solution without std::function look over here: // DEPRECATED - use the below form now.
 * http://stackoverflow.com/questions/9535680/functions-functors-as-template-parameters-can-they-be-stored
 * FIXME: Use the construct from the standard
 */
template <typename... Args> struct TupleUnpackToParameterByReference
{ /* CHECK ME: Shouldn't we use a reference here? Hmm.. There is some discussion about on stackoverflow.
   */
    const std::function<void( const Args&... )>& func;

    template <int... S> void callFunc( const std::tuple<Args...>& params, Seq<S...> )
    {
        func( std::get<S>( params )... );
    }

    void operator( )( const std::tuple<Args...>& params )
    {
        callFunc( params, typename Gens<sizeof...( Args )>::type( ) );
    }

    TupleUnpackToParameterByReference( const std::function<void( const Args&... )>& func ) : func( func )
    {} // constructor
}; // struct


/** @brief Unpacks a tuple for function application.
 */
template <typename... Args> struct TupleUnpackAndCallFunctor
{
    const std::tuple<Args...>& _params;

    template <typename Functor, int... S> void callFunc( Functor&& f, Seq<S...> )
    {
        f( std::get<S>( _params )... );
    } // method

    template <typename Functor> void operator( )( Functor&& f )
    {
        callFunc( std::forward<Functor>( f ), typename Gens<sizeof...( Args )>::type( ) );
    } // operator ()

    TupleUnpackAndCallFunctor( const std::tuple<Args...>& params ) : _params( params )
    {} // constructor

    /* Constructor that involves directly the functor call. */
    template <typename Functor>
    TupleUnpackAndCallFunctor( Functor&& f, const std::tuple<Args...>& params ) : _params( params )
    {
        callFunc( std::forward<Functor>( f ), typename Gens<sizeof...( Args )>::type( ) );
    } // constructor
}; // struct


/** @brief Application of some function f to all elements of some tuple.
 * Iteration starts with index 0.
 * Works with the empty tuple as well!
 * Taken from http://pastebin.com/h0Je0453
 */
template <int I, int TYPE_SIZE, typename Tuple>
struct iterate_over_tuple_impl : public iterate_over_tuple_impl<I + 1, TYPE_SIZE, Tuple>
{
    typedef typename std::tuple_element<I, Tuple>::type tp;

    template <typename Function> void operator( )( Function&& f, Tuple& t )
    {
        /* Application of the function to the i-th element of the tuple.
         * Improvement: Deliver the I and TYPE_SZIE as integral template parameter.
         */
        f( std::get<I>( t ), I, TYPE_SIZE );

        iterate_over_tuple_impl<I + 1, TYPE_SIZE, Tuple>::operator( )( std::forward<Function>( f ), t );
    } // operator
}; // struct


/* METAPROGRAMMING */
template <int I, typename Tuple> struct iterate_over_tuple_impl<I, I, Tuple>
{
    template <typename Function> void operator( )( Function&& f, Tuple& t )
    {}
}; // struct


/** @brief Fills the element of the given tuple t by repeatedly calling the function f.
 * Starts with tuple element 0 and wok up to tuple element sizeof...(TupleTypes).
 */
template <typename Functor, typename... TupleTypes>
void iterateOverTuple( Functor&& functor, std::tuple<TupleTypes...>& tuple )
{
    iterate_over_tuple_impl<0, sizeof...( TupleTypes ), std::tuple<TupleTypes...>>( )( std::forward<Functor>( functor ),
                                                                                       tuple ); // function call
}; // struct


/** @brief Like iterateOverTuple but with an additional currying of 2 arguments.
 * (So, iterateOverTupleCurry2 requires a functor that expects two arguments.)
 */
template <typename Functor, typename... TupleTypes>
void iterateOverTupleCurry2( Functor&& functor, std::tuple<TupleTypes...>& tuple )
{
    iterate_over_tuple_impl<0, sizeof...( TupleTypes ), std::tuple<TupleTypes...>>( )(
        std::bind( std::forward<Functor>( functor ), std::placeholders::_1, std::placeholders::_2 ),
        tuple ); // function call
}; // struct

/* Delivers the n-th type in a type-list
 */
template <class... Args> struct type_list
{
	template <std::size_t N> using typex = typename std::tuple_element<N, std::tuple<Args...>>::type;
};

// Generic version of next power 2
// template <typename Type> Type next_power2( Type value )
// {
//     --value;
//     for( size_t uiCount = 1; uiCount < sizeof( Type ) * CHAR_BIT; uiCount *= 2 )
//         value |= value >> uiCount;
//     return value + 1;
// } // method
// Type uiMaxAlloc = static_cast<size_t>( std::numeric_limits<Type>::max( ) );
