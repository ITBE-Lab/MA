#include "graphicalMethod.h"


void exportGraphicalMethod()
{
	//export the StripOfConsideration class
	boost::python::class_<
        StripOfConsideration, 
        boost::python::bases<Container>, 
        std::shared_ptr<StripOfConsideration>
    >(
			"StripOfConsideration",
			"Holds the matches close to a selected anchor match\n"
		)
		.def(
			"get_score",
			&StripOfConsideration::getValue,
			"arg1: self\n"
			"returns: value of the content within the strip\n"
			"\n"
			"Note that the actual value of the content gets calculated over multiple steps.\n"
			"Therefore the returned value is not always correct.\n"
			"The returned value is however always higher than the actual score.\n"
		);

	//register a pointer to StripOfConsideration as return value to boost python
    boost::python::register_ptr_to_python< std::shared_ptr<StripOfConsideration> >();

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<StripOfConsideration>, 
			std::shared_ptr<Container>
		>(); 

	//register return values of vectors of strips
	boost::python::class_<std::vector<std::shared_ptr<StripOfConsideration>>>("VecRetStrip")
		.def(boost::python::vector_indexing_suite<
				std::vector<std::shared_ptr<StripOfConsideration>>,
				/*
				*	true = noproxy this means that the content of the vector is already exposed by
				*	boost python. 
				*	if this is kept as false, StripOfConsideration would be exposed a second time.
				*	the two StripOfConsiderations would be different and not intercastable.
				*	=> keep this as true
				*/
				true
			>());

	//export the StripOfConsiderationVector class
	boost::python::class_<
			StripOfConsiderationVector, 
			boost::python::bases<Container>, 
			std::shared_ptr<StripOfConsiderationVector>
		>(
			"StripOfConsiderationVector",
			"	x: the vector holding the strips\n."
			"\n"
			"Contains multiple strips of consideration.\n"
		)
		.def_readwrite("x", &StripOfConsiderationVector::x);
	
	//register a pointer to StripOfConsideration as return value to boost python
	boost::python::register_ptr_to_python< std::shared_ptr<StripOfConsiderationVector> >();

	//tell boost python that pointers of these classes can be converted implicitly
	boost::python::implicitly_convertible< 
			std::shared_ptr<StripOfConsiderationVector>, 
			std::shared_ptr<Container>
		>(); 
}//function