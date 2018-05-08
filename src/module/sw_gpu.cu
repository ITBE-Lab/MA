/***********************************************
* # Copyright 2018. Arne Kutzner
* # Contact: Arne Kutzner
* #          kutzner@hanyang.ac.kr
* #
* # GPL 2.0 applies.
* **********************************************/

// NVIDIA card architecture: https://www.anandtech.com/show/3809/nvidias-geforce-gtx-460-the-200-king/2
// https://stackoverflow.com/questions/6647915/cuda-texture-memory-space
// Throw assert code: https://www.softwariness.com/articles/assertions-in-cpp/

#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <iomanip> // std::setw( )
#include <vector>
#include <memory> // for smart pointer
#include <algorithm>
#include <cstring>

#include <util/system.h> // only for measuring runtime

#ifdef __INTELLISENSE__
	#include "cuda_runtime_api.h"
#endif

#include "module/sw_gpu.h"
#include "module/sw_gpu_config.h"
#include "module/sw_gpu_defines.h"
#include "module/sw_gpu_kernel.h"

#if ( USE_THREADPOOL == 1)
	#include "util/threadPool.h"
#endif

#if ( DO_TESTS == 1 ) // define in config.h
	#include "../SW_naive.h"
#endif

/* Random nucleotide sequence of length uiLen, represented as codes.
 */
std::vector<char> randomNucleoitdeCodesSeq( const size_t uiLen )
{
	static const char nucleotides[] = { 0, 1, 2, 3 };

	std::vector<char> vNucSeq( uiLen );
	for (size_t i = 0; i < uiLen; ++i)
	{
		vNucSeq[i] = nucleotides[rand() % (sizeof( nucleotides ) - 1)];
	} // for

	return vNucSeq;
} // function

int xScoreMatrix[NUM_OF_SYMBOLS][NUM_OF_SYMBOLS] = { {10, -3, -3, -3},
													 {-3, 10, -3, -3},
													 {-3, -3, 10, -3},
													 {-3, -3, -3, 10} };

/* See http://en.cppreference.com/w/cpp/types/aligned_storage
 */
template<class T, std::size_t N>
class static_vector
{
	// properly aligned uninitialized storage for N T's
	typename std::aligned_storage<sizeof(T), alignof(T)>::type data[N];
	std::size_t m_size = 0;

public:
	// Create an object in aligned storage
	template<typename ...Args> void emplace_back(Args&&... args) 
	{
		if( m_size >= N ) // possible error handling
			throw std::bad_alloc{};
		new(data+m_size) T(std::forward<Args>(args)...);
		++m_size;
	}

	// Access an object in aligned storage
	const T& operator[](std::size_t pos) const 
	{
		// note: needs std::launder as of C++17
		return *reinterpret_cast<const T*>(data+pos);
	}

	// Delete objects from aligned storage
	~static_vector() 
	{
		for(std::size_t pos = 0; pos < m_size; ++pos) {
			// note: needs std::launder as of C++17
			reinterpret_cast<T*>(data+pos)->~T();
		}
	}
}; // class


//// // Texture reference for the query profile.
//// texture<int4, cudaTextureType2D, cudaReadModeElementType> texRef;

/* Reduce function on device side implementing OR.
 * This version uses contiguous threads, but its interleaved
 * addressing results in many shared memory bank conflicts.
 * For efficiency reasons we use int for boolean values.
 */
__device__ void
reduceOr( int *sdata )
{	/* The outcome is sdata[0]
	 */
	const unsigned int &tid = threadIdx.x;

	__syncthreads();

	// do reduction in shared mem
	for( unsigned int s = 1; s < blockDim.x; s *= 2 )
	{
		int index = 2 * s * tid;

		if( index < blockDim.x )
		{
			sdata[index] = sdata[index] || sdata[index + s];
		} // if

		__syncthreads();
	} // for
} // kernel

/* Reduce function on device side implementing MAX.
 * This version uses contiguous threads, but its interleaved
 * addressing results in many shared memory bank conflicts.
 */
__device__ void
reduceMax( int *sdata )
{	/* The outcome is sdata[0]
	*/
	const unsigned int &tid = threadIdx.x;

	__syncthreads();

	// do reduction in shared mem
	for( unsigned int s = 1; s < blockDim.x; s *= 2 )
	{
		int index = 2 * s * tid;

		if( index < blockDim.x )
		{
			sdata[index] = max(sdata[index], sdata[index + s]);
		} // if

		__syncthreads();
	} // for
} // kernel



/* Non-branching max(x, 0)
 * Instead of branching we have a comparison and multiplication.
 */
__device__ int
limitToZero( int iValue )
{
	//// return (iValue > 0) * iValue;
	return max( iValue, 0 );
} // device function



/* Non-branching maximum.
 * Instead of branching we have a comparison and multiplication.
 * Observation: Slower than standard max
 */
__device__ int
artihmMax( int iFirst, int iSecond )
{
	return iFirst + (iSecond > iFirst) * (iSecond - iFirst);
} // device function

typedef unsigned int int4hw;

/* Per halfword maximum computation*/
__device__ unsigned int
iv2_max( unsigned int uiFirstPair, unsigned int uiSecondPair )
{	/* per-halfword signed comparison: a > b ? 0xffff : 0. */
	auto uiBitVector = __vcmpgts2( uiFirstPair, uiSecondPair );
	return (uiBitVector & uiFirstPair) | (~uiBitVector & uiSecondPair);
}


#if ( SCORES_IN_CONSTANT_MEMORY == 1 )
__constant__ static int4 i4_Scores[256];
#endif

#if ( DO_TESTS == 1 )
//// -- NaiveSW<REFERENCE_SIZE, QUERY_SIZE> pNaiveSW;
#endif


/* Informs about the capabilities of a GPU card */
class GPU_CardInformer {
public:
	const unsigned int uiDeviceId; // DeviceId of the card
	cudaDeviceProp devProp; // Nvidia properties of the card


	GPU_CardInformer( unsigned int uiDeviceId ) :
		uiDeviceId( uiDeviceId )
	{
		cudaGetDeviceProperties( &this->devProp, uiDeviceId );
	} // constructor


	/* Delivers the available free memory for the device */
	size_t freeMemory()
	{	
		cudaSetDevice( this->uiDeviceId );
		size_t uiFreeMem, uiTotalMem;
		cudaMemGetInfo( &uiFreeMem, &uiTotalMem );
		CUERR

			return uiFreeMem;
	} // method


	/* Prints Info with respect to the card. */
	void printDevInfo()
	{
		std::cout 
			<< "Major revision number:         " << this->devProp.major << "\n" // compute capabilities
			<< "Minor revision number:         " << this->devProp.minor << "\n" // compute capabilities
			<< "Name:                          " << this->devProp.name << "\n"
			<< "Total global memory:           " << this->devProp.totalGlobalMem << "\n"
			<< "Total shared memory per block: " << this->devProp.sharedMemPerBlock << "\n"
			<< "Total registers per block:     " << this->devProp.regsPerBlock << "\n"
			<< "Warp size:                     " << this->devProp.warpSize << "\n"
			<< "Maximum memory pitch:          " << this->devProp.memPitch << "\n"
			<< "Maximum threads per block:     " << this->devProp.maxThreadsPerBlock << "\n"
			<< "Maximum threads per Multiproc.:" << this->devProp.maxThreadsPerMultiProcessor << "\n";
		for( int i = 0; i < 3; ++i )
			std::cout << "Maximum dimension " << i << " of block:   " << this->devProp.maxThreadsDim[i] << "\n";
		for( int i = 0; i < 3; ++i )
			std::cout << "Maximum dimension " << i << " of grid:    " << this->devProp.maxGridSize[i] << "\n";
		std::cout
			<< "Clock rate:                    " << devProp.clockRate << "\n"
			<< "Total constant memory:         " << devProp.totalConstMem << "\n"
			<< "Texture alignment:             " << devProp.textureAlignment << "\n"
			<< "Concurrent copy and execution: " << (devProp.deviceOverlap ? "Yes" : "No") << "\n"
			<< "Number of multiprocessors:     " << devProp.multiProcessorCount << "\n"
			<< "Kernel execution timeout:      " << (devProp.kernelExecTimeoutEnabled ? "Yes" : "No") << "\n"
			<< "Number of Concurrent Kernels:  " << devProp.concurrentKernels << "\n"
			<< "Compute Mode                   " << devProp.computeMode << "\n";
	} // method
}; // GPU_CardInformer


template<typename ELEMENT_TYPE>
class AlignedHostVector {
public:
	ELEMENT_TYPE *pHost; // The anchor of the host vector
	size_t uiSize; // size of host array in elements
	size_t uiSizeInBytes; // size of host array in bytes

	/* Constructor */
	AlignedHostVector( size_t uiSize ) :
		uiSize( uiSize ),
		uiSizeInBytes( uiSize * sizeof(ELEMENT_TYPE) )
	{
		cudaMallocHost( (void**)&this->pHost, uiSizeInBytes );
		CUERR
	} // constructor

	/* Delete copy constructor */
	AlignedHostVector(const AlignedHostVector& that) = delete;

	/* Move constructor */
	AlignedHostVector( AlignedHostVector&& rOther ) noexcept 
		: pHost ( rOther.pHost ),
		uiSize( rOther.uiSize )
	{
		std::cout << "AlignedHostVector move constructor" << std::endl;
		pHost = NULL;
		uiSize = 0;
	} // Move Constructor

	  /* Destructor */
	~AlignedHostVector()
	{
		if( pHost != NULL )
			cudaFreeHost( this->pHost );
	} // destructor
}; // class


/* Wrapper for device vector.
 * Introduction to unified memory: 
 * http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#um-unified-memory-programming-hd
 */
template <typename T_ELEMENT_TYPE>
class DeviceVector {
public:
	const size_t uiCapacity; // size of vector on device in elements

	/* Cache for the device vector */
	std::shared_ptr<std::vector<T_ELEMENT_TYPE>> pvSharedBackup = nullptr;

private:
	const size_t uiCapacityInBytes; // capacity of the vector in bytes
	bool bInUnifiedMemory; // flag that tells about unified memory
	size_t uiBytesUsed; // number of bytes filled by the last update of the vector
	std::unique_ptr<std::vector<T_ELEMENT_TYPE>> pHostCopy; // buffer for an host copy of device vector

	void copyDeviceToHost( void *pDestAddr )
	{
		if( this->bInUnifiedMemory )
		{	/* We can directly copy ... */
			std::memcpy( pDestAddr, this->pvAnchor, this->uiCapacityInBytes );
		} // if
		else
		{	/* We have to copy via cudaMemcpyHostToDevice */
			metaMeasureAndLogDuration<false>
				(	"cudaMemcpy get device requires",	// text message
					[&] () // lambda by reference
				{
					cudaMemcpy( pDestAddr, this->pvAnchor, this->uiCapacityInBytes, cudaMemcpyDeviceToHost );
				} // lambda
			); // function call
			CUERR
		} // else
	} // method

	//// /* Indexes are expressed in T_ELEMENT_TYPE */
	//// void copyFromTo( void *pDestAddr, size_t uiStart, size_t uiEnd )
	//// {
	//// 	assert( this->bInUnifiedMemory == false );
	//// 	assert( (uiStart < uiCapacity) && (uiEnd <= uiCapacity) && (uiStart <= uiEnd );
	//// 
	//// 	cudaMemcpy( pDestAddr, this->pvAnchor , this->uiCapacityInBytes, cudaMemcpyDeviceToHost );
	//// 	CUERR
	//// } // method 

public:
	T_ELEMENT_TYPE* pvAnchor; // anchor of vector on device

	/* Constructor */
	DeviceVector( size_t uiSize, // size of requested vector 
				  bool bAskForUnifiedMemory ) : // use unified memory if available (only available on Linux)
		pvAnchor( NULL ), // anchor of vector on device
		uiCapacity( uiSize ), // size of vector on device in elements
		uiCapacityInBytes( uiSize * sizeof(T_ELEMENT_TYPE) ), // size of vector on device in bytes
		bInUnifiedMemory( false ), // set unified memory flag (later use bAskForUnifiedMemory and check the card properties)
		uiBytesUsed( 0 )
	{
		//// size_t uiFreeMemBefore, uiFreeMemAfter, uiTotalMem;
		//// cudaMemGetInfo( &uiFreeMemBefore, &uiTotalMem );
		
		if( bInUnifiedMemory )
		{	/* Allocate unified memory */
			cudaMallocManaged( reinterpret_cast<void **>(&pvAnchor), uiCapacityInBytes );
		} // if
		else
		{	/* Allocate standard device memory */
			cudaMalloc( reinterpret_cast<void **>(&pvAnchor), uiCapacityInBytes );
		} // else
		CUERR //  TO DO: Raise exception

		/* Continue here be creating this vector if required */
		//// pHostCopy = std::make_unique<std::vector<T_ELEMENT_TYPE>>( this->uiCapacity );

		//// cudaMemGetInfo( &uiFreeMemAfter, &uiTotalMem );
		//// std::cout << "RAW: " << uiCapacityInBytes << " Real: " << uiFreeMemBefore - uiFreeMemAfter 
		//// 		  << " DIFF: " << (long)(uiFreeMemBefore - uiFreeMemAfter) - (long)uiCapacityInBytes << std::endl;
	} // constructor

	/* Constructor
	 * Vector is initialized by host-vector.
	 */
	DeviceVector( const std::vector<T_ELEMENT_TYPE> &rvSequence )
		: DeviceVector( rvSequence.size(), true ) // call basic constructor
	{
		this->updateDeviceVector( rvSequence );
	} // constructor

	/* Updates the Device-Vector with the content of the argument vector.
	 * Please note: Normally this is quite expansive, because the vector is not pinned.
	 */
	void updateDeviceVector( const std::vector<T_ELEMENT_TYPE> &rvHostVector )
	{
		assert( rvHostVector.size() <= uiCapacity ); // reject oversized inputs
		
		this->uiBytesUsed = rvHostVector.size() * sizeof( T_ELEMENT_TYPE );
		auto &rvSequenceNoConst = const_cast<std::vector<T_ELEMENT_TYPE> &>(rvHostVector);
		if( this->bInUnifiedMemory )
		{	/* We can directly copy ... */
			std::memcpy( this->pvAnchor, &rvSequenceNoConst[0], this->uiBytesUsed );
		} // if
		else
		{	/* We have to copy via cudaMemcpyHostToDevice */
			metaMeasureAndLogDuration<false>
			(	"cudaMemcpy update device requires",	// text message
					[&] () // lambda by reference
				{
					cudaMemcpy( this->pvAnchor, &rvSequenceNoConst[0], this->uiBytesUsed, cudaMemcpyHostToDevice );
				} // lambda
			); // function call
			
		} // else
		CUERR
	} // constructor


	/* Clear the vector on device side */
	void clear()
	{
		cudaMemset( this->pvAnchor, 0, this->uiCapacityInBytes );
		CUERR
	} // method

	/* Get a shared pointer to the backup vector */
	inline std::shared_ptr<std::vector<T_ELEMENT_TYPE>> getSharedBackup()
	{
		assert( this->pvSharedBackup != nullptr );
		return this->pvSharedBackup;
	} // method


	/* Update the backup vector */
	void updateSharedBackup()
	{
		if( this->pvSharedBackup == nullptr )
			this->pvSharedBackup = std::make_shared<std::vector<T_ELEMENT_TYPE>>( this->uiCapacity );

		this->copyDeviceToHost( &(*(this->pvSharedBackup))[0] );
	} // method


	/* The Host-vector should be aligned for maximum performance */
	void getCopyIntoVector( std::vector<T_ELEMENT_TYPE> &rvHostVector )
	{
		rvHostVector.resize( this->uiCapacity );
		this->copyDeviceToHost( &rvHostVector[0] );
	} // method


	/* Get a copy of the device-vector on host-side as STL-vector.
	 * This is quite expensive, because the vector must not be pinned.
	 */
	std::vector<T_ELEMENT_TYPE> getCopyAsVector() 
	{	
		std::vector<T_ELEMENT_TYPE> vReturnedVector( this->uiCapacity ); // time expensive
		this->copyDeviceToHost( &vReturnedVector[0] );
		return vReturnedVector;
	} // method


	/* Get a copy of the device-vector on host-side as aligned Host Vector.
	 * Copy for these vectors is a bit faster than for STL-vectors, but allocation takes more time.
	 */
	AlignedHostVector<T_ELEMENT_TYPE> getCopyAsAlignedHostVector()
	{
		AlignedHostVector<T_ELEMENT_TYPE> vReturnedVector( this->uiCapacity ); // time expensive
		this->copyDeviceToHost( vReturnedVector.pHost );
		return vReturnedVector;
	} // method


	/* Get a copy of the device-vector on host-side as array */
	std::unique_ptr<T_ELEMENT_TYPE[]> getCopyAsArray()
	{
		std::unique_ptr<T_ELEMENT_TYPE[]> puArray( new T_ELEMENT_TYPE[this->uiCapacity] );
		auto pArrayOnHost = puArray.get(); // get the inner pointer
		this->copyDeviceToHost( pArrayOnHost );
		return puArray;
	} // method


	/* Copy the content of the complete device vector to a host array */
	void copyToArray( T_ELEMENT_TYPE pArrayOnHost[] )
	{
		this->copyDeviceToHost( pArrayOnHost );
	} // method


	/* Dumps vectors that are holding pairs for debugging purposes. */
	void dump()
	{
		auto vCopyAtHost = this->getCopyAsVector();
		for( auto const &pair : vCopyAtHost )
		{
			std::cout << "(" << pair.x << ", " << pair.y << ") ";
		} // for
		std::cout<< std::endl;
	} // method


	/* Destructor*/
	~DeviceVector()
	{
		if( this->pvAnchor )
		{
			cudaFree( this->pvAnchor );
			CUERR
		} // if
	} // destructor
}; // class

/* Segmented device vector on the foundation of pitched CUDA-memory.
 * Here we keep the complete vector on device side.
 */
template <typename T_ELEMENT_TYPE>
class DeviceVector2D {
private:
	uint8_t* pCudaAnchor;
	size_t uiHeight; // height of the 2D matrix
	size_t uiWidth; // width (number of columns) expressed in elements
	size_t uiWidthInBytes; // width expressed in bytes (non-pitched value)
	size_t uiPitch; // expressed in bytes (pitched width)

public:
	DeviceVector2D( size_t uiHeight, size_t uiWidth ) :
		pCudaAnchor( NULL ),
		uiHeight( uiHeight ),
		uiWidth( uiWidth ),
		uiWidthInBytes( uiWidth * sizeof(T_ELEMENT_TYPE) )
	{	/* Allocate pitched memory */
		//// std::cout << "uiHeight " << uiHeight << " uiWidth " << uiWidth << std::endl;
		cudaMallocPitch( (void**)&pCudaAnchor, &this->uiPitch, this->uiWidthInBytes, this->uiHeight );
		CUERR
	} // constructor


	/* Delivers device anchor for row with index uiIndex */
	T_ELEMENT_TYPE* operator[] (size_t uiIndex )
	{	
		assert( uiIndex < this->uiHeight );
		return reinterpret_cast<T_ELEMENT_TYPE *>(this->pCudaAnchor + (uiIndex * this->uiPitch));
	} //method


	/* Fill the vector array on the foundation of the argument */
	void fill( const std::vector<T_ELEMENT_TYPE> &rvVector )
	{	/* Check correct size of input vector */
		assert( this->uiWidth * this->uiHeight == rvVector.size() );

		cudaMemcpy2D( this->pCudaAnchor, // destination anchor
					  this->uiPitch, // width on device (is device pitch)
					  &rvVector[0], // source anchor
					  this->uiWidthInBytes, // pitch of source memory
					  this->uiWidthInBytes, // width of matrix transfer in bytes
					  this->uiHeight, // height
					  cudaMemcpyHostToDevice ); // kind of copy
	} // method


	/* Destructor. Releases all resources */
	~DeviceVector2D()
	{	/* Release allocated CUDA memory. */
		if( pCudaAnchor != NULL )
			cudaFree( pCudaAnchor );
		CUERR
	} // destructor
}; // class


/* Query Profile for GPU SW.
 */
template<typename SCORE_TP4, typename SCORE_TP, int STRIPE_WIDTH>
class QueryProfile
{
private:
	/* Correctly casted local copy of scoring matrix */
	SCORE_TP4 xLocalScoreMatrix[NUM_OF_SYMBOLS]; // is aligned to 16 because of TP_4
	
	/* Returns a capacity that is multiple of STRIPE_WIDTH */
	size_t adjustCapacityToStripeWidth( size_t uiSuggestedCapacity )
	{
		assert( uiSuggestedCapacity > 0 );
		auto uiAdjustedCapacity = (((uiSuggestedCapacity - 1) / STRIPE_WIDTH) + 1) * STRIPE_WIDTH;
		return uiAdjustedCapacity;
	} // method

public:
	/* Public attributes */
	const size_t uiCapacity; // capacity of query profile (must be a multiple of STRIPE_WIDTH)

private:
	DeviceVector2D<SCORE_TP4> xDeviceVector; // device vector that keeps the profile
	size_t uiSize; // size of actual profile
	size_t uiNumberOfStripes; // number of stripes in actual profile

public:
	/* Constructor.
	 * uiCapacity is expressed in number of symbols maximally in the query.
	 */
	QueryProfile( const size_t uiSuggestedCapacity ) :
		uiCapacity( adjustCapacityToStripeWidth( uiSuggestedCapacity ) ),
		xDeviceVector( uiCapacity / STRIPE_WIDTH, STRIPE_WIDTH ), // allocated the device vector
		uiSize( 0 ), // profile is initially empty
		uiNumberOfStripes( 0 ) // profile is initially empty
	
	{	/* uiCapacity of stripe must be a multiple of query size */
		assert( this->uiCapacity % STRIPE_WIDTH == 0 );

		/* Initialize local scoring matrix by applying appropriate type cast */
		for( size_t uiRow = 0; uiRow < NUM_OF_SYMBOLS; uiRow++ )
			for( size_t uiColumn = 0; uiColumn < NUM_OF_SYMBOLS; uiColumn++ )
			{
				(reinterpret_cast<SCORE_TP *>(this->xLocalScoreMatrix))[uiRow * NUM_OF_SYMBOLS + uiColumn] 
					= static_cast<SCORE_TP>(xScoreMatrix[uiRow][uiColumn]);
			} // for
	} // constructor

	
	/* Delivers device anchor for row with index uiIndex.
	 * Delivers a pair consisting of a device-vector and stripe-size.
	 */
	std::pair<SCORE_TP4*, size_t> operator[] ( size_t uiStripeIndex )
	{	
		assert( uiStripeIndex < this->uiNumberOfStripes );
		/* Size of requested stripe (of stripe uiIndex) */
		size_t uiRemainder = this->uiSize % STRIPE_WIDTH;
		size_t uiStripeSize = uiStripeIndex >= this->uiNumberOfStripes - 1
									? (uiRemainder == 0 ? STRIPE_WIDTH : uiRemainder) // is last stripe 
									: STRIPE_WIDTH; // is some inner stripe
		assert( uiStripeSize > 0 );
		return std::make_pair( xDeviceVector[uiStripeIndex], uiStripeSize );
	} //method


	inline auto getNumberOfStripes()
	{
		return this->uiNumberOfStripes;
	} // method


	/* Computes the scoring profile on the foundation of the given query.
	 */
	void set( const std::vector<char> &rvQuerySeq )
	{	
		assert( rvQuerySeq.size() <= this->uiCapacity );
		assert( rvQuerySeq.size() > 0 );

		this->uiSize = rvQuerySeq.size();
		this->uiNumberOfStripes = ( ((rvQuerySeq.size() - 1) / STRIPE_WIDTH) + 1 );
		//// std::cout << "rvQuerySeq.size(): " << rvQuerySeq.size() << std::endl;
		//// std::cout << "this->uiNumberOfStripes: " << this->uiNumberOfStripes << std::endl;
		
		/* Compute profile vector and write it to the device */
		std::vector<SCORE_TP4> vHostVector( this->uiCapacity ); // auxiliary host vector
		for( size_t uiIndex = 0; uiIndex < this->uiSize; uiIndex++ )
		{
			vHostVector[uiIndex] = this->xLocalScoreMatrix[rvQuerySeq[uiIndex]];
		} // for
		xDeviceVector.fill( vHostVector );
	} // method
}; // class


//// template<typename SCORE_TP4, typename SCORE_TP2, typename SCORE_TP, int STRIPE_WIDTH>
//// class QueryDescriptor : public QueryProfile<SCORE_TP4, SCORE_TP, STRIPE_WIDTH> {
//// public:
//// 	std::vector<SCORE_TP2> vHostBackup; // Backup of the HE-vector on host side
//// 
//// 
//// 	/* constructor */
//// 	QueryDescriptor( const size_t uiCapacity ) :
//// 		QueryProfile<SCORE_TP4, SCORE_TP, STRIPE_WIDTH>( uiCapacity ),
//// 		vHostBackup( uiCapacity )
//// 	{} // constructor
//// 
//// 
//// 	void storeHE_Vector( DeviceVector<SCORE_TP2> &rxDeviceVector,
//// 						 size_t uiChunkId )
//// 	{	
//// 		assert( uiChunkId < this->uiNumberOfSegments );
//// 		assert( uiChunkId < this->uiNumberOfSegments );
//// 		
//// 		/* You did forget the number of reference segments over here.
//// 		 * We must save the last vector merely, because the inner vector will be uninteresting.
//// 		 */
//// 		//// std::cout << rxDeviceVector.uiCapacity << "   " << this->uiCapacity << std::endl;
//// 		assert( rxDeviceVector.uiCapacity == this->uiCapacity );
//// 
//// 		rxDeviceVector.copyToArray( &vHostBackup[0] );
//// 	} // method
//// }; // class

/* - Can transpose sequence (reference) for efficient GPU-processing
 * - Can pack/unpack sequences
 */
template<typename ELEMENT_TYPE>
class SequenceTransformer {
public:
	/* Copy source to destination and reorganizes the elements according to uiSegmentSize */
	static void copyTransposedTo( const ELEMENT_TYPE aSource[], // copy source
								  ELEMENT_TYPE aDestination[], // copy destination
								  size_t uiChunkSize, // size of chunks
								  size_t uiSourceSize, // size of source
								  size_t uiSegmentSize, // size of segment within chunk
								  size_t uiSourceOffset ) // offset with source
	{	/* Check whether destination size is sound */
		assert( uiChunkSize % uiSegmentSize == 0 ); // is a must!
		assert( uiSourceSize <= uiChunkSize );

		size_t uiNumberOfSegments = uiChunkSize / uiSegmentSize;
		/* TO DO: Could be implemented more efficiently via two nested loops */
		for (size_t uiItr = 0; uiItr < uiChunkSize; uiItr++)
		{
			size_t uiSegmentPos = uiItr / uiNumberOfSegments;
			size_t uiSegmentId = uiItr % uiNumberOfSegments;

			size_t uiRelativeSrcPos = uiSegmentId * uiSegmentSize + uiSegmentPos;
			aDestination[uiItr] = uiRelativeSrcPos < uiSourceSize ? aSource[uiSourceOffset + (uiRelativeSrcPos)] : 0;
		} // for
	} // method


	/* Reorganizes the sequence for optimal memory bandwidth on GPU (coalescing):
	 * TO DO: Do this in place, because it can be done by a simple swapping.
	 * Behavior for segment size 2:
	 * c0, c1, c2, c3, c4, c5, c6, c7 becomes:
	 * c0, c2, c4, c6, c1, c3, c5, c7
	 */
	static std::vector<ELEMENT_TYPE> transposedSeq( const std::vector<ELEMENT_TYPE> &rvSeq, // input vector
												    const size_t uiSegmentSize ) // size of segment within chunk
	{	/* TO DO: Check rvSeq.size() % uiSegmentSize == 0 */
		std::vector<ELEMENT_TYPE> vSegmentedSeq( rvSeq.size() );

		copyTransposedTo( &rvSeq[0], // copy source
						  &vSegmentedSeq[0], // copy destination
						  vSegmentedSeq.size(), // size of chunks
						  vSegmentedSeq.size(),	// length of source
						  uiSegmentSize, // segment size used for reordering
						  0 ); // reference offset

		return vSegmentedSeq;
	} // method


	/* Reverts a transposed vector to its original.
	 * Returns a copy; does not work in place.
	 * TO DO: Could be done in-place by using a juggling approach.
	 */
	static std::vector<ELEMENT_TYPE> inverseTransposedSeq( const std::vector<ELEMENT_TYPE> &rvSeq, // input vector
														   const size_t uiSegmentSize ) // segments size
	{	/* This could be done in-place by using the juggling-principle.
		 * See merging-papers for more information.
		 */
		std::vector<ELEMENT_TYPE> vReturnedSeq( rvSeq.size() );
		size_t uiNumberOfSegments = rvSeq.size() / uiSegmentSize;
		size_t uiRow = 0;
		for( size_t uiSegmentPos = 0; uiSegmentPos < uiSegmentSize; uiSegmentPos++ )
		{	
			for( size_t uiSegmentId = 0; uiSegmentId < uiNumberOfSegments; uiSegmentId++ )
			{
				vReturnedSeq[(uiSegmentId * uiSegmentSize) + uiSegmentPos ] = rvSeq[uiRow++];
			} // for
		} // for

		return vReturnedSeq;
	} // method
}; // class


/* Reference on GPU side.
 * On GPU side references are stored in transposed form in order to use coalescing.
 */
class ChunkedTransposedReference {
private:
	std::vector<std::shared_ptr<std::vector<char>>> vChunks;

public:
	const size_t uiSeqSize; // overall size of the input sequence
	const size_t uiChunkSize; // size of a single chunk consumed by the kernel
	const size_t uiNumOfSegments; // should be equal to BLOCK_SIZE * GRID_SIZE and comes directly form hardware
	const size_t uiNumberOfChunks; // number of big chunks. 
	const size_t uiSegmentSize; // size of segments inside the kernel (used for element transposing)

	/* Constructor */
	ChunkedTransposedReference( std::vector<char> &rvSequence, // input sequence
								const size_t uiChunkSize, // size of the chunks fed to the kernel
								const size_t uiNumOfSegments ) : // should be equal to BLOCK_SIZE * GRID_SIZE and chosen according to hardware
	    /* Prepare the segmented data-structure all chars zero-initialized */
		vChunks(),
		uiSeqSize( rvSequence.size() ),
		uiChunkSize( uiChunkSize ),
		uiNumOfSegments( uiNumOfSegments ),
		//// uiSegmentSize( uiSegmentSize ),
#if ( ALLOW_ABITRARY_REF_SIZE != 1 )
		uiNumberOfPrimarySegments( uiFullSize / uiPrimarySegmentSize ),
#else
		/* Compute number of required chunks */
		//// uiNumberOfChunks( ((uiSeqSize - 1) / uiChunkSize) + 1 ), // old calculation
		uiNumberOfChunks( (uiSeqSize + (uiChunkSize - 1)) / uiChunkSize ),
#endif
		//// uiNumOfSegments( uiChunkSize / uiSegmentSize )
		uiSegmentSize( uiChunkSize / uiNumOfSegments )
		
	{
		assert( uiChunkSize % uiNumOfSegments == 0 ); // this is an absolute must !
		//// assert( uiChunkSize % uiSegmentSize == 0 ); // this is an absolute must !
#if ( ALLOW_ABITRARY_REF_SIZE != 1 )
		assert( uiFullSize % uiPrimarySegmentSize == 0 );
#endif	

		for( size_t uiChunkId = 0; uiChunkId < this->uiNumberOfChunks; uiChunkId++ )
		{	/* Allocate memory for a single chunk.
			 * All chunks are of equal size, there are no undersized chunks 
			 */
			vChunks.push_back( std::make_shared<std::vector<char>>( uiChunkSize ) );
			//// for( auto &cSym : *vChunks.back() )
			//// {
			//// 	cSym = 0;
			//// } // for
			//// std::cout << "COPY SEG:" << uiChunkId << " ADDR: " << size_t(&(*(vChunks[uiChunkId]))[0]) << std::endl;
			SequenceTransformer<char>::copyTransposedTo( &rvSequence[0], // source anchor
														 &(*(vChunks[uiChunkId]))[0], // destination anchor
														 uiChunkSize, // size of chunks
														 numOfSymUsedInChunk( uiChunkId ), // length of source
														 uiSegmentSize, // size of segment within chunk
														 uiChunkId * uiChunkSize ); // offset within source
			
			//// std::cout << "SIZE:" << (*(vChunks[uiChunkId])).size() << std::endl;
			//// for( auto cSymbol : *(vChunks[uiChunkId]) )
			//// {
			//// 	assert( cSymbol >= 0 && cSymbol < NUM_OF_SYMBOLS );
			//// } // for
		} // for
	} // constructor

	/* Delivers the size of the chunk uiChunkId */
	size_t numOfSymUsedInChunk( size_t uiChunkId ) const
	{
		assert( this->uiNumberOfChunks > 0 );
		assert( uiChunkId < this->uiNumberOfChunks );

		if( uiChunkId < (this->uiNumberOfChunks - 1) )
		{	// return size of inner segment
			return this->uiChunkSize;
		} // if
		else
		{	// return size of last segment 
			assert( uiChunkId == uiNumberOfChunks - 1 );
			return this->uiSeqSize % this->uiChunkSize == 0 ? this->uiChunkSize
															: this->uiSeqSize % this->uiChunkSize;
		} // else
	} // method

	std::shared_ptr<std::vector<char>> getChunk( size_t uiIndex ) const
	{
		return vChunks.at( uiIndex ); // we do a range check
	} // method
}; // class

//// /* Encapsulates a reference for SW-kernel purposes */
//// class ReferenceSequenceHolder : public SequenceTransformer<char> {
//// private:
//// 	/* Helper function of the constructor */
//// 	std::vector<char> construct( const std::vector<char> &rvSeq,
//// 								 const size_t uiSegmentSize,
//// 								 bool bPacked )
//// 	{	/* Create a copy of the input-vector */
//// 		std::vector<char> vReturnedSeq( rvSeq );
//// 		
//// 		if( bPacked )
//// 		{	/* Do the packing of the vector */
//// 			pack2Vertically( vReturnedSeq );
//// 		} // of
//// 
//// 		/* Return transposed sequence */
//// 		return transposedSeq( vReturnedSeq, uiSegmentSize );
//// 	} // method
//// 
//// public:
//// 	const std::vector<char> vSequence; // keeps the transposed input sequence
//// 	const bool bPacked; // flag that indicates whether the vector is packed or not
//// 
//// 	/* Constructor */
//// 	ReferenceSequenceHolder( const std::vector<char> &rvSeq, // input sequence
//// 							 const size_t uiSegmentSize, // segment size for coalescing
//// 							 bool bPacked ) : // create packed vector
//// 		vSequence( construct( rvSeq, uiSegmentSize, bPacked ) ),
//// 		bPacked( bPacked )
//// 	{} // constructor
//// }; // class

//// /* Scoring Profile for a single stripe.
////  * The scoring profile should be updated before the stripe kernel is called.
////  */
//// template<typename SCORE_TP4, typename SCORE_TP, int STRIPE_WIDTH>
//// class StripeScoringProfile {
//// public:
//// 	/* (A -> int, C -> int, G -> int, T -> int */
//// 	DeviceVector<SCORE_TP4> xScoresOnDevice; 
//// 	std::vector<SCORE_TP4> vScoresOnHost;
//// 	/* Important observation: Never cast towards vector types like int4, because they are always aligned */
//// 	SCORE_TP4 xLocalScoreMatrix[NUM_OF_SYMBOLS]; // is aligned to 16 !
//// 	
//// 	/* Constructor */
//// 	StripeScoringProfile() :
//// 		xScoresOnDevice( STRIPE_WIDTH, true ), // initialize device vector
//// 		vScoresOnHost( STRIPE_WIDTH ) // initialize host vector
//// 	{	/* Initialize local scoring matrix by applying appropriate type cast */
//// 		for( size_t uiItr = 0; uiItr < NUM_OF_SYMBOLS; uiItr++ )
//// 			for( size_t uiColumn = 0; uiColumn < NUM_OF_SYMBOLS; uiColumn++ )
//// 			{
//// 				(reinterpret_cast<SCORE_TP *>(this->xLocalScoreMatrix))[uiItr * NUM_OF_SYMBOLS + uiColumn] 
//// 					= static_cast<SCORE_TP>(xScoreMatrix[uiItr][uiColumn]);
//// 			} // for
//// 	} // constructor
//// 
//// 	/* Updates the profile on device with the section of the query 
//// 	 * that starts at uiStartIndex 
//// 	 */
//// 	void updateDevice( const std::vector<char> &rvQuerySeq,
//// 					   const size_t uiStartIndex )
//// 	{	/* prepare the profile on host side */
//// 		for( size_t uiColIndex = 0; uiColIndex < STRIPE_WIDTH; uiColIndex++ )
//// 		{
//// 			vScoresOnHost[uiColIndex] = this->xLocalScoreMatrix[rvQuerySeq[uiStartIndex + uiColIndex]];
//// 		} // for
//// 		
//// 		/* Bring the profile to the device.
//// 		 * Update the device vector with the content of the host vector
//// 		 */
//// 		xScoresOnDevice.updateDeviceVector( vScoresOnHost );
//// 	} // method
//// }; // class


/* Works as manager for a collection of device vectors.
 * The class is not managing the device vectors.
 */
template<typename SCORE_TP4, // (int4 for 32 bit, int4 for 2x16 bit
		 typename SCORE_TP2, // (int2 for 32 bit, uint2 for 2x16 bit)
		 typename SCORE_TP, // (int for 32 bit, unsigned int for 2x16 bit)
		 typename CHECKSUM_TP>
class SW_GPU_MemoryCalculator : GPU_CardInformer {
public :
	const size_t uiCostPerRefRow;

	SW_GPU_MemoryCalculator( unsigned int uiDeviceId ) :
		GPU_CardInformer( uiDeviceId ),
		uiCostPerRefRow(   sizeof( SCORE_TP2 ) // HF_VectorOne
						 + sizeof( SCORE_TP2 ) // HF_VectorTwo
						 + sizeof( SCORE_TP ) // M_Vector
						 + sizeof( CHECKSUM_TP ) // C_Vector
						 + sizeof( char ) ) // Reference itself
		
	{} // constructor


	/* Fix-cost independent of any per row-cost.
	 * The alignment-cost per vector can be up to 
	 * Older cards: 1048576
	 * Newer cards: 2097152 
	 */
	size_t fixCost( size_t uiStripeWidth,
					size_t uiNumberOfSegments )
	{
		return   uiStripeWidth * sizeof( SCORE_TP4 ) // Scoring Profile (Query Profile)
				+ uiNumberOfSegments * sizeof( char ) // LazyFixedVector
				+ 2 * uiNumberOfSegments * uiStripeWidth * sizeof( SCORE_TP2 ) // HE_Cache
				+ uiStripeWidth * sizeof( SCORE_TP2 ); // HE_CarryOverVector
	} // method


	/* Computes the maximal number of rows for SW-computation with the current free memory */
	size_t maximalNumberOfRows( size_t uiStripeWidth,
								size_t uiNumberOfSegments )
	{
		auto uiFreeMem = this->freeMemory();
		auto uiFixCost = this->fixCost( uiStripeWidth, uiNumberOfSegments );
		/* 9 vectors are allocated */
		uiFixCost += 2097152 * 9; // estimated maximum costs of all alignments. Older card: 1048576 * 9

		if( uiFreeMem < uiFixCost )
			return 0;

		size_t uiMaxRows = (uiFreeMem - uiFixCost) / uiCostPerRefRow;
		
		/* The number of rows must be a multiple of the uiNumberOfSegments */
		return (uiMaxRows / uiNumberOfSegments) * uiNumberOfSegments;
	} // method
}; // class


 /* Host wrapper for a vector on device.
  * Important in the context of this data-structure is the alignment in the memory.
  * TO DO: Check, that vector sizes on the GPU do not reach beyond the 32-bit world.
  */
template<typename SCORE_TP4, // (int4 for 32 bit, int4 for 2x16 bit
		 typename SCORE_TP2, // (int2 for 32 bit, uint2 for 2x16 bit)
		 typename SCORE_TP, // (int for 32 bit, unsigned int for 2x16 bit)
		 typename CHECKSUM_TP, // 
		 unsigned int STRIP_WIDTH> // width of the strips used for GPU kernel calls
class SW_GPU_Processor : SW_GPU_MemoryCalculator<SCORE_TP4, SCORE_TP2, SCORE_TP, CHECKSUM_TP>
{
public:
	/* This value should be assigned depending on the GPU-hardware.
	 * Typically it is BLOCK_SIZE * GRID_SIZE.
	 * Represents the number of parallel threads.
	 */
	const size_t uiNumberOfSegments;

	/* The reference size is fix after construction. 
	 * The number of segments plus the reference size decide the structure of the outlay reference. 
	 * Warning the value should not be larger than 2^32, or we get problems with the kernel.
	 */
	const size_t uiRefCapicity;
	
	size_t uiSegmentSize; // size of individual segment in device

	/* Device pointer anchors.
	 * Possible improvement: Pitched HE-vector. (The reading of some columns might happen unaligned now.)
	 */
	DeviceVector<SCORE_TP2> xHF_VectorOne; // Size: reference size
	DeviceVector<SCORE_TP2> xHF_VectorTwo; // Size: reference size
	DeviceVector<SCORE_TP> xM_Vector; // Size: reference size
	DeviceVector<CHECKSUM_TP> xC_Vector; // Size: reference size
	DeviceVector<SCORE_TP2> xHE_CarryOverVector; // size: STRIPE_WIDTH
	DeviceVector<SCORE_TP2> xHE_CacheOne; // size: uiNumberOfSegments * STRIPE_WIDTH
	DeviceVector<SCORE_TP2> xHE_CacheTwo; // size: uiNumberOfSegments * STRIPE_WIDTH

	/* Reference vector on device. This vector has size uiRefCapicity.
	 */
	DeviceVector<char> xRefSeqTransposed; 
	
	/* Vector for inter segment communication with respect to the lazy E-loop. 
	 */
	DeviceVector<char> pLazyFixedVector; // size: uiNumberOfSegments
	// Host-copy of the LazyFixedVector
	std::unique_ptr<char[]> pContinuationFlags; // size: uiNumberOfSegments

#if ( USE_THREADPOOL == 1)
	ThreadPool xThreadPool;
	std::future<bool> xMaxExtractSuccess;
#endif

#if ( DO_TESTS == 1 )
	std::vector<char> vRefSeqDebug; // Pointer to reference sequence for debugging
#endif

	/* Clears all device vectors.
	 */
	void clearDeviceVectors()
	{	/* Set all vectors to fully 0 */
		xHF_VectorOne.clear();
		xHF_VectorTwo.clear();
		xHE_CacheOne.clear();
		xHE_CacheTwo.clear();
		xM_Vector.clear();
		xC_Vector.clear();

		pLazyFixedVector.clear();
	} // method

	/* Constructor.
	 * The SW-GPU-Processor gets a reference size as input and tries 
	 * to allocate the required amount of device memory.
	 */
	SW_GPU_Processor( unsigned int uiDeviceId, // device id of GPU
					  size_t uiRequestedSize ) : // requested length of reference sequence
		/* Get the memory calculator initialized */
		SW_GPU_MemoryCalculator<SCORE_TP4, SCORE_TP2, SCORE_TP, CHECKSUM_TP>( uiDeviceId ),

		/* Set the core capacities of the SW-GPU processor */
		uiNumberOfSegments( BLOCK_SIZE * GRID_SIZE ),
		uiRefCapicity( uiRequestedSize > 0 ? uiRequestedSize // 0 indicates user wishes specific size.
									   : this->maximalNumberOfRows( STRIP_WIDTH, uiNumberOfSegments ) ),
		
		/* Segments size is initially 0. This indicates that no reference has been loaded so far */
		uiSegmentSize( 0 ),

		/* Initialize the vectors and scoring profile */
		xHF_VectorOne( uiRefCapicity, true ),
		xHF_VectorTwo( uiRefCapicity, true ),
		xM_Vector( uiRefCapicity, true ),
		xC_Vector( uiRefCapicity, false ),

		xHE_CarryOverVector( STRIP_WIDTH, false ),
		xHE_CacheOne( uiNumberOfSegments * STRIP_WIDTH, true ),
		xHE_CacheTwo( uiNumberOfSegments * STRIP_WIDTH, true ),
		
		xRefSeqTransposed( uiRefCapicity, true ),
		
		/* Initialize the continuation vectors */
		pLazyFixedVector( uiNumberOfSegments, true ),
		pContinuationFlags( new char[this->uiNumberOfSegments] ) // unique pointer
#if ( DO_TESTS == 1 )
		//// , vRefSeqDebug( NULL )
#endif
#if ( USE_THREADPOOL == 1)
		, xThreadPool( 1 )
#endif

	{	/* Reference size must be a multiple of number of segments */
		assert( uiRefCapicity % uiNumberOfSegments == 0 );

		/* Demands the logic of Smith-Waterman's algorithm */
		assert( BLOCK_SIZE >= STRIP_WIDTH );
		
		//// /* The lazy-fixed loop cast from char to long for optimization purposes.
		////  * This cast would be no longer OK if the segment size is no multiple of sizeof(long).
		////  */
		//// assert( uiNumberOfSegments % sizeof( long ) == 0 );

		//// std::cout << "Reference capacity: " << this->uiRefCapicity << std::endl;		
		clearDeviceVectors();
	} // constructor

	/* Load a reference to the SW-GPU processor.
	 * Here we have to deliver a transposed sequence. 
	 */
	void loadReference( const std::vector<char> &rvRefSeq )
	{
		assert( this->uiRefCapicity == rvRefSeq.size() );
		assert( rvRefSeq.size() % this->uiNumberOfSegments == 0 );
		
#if ( DO_TESTS == 1 )
		/* In debug-mode keep a pointer to the loaded reference */
		//// for( auto cSymbol : rvRefSeq )
		//// {
		//// 	assert( cSymbol >= 0 && cSymbol < NUM_OF_SYMBOLS );
		//// } // for

		vRefSeqDebug = SequenceTransformer<char>::inverseTransposedSeq( rvRefSeq, this->uiSegmentSize );
		//// for( auto cSymbol : vRefSeqDebug )
		//// {
		//// 	assert( cSymbol >= 0 && cSymbol < NUM_OF_SYMBOLS );
		//// } // for
		//// vRefSeqDebug = &rvRefSeq;
#endif
		
#if ( USE_PACKED == 1 )
		/* Vector must be packed before transposing */
		pack2Vertically( vReturnedSeq );
#endif

		xRefSeqTransposed.updateDeviceVector( rvRefSeq );
	} // method


	/* Backup function for the "last row" of the active HE-Cache.
	 * Possible Improvement: In the case of unified memory the vector-copy can be avoided.
	 */
	void storeHE_Backup( const size_t uiStripeId, // id of the stripe which receives a backup
						 std::vector<SCORE_TP2> &rvBackupVector, // vector that receives the backup; has size of query
						 DeviceVector<SCORE_TP2>* rxHE_Cache ) // current active HE-cache
	{
		assert( (uiStripeId + 1) * STRIP_WIDTH <= rvBackupVector.size() );
		assert( rxHE_Cache->uiCapacity == STRIP_WIDTH * this->uiNumberOfSegments ); 

		/* Get a host-backup of the device vector */
		std::vector<SCORE_TP2> vHE_CacheCopy = rxHE_Cache->getCopyAsVector();

		size_t uiIndex = this->uiNumberOfSegments - 1; // HE-vector-pitch - 1
		size_t uiOffset = uiStripeId * STRIP_WIDTH;
		for( size_t uiCounter = 0; uiCounter < STRIP_WIDTH; uiCounter++ )
		{
			//// assert( uiOffset + uiCounter < rvBackupVector.size() );
			rvBackupVector[uiOffset + uiCounter] = vHE_CacheCopy[uiIndex];
			uiIndex += uiNumberOfSegments; // HE-vector-pitch
		} // for

		//// /* Get the H-Value of the last row of the HF-vector */
		//// std::vector<SCORE_TP2> vHF_VectorCopy = rxHF_Vector->getCopyAsVector();
		//// rH_ValueBackup = vHF_VectorCopy[vHF_VectorCopy.size() - 1].x;
	} // method


	/* Restores the first row of the HE-vector on the foundation of the backup-vector.
	 * Possible Improvement: In the case of unified memory the vector-copy twice can be avoided.
	 */
	void restoreHE_Backup( const size_t uiStripeId, // id of the stripe which receives a backup
						   const std::vector<SCORE_TP2> &rvBackupVector, // vector that receives the backup; has size of query
						   DeviceVector<SCORE_TP2>* rxHE_CacheSegment0 ) // current active HE-vector
	{
		assert( (uiStripeId + 1) * STRIP_WIDTH <= rvBackupVector.size() );
		assert( rxHE_CacheSegment0->uiCapacity == STRIP_WIDTH );
#if ( 1 )
		/* Do appropriate sub-vector copy */
		auto pFirst = rvBackupVector.begin() + (uiStripeId * STRIP_WIDTH);
		auto pLast = rvBackupVector.begin() + ((uiStripeId + 1) * STRIP_WIDTH);
		std::vector<SCORE_TP2> vStripePart( pFirst, pLast );
		this->xHE_CarryOverVector.updateDeviceVector( vStripePart );

		//// /* Extract the Left-Up-H-Value */
		//// if( uiStripeId == 0 )
		//// {	/* Here we do not have any preceding stripe ... */
		//// 	rLeftUpH_Value = 0;
		//// } // if
		//// else
		//// {
		//// 	auto pFirst = rvBackupVector.begin() + ((uiStripeId - 1) * TILE_WIDTH);
		//// 	auto pLast = rvBackupVector.begin() + (uiStripeId * TILE_WIDTH);
		//// 	std::vector<SCORE_TP2> vStripePart( pFirst, pLast );
		//// 	rLeftUpH_Value = vStripePart[vStripePart.size() - 1].x;
		//// } // else

#else // different HE-vector implementation
		/* Get a host-backup of the device vector.
		 * In the case of unified memory we can avoid this copy.
		 */
		auto vHE_CacheCopy = rxHE_CacheSegment0->getCopyAsVector();

		size_t uiIndex = 0; // We work now with the first row.
		size_t uiOffset = uiStripeId * TILE_WIDTH;
		for( size_t uiCounter = 0; uiCounter < TILE_WIDTH; uiCounter++ )
		{
			vHE_CacheCopy[uiIndex] = rvBackupVector[uiOffset + uiCounter];
			uiIndex += uiNumberOfSegments; // HE-vector-pitch
		} // for

		/* In the case of unified memory we can avoid this copy */
		rxHE_CacheSegment0->updateDeviceVector( vHE_CacheCopy );
#endif
	} // method


	/* The kernel is either in the normal mode or lazy-mode.
	 * The lazy-mode is similar to the lazy F-loop, but works with a checksum.
	 */
	template<bool LAZY_MODE>
	void callTileKernel( bool bHF_VectorsReversed,
						 bool bHE_CacheReversed,
						 SCORE_TP4* pQueryProfileOfStripe, // query profile for stripe
						 SCORE_TP iLeftUpH_Value,
						 const unsigned int uiWidth )  // if FULL_STRIPE_WIDTH is false we get the width of the stripe here 
	{
		dim3 xDimBlock( BLOCK_SIZE, 1 ); // (x, y)
		dim3 xDimGrid( GRID_SIZE, 1 ); // (x, y)
		assert( BLOCK_SIZE * GRID_SIZE * this->uiSegmentSize == this->uiRefCapicity );

		/* Set read and write vector according to the reversed flag */
		auto pvHF_VectorRead = bHF_VectorsReversed ? this->xHF_VectorTwo.pvAnchor : this->xHF_VectorOne.pvAnchor;
		auto pvHF_VectorWrite = bHF_VectorsReversed ? this->xHF_VectorOne.pvAnchor : this->xHF_VectorTwo.pvAnchor;

		auto pvHE_CacheRead = bHE_CacheReversed ? this->xHE_CacheTwo.pvAnchor : this->xHE_CacheOne.pvAnchor;
		auto pvHE_CacheWrite = bHE_CacheReversed ? this->xHE_CacheOne.pvAnchor : this->xHE_CacheTwo.pvAnchor;
		
		/* Call the tile kernel .
		 * In the true branch the compiler unrolls the column iteration in the kernel,
		 * which results in a strong performance boost.
		 */
		if( uiWidth == STRIP_WIDTH )
		{
			tileKernel<SCORE_TP4, SCORE_TP2, SCORE_TP, CHECKSUM_TP, BLOCK_SIZE, STRIP_WIDTH, LAZY_MODE, 
					  true> // call with FULL_STRIPE_WIDTH == true -> This results in unrolling the column iteration 
					  <<<xDimGrid, xDimBlock >>>
				( pvHF_VectorRead, // Vector has reference size
				  pvHF_VectorWrite, // Vector has reference size
				  this->xM_Vector.pvAnchor, // Vector has reference size
				  this->xC_Vector.pvAnchor, // Vector has reference size
				  this->xHE_CarryOverVector.pvAnchor,
				  iLeftUpH_Value,
				  pvHE_CacheRead, // Ingoing row-cache
				  pvHE_CacheWrite, // Outgoing row-cache
				  this->xRefSeqTransposed.pvAnchor, // Address of reference on device
				  pQueryProfileOfStripe, // Query-profile for stripe
				  this->pLazyFixedVector.pvAnchor, // Communicates the need of a continuation
				  (unsigned int)this->uiSegmentSize, // Size of each segment (TILE_HEIGHT)
				  STRIP_WIDTH ); // The argument is only active in the case of FULL_STRIPE_WIDTH is false
		} // if
		else
		{
			tileKernel<SCORE_TP4, SCORE_TP2, SCORE_TP, CHECKSUM_TP, BLOCK_SIZE, STRIP_WIDTH, LAZY_MODE, 
					  false> // call with FULL_STRIPE_WIDTH == false -> Without optimization 
					  <<<xDimGrid, xDimBlock >>>
				( pvHF_VectorRead, // Vector has reference size
				  pvHF_VectorWrite, // Vector has reference size
				  this->xM_Vector.pvAnchor, // Vector has reference size
				  this->xC_Vector.pvAnchor, // Vector has reference size
				  this->xHE_CarryOverVector.pvAnchor,
				  iLeftUpH_Value,
				  pvHE_CacheRead, // Ingoing row-cache
				  pvHE_CacheWrite, // Outgoing row-cache
				  this->xRefSeqTransposed.pvAnchor, // Address of reference on device
				  pQueryProfileOfStripe, // Query-profile for stripe
				  this->pLazyFixedVector.pvAnchor, // Communicates the need of a continuation
				  (unsigned int)this->uiSegmentSize, // Size of each segment (TILE_HEIGHT)
				  uiWidth ); // The argument is only active in the case of FULL_STRIPE_WIDTH is false
		} // else
		CUERR

		/* Change: Work with two maximum-vectors, One active vector and the one of the previous computation.
		 * Use the CPU-time that we wait for cudaDeviceSynchronize() for processing the previous maximum vector.
		 */

		/* Wait until the kernel finished its job */
		cudaDeviceSynchronize();
	} // method


	 /* Updates all vectors with a fresh symbol of the query.
	  * the resulting kernel call computes a single column of the matrix.
	  * WARNING: Without cudaDeviceSynchronize(); we get CUDA-errors.
	  */
	void doStrip( bool bHF_VectorsReversed,
				  const size_t uiStripeId, // id of current stripe
				  std::vector<SCORE_TP2> &rvBackupVector, // HE-backup vector
				  SCORE_TP &riLeftUpH_Value, // second component of backup, the H value-backup
				  SCORE_TP4* pProfileOnDevice, // query profile for stripe
				  const unsigned int uiWidth ) // number of inspected elements 
	{
		assert( uiWidth <= STRIP_WIDTH );
		bool bHE_CacheReversed = false;

		/* initialize incoming HE-vector for the first segment for the current stripe.
		 * Clear the primary HE-vector.
		 * Working with two vectors is a fix, because if we work with one vector we need a pitched data-structure.
		 */
		this->restoreHE_Backup( uiStripeId, rvBackupVector, &(this->xHE_CarryOverVector) ); // current active HE-vector

		this->xHE_CacheOne.clear();

		/* First call of the tile kernel happens in the "init mode" */
		this->callTileKernel<false>( bHF_VectorsReversed, bHE_CacheReversed, pProfileOnDevice, riLeftUpH_Value, uiWidth );

		/* TO DO: Use a counter. Just for paranoid programmers ...*/
		while( true )
		{	/* With each iteration the cache-vectors have to be swapped.
			 * In the beginning the cache vectors have to be cleared.
			 */
			bHE_CacheReversed = !bHE_CacheReversed;

			/* Call the stripe kernel using the lazy-mode, where the checksum decides about continuation */
			this->callTileKernel<true>( bHF_VectorsReversed, bHE_CacheReversed, pProfileOnDevice, riLeftUpH_Value, uiWidth );
	
			if( this->lazyFixed() )
			{
				break;
			} // if
		} // while
#if ( DO_TESTS == 1 )
		//// std::cout << "**** Check Last column:" << std::endl;
		//// this->debugCheckHF_VectorTwo( QUERY_SIZE - 1, bHF_VectorsReversed );
#endif
		/* Save the up-left H-Value for the next call of doStrip */
		riLeftUpH_Value = rvBackupVector[((uiStripeId + 1) * STRIP_WIDTH) - 1].x;
		
		/* Backup the final values of the HE-vector for the current stripe */
		storeHE_Backup( uiStripeId, // id of the stripe which receives a backup
						rvBackupVector, // vector that receives the backup; has size of query
						//// rH_ValueBackup, // H value backup
						bHE_CacheReversed ? &(this->xHE_CacheTwo) : &(this->xHE_CacheOne) ); // current active HE-vector
						//// bHF_VectorsReversed ? &(this->xHF_VectorTwo) : &(this->xHF_VectorOne) ); // current active HF-vector
	} // method


	/* Do the GPU-SW for the given query on the given reference segment.
	 * Changes: Instead of the query deliver the complete query profile.
	 *		    Deliver a reference id, that indicates what reference shall be processed.
	 */
	void doQueryForChunk( QueryProfile<SCORE_TP4, SCORE_TP, STRIP_WIDTH> &rxQueryProfile,
#if ( DO_TESTS == 1 )
							   std::vector<char> &rvQuerySeq, 
#endif
							   std::vector<SCORE_TP2> &rvHE_BackupVector)
							   //// SCORE_TP &rH_ValueBackup, // second component of backup, the H value-backup) // std::vector<char> &rvQuerySeq )
	{	
		/* According to the standard the vector should be zero initialized.
		 * See: https://stackoverflow.com/questions/25198405/are-members-of-structs-in-a-vector-initialized-with-zero-in-c
		 */
		//// ** std::vector<SCORE_TP2> vHE_BackupVector( rvQuerySeq.size() );

#if ( DO_TESTS == 1 )
		/* Fill the SW-matrix for later debugging */
		for( auto cSymbol : this->vRefSeqDebug )
		{
			assert( cSymbol >= 0 && cSymbol < NUM_OF_SYMBOLS );
		} // for
		//// std::shared_ptr< NaiveSW<CHUNK_SIZE, QUERY_SIZE> > pNaiveSW 
		//// 	= std::make_shared<NaiveSW<CHUNK_SIZE, QUERY_SIZE>>();
		//// pNaiveSW->fillMatrix( this->vRefSeqDebug, rvQuerySeq, xScoreMatrix );
#endif

		/* The core loop that iterates over the query in steps of size TILE_WIDTH */
		bool bHF_VectorsReversed = false;
		SCORE_TP iLeftUpH_Value = 0;
		
		for( size_t uiStripeId = 0; uiStripeId < rxQueryProfile.getNumberOfStripes() ; uiStripeId++ )
		{	/* Create the scoring profile in the local vector.
			 * HINT: This is quite expensive.
			 * The profile could be fully prepared by the CPU before calling the GPU.
			 */
			std::pair<SCORE_TP4*, size_t> xStripeData = rxQueryProfile[uiStripeId];
			unsigned int uiStripWidth = static_cast<unsigned int>(xStripeData.second);
			/* Compute matrix for the stripe*/
			this->doStrip( bHF_VectorsReversed,
							uiStripeId, // id of current stripe
							rvHE_BackupVector, // HE-backup vector
							iLeftUpH_Value, // set by reference in method
							xStripeData.first, //// rxQueryProfile.xDeviceVector[uiStripeId],
						    uiStripWidth ); //// TILE_WIDTH ); // number of elements in stripe
			
			/* The HF-vectors must be swap after each call of updateTile */
			bHF_VectorsReversed = !bHF_VectorsReversed;
		} // for

#if ( DO_TESTS == 1 )
		//// std::cout << "**** Check last column:" << std::endl;
		//// /* If the divisor is greater than 1, we get problem here for all segments following the first one */
		//// this->debugCheckHF_VectorTwo( pNaiveSW, QUERY_SIZE - 1, !bHF_VectorsReversed );
#endif
	} // method


	void doQueryForChunkedReference( std::vector<char> &rvQuerySeq, // query
									 const ChunkedTransposedReference &rxChunkedRef, // reference in segmented form
									 std::vector<size_t> &rvMaxScorePositions, // vector receives maximum positions
									 SCORE_TP &iOverallMaxScore )  // vector that logs the maximum-score positions
	{	
		assert( rxChunkedRef.uiChunkSize == this->uiRefCapicity );
		assert( rxChunkedRef.uiSegmentSize == this->uiSegmentSize );
		
		/* Create the query profile */
		QueryProfile<SCORE_TP4, SCORE_TP, STRIP_WIDTH> xQueryProfile( rvQuerySeq.size() );
		xQueryProfile.set( rvQuerySeq );

		/* According to the standard the vector is zero initialized.
		 * See: https://stackoverflow.com/questions/25198405/are-members-of-structs-in-a-vector-initialized-with-zero-in-c
		 */
		std::vector<SCORE_TP2> vHE_BackupVector( xQueryProfile.getNumberOfStripes() * STRIP_WIDTH );

		iOverallMaxScore = 0;

		//// /* Receives and keeps the maximum scores after each kernel call */
		//// std::vector<SCORE_TP> vM_Buffer_Vector;

		/* Iterate over all reference chunks. */
		for( size_t uiChunkId = 0; uiChunkId < rxChunkedRef.uiNumberOfChunks; uiChunkId++ )
		{	/* Load the required reference segment into device */
			/* TO DO: We clear to many vectors over here. */
			this->clearDeviceVectors();

			/* Get the appropriate section of the reference in the device memory */
			metaMeasureAndLogDuration<false>
			(	"Time for loading reference",	// text message
				[&] () // lambda by reference
				{
					//// this->loadReference( *(rxChunkedRef.vChunks[uiChunkId]) );
					this->loadReference( *(rxChunkedRef.getChunk( uiChunkId ) ) );
				} // lambda
			); // function call

			/* Process the chunk via the kernel */
			metaMeasureAndLogDuration<false>
			(	"Time for kernel execution", // text message
				[&] () // lambda by reference
				{
					this->doQueryForChunk( xQueryProfile,
#if ( DO_TESTS == 1 )
												rvQuerySeq, // for test we require the original query sequence
#endif
												vHE_BackupVector );
				} // lambda
			); // function call

			/* Extract the maximum score and positions of the maximum score */
			this->xM_Vector.updateSharedBackup();
			
#if (USE_THREADPOOL == 1)
			if( this->xMaxExtractSuccess.valid() )
			{	/* The future is active and we have to block until 
				 * work of the previous maximum extraction is done 
				 */
				bool wait = this->xMaxExtractSuccess.get();
			} // if
			
			/* Set the future by queuing in the Threadpool */
			this->xMaxExtractSuccess
				= this->xThreadPool.enqueue( [&] ( size_t uiThreadId, size_t uiChunkIdLoc )
				{
					//// std::cout << "Thread " << uiThreadId << " starts max extraction" << std::endl;
					this->extractMaxima( rvMaxScorePositions,
										 iOverallMaxScore,
										 rxChunkedRef.numOfSymUsedInChunk( uiChunkIdLoc ),
										 uiChunkIdLoc * rxChunkedRef.uiChunkSize );
					//// std::cout << "Thread " << uiThreadId << " ends max extraction" << std::endl;
					return true;
				}, // lambda
				
				uiChunkId
			); // enqueue
#else			
			metaMeasureAndLogDuration<false>
			(	"Time for extracting maximum scores",	// text message
				[&] () // lambda by reference
				{
					this->extractMaxima( rvMaxScorePositions,
										 iOverallMaxScore,
										 rxChunkedRef.numOfSymUsedInChunk( uiChunkId ),
										 uiChunkId * rxChunkedRef.uiChunkSize );
				} // lambda
			); // function call
		
#endif
		} // for

#if (USE_THREADPOOL == 1)
		/* We have to wait until the final maximum extraction finished */
		if( this->xMaxExtractSuccess.valid() )
		{	/* The future is active and we have to block until work is done */
			//// std::cout << "WAIT FOR FINISHING MAX EXTRACTION" << std::endl;
			bool wait = this->xMaxExtractSuccess.get();
		} // if
		//// std::cout << "MAX EXTRACTION DONE" << std::endl;
#endif

		/* Sort the vector that comprises all maximum positions */
		std::sort( rvMaxScorePositions.begin(), rvMaxScorePositions.end() );
	} // method


	/* Collects the maxima contained in the maximum-vector.
	 * Optimization: Do the maximum collection of the previous kernel call during the next kernel call.
	 */
	void extractMaxima( std::vector<size_t> &rvMaxScorePositions, // vector keeping max positions
						SCORE_TP &iOverallMaxScore, // the maximum values
						const size_t uiNumOfSymUsedInChunk, // number of symbols actually used in chunk
						const size_t offset ) // reference offset for calculating the original position
	{
		assert( uiNumOfSymUsedInChunk > 0 && uiNumOfSymUsedInChunk <= this->uiRefCapicity );

		/* Fetch the ME-vector from the device */
		auto pME_VectorBackup = this->xM_Vector.getSharedBackup();
		auto &rvME_VectorBackup = *pME_VectorBackup;

		/* Find all maximum positions (done on CPU side).
		* (The maximum search could be done efficiently by a reduce kernel on GPU side)
		*/
		for( size_t uiItr = 0; uiItr < this->uiRefCapicity; uiItr++ )
		{
			if( rvME_VectorBackup[uiItr] >= iOverallMaxScore )
			{
				size_t uiSymPosInChunk = rowIdOptimizedToStandard( uiItr );
				if( uiSymPosInChunk >= uiNumOfSymUsedInChunk )
				{	/* Position is out of range */
					continue;
				} // if
				
				if( rvME_VectorBackup[uiItr] > iOverallMaxScore )
				{	// fresh overall maximum detected 
					rvMaxScorePositions.clear();
					iOverallMaxScore = rvME_VectorBackup[uiItr];
				} // if
#if	( OPTIMIZED_INDEXING == 0 )
				rvMaxScorePositions.push_back( this->rowIdOptimizedToStandard( uiRow ) );
#else
				//// rvMaxScorePositions.push_back( offset + uiItr );
				rvMaxScorePositions.push_back( offset + uiSymPosInChunk );
#endif // OPTIMIZED_INDEXING
			} // if
		} // for
	} // method

#if	( OPTIMIZED_INDEXING == 1 )
	/* Translates a optimized row coordinate to the real value */
	inline size_t rowIdOptimizedToStandard( size_t uiRowOptimized )
	{
		size_t uiSegmentPos = uiRowOptimized / this->uiNumberOfSegments;
		size_t uiSegmentId = uiRowOptimized % this->uiNumberOfSegments;
		return uiSegmentId * uiSegmentSize + uiSegmentPos;
	} // method

	/* Translates a real row coordinate to the optimized value */
	inline size_t rowIdStandardToOptimzied( size_t uiRowStandard )
	{
		size_t uiSegmentId = uiRowStandard % this->uiSegmentSize;
		size_t uiSegmentPos = uiRowStandard / this->uiSegmentSize;
		return uiSegmentId * uiNumberOfSegments + uiSegmentPos;
	} // method
#endif

#if ( DO_TESTS == 1 )
	/* Fetch a HF-vector from device */
	auto get_HF_Vector( bool bChooseVectorOne )
	{
		return bChooseVectorOne ? this->xHF_VectorOne.getCopyAsVector() : this->xHF_VectorTwo.getCopyAsVector();
	} // method

	void debugCheckHF_VectorTwo( std::shared_ptr<NaiveSW> pNaiveSW,
								 size_t uiColumn,
								 bool bChooseVectorOne )
	{
#if ( USE_PACKED == 1 )
		auto vDeviceVector = reinterpret_cast<std::vector<int2> &>(this->get_HF_Vector( bChooseVectorOne ));
		vDeviceVector = SequenceTransformer<int2>::inverseTransposedSeq( vDeviceVector, this->uiSegementSize );
		/* Doubles the size of the vector, so that it is equal to original reference size*/
		vDeviceVector = SequenceTransformer<int2>::unpackint2_2Vertically( vDeviceVector );
#else
		auto vDeviceVector = this->get_HF_Vector( bChooseVectorOne );
		assert( vDeviceVector.size() == CHUNK_SIZE );
		assert( this->uiRefCapicity == CHUNK_SIZE );
		vDeviceVector = SequenceTransformer<SCORE_TP2>::inverseTransposedSeq( vDeviceVector, this->uiSegementSize );
#endif

		auto vNaiveSW_HRow = pNaiveSW->getH_Row( uiColumn );
		for( size_t uiRow = 0; uiRow < this->uiRefCapicity; uiRow++ )
		{
#if	( OPTIMIZED_INDEXING == 0 )
			auto uiRowOnDeviceVector = this->rowIdStandardToOptimzied( uiRow );
#else
			auto uiRowOnDeviceVector = uiRow;
#endif // OPTIMIZED_INDEXING
			if( vDeviceVector[uiRowOnDeviceVector].x != vNaiveSW_HRow[uiRow] )
			{
				std::cout << "Diff at row " << uiRow << " column " << uiColumn << " value SW " 
					<< vNaiveSW_HRow[uiRow] << " value GPU " << static_cast<int>(vDeviceVector[uiRowOnDeviceVector].x)
					<< " relative_in_block " << (uiRow % this->uiSegementSize) << std::endl;
				exit( 0 );
			} // if
		} // for
	} // method
#endif // DO_TESTS

	/* Returns if the kernel in lazy mode could fix all inaccuracies.
	 * If we get true over here, the lazy kernel does not need to be called once more. 
	 */
	bool lazyFixed()
	{	/* Copy the lazy fixed-vector to host memory */
		this->pLazyFixedVector.copyToArray( this->pContinuationFlags.get() );
		auto pFlagsAsLong = this->pContinuationFlags.get();

		/* WARNING!: In order to get the below code working, the vector has to be defined of type long
		 * or you risk misalignments
		 */
		//// long* pFlagsAsLong = reinterpret_cast<long *>(this->pContinuationFlags.get());
		
		for( size_t uiCounter = 0; uiCounter < this->uiNumberOfSegments; ++uiCounter )
		{
			if( pFlagsAsLong[uiCounter] != 0 )
			{	/* Some segment indicates the need for recall of lazy-Checksum-loop */
				return false;
			} // if
		} // for
		return true;
	} // method	

	/* Destructor */
	~SW_GPU_Processor()
	{} // destructor
}; // class


/* Main function for ongoing CUDA work.
 * Returned value: Maximum score.
 * 32 bit scoring: <int4, int2, int, int>
 * 24 bit scoring: <float4, float2, float, float>
 */
template<typename SCORE_TP4, typename SCORE_TP2, typename SCORE_TP, typename CHECKSUM_TP>
std::vector<GPUReturn> cudaAlignTmpl
(   std::vector<char> &rvRefSeq, // reference sequence
	std::vector<std::vector<char>> &rvQuerySeqs, // vector of query sequences
    unsigned int uiDeviceId
) 
{	/* Do checks and reset device */
	assert( sizeof( SCORE_TP4 ) == 4 * sizeof( SCORE_TP ) && sizeof( SCORE_TP2 ) == 2 * sizeof( SCORE_TP ) );

    //// for(auto& rQ : rvQuerySeqs)
        for(char c :rvRefSeq)
            if(c >= 4)
                std::cout << "WARNING " << c << std::endl;
	
	cudaSetDevice( uiDeviceId );
	/* Reset all GPU devices */
	cudaDeviceReset();
	
	/* The raw size of the reference.*/
	size_t uiRefSize = rvRefSeq.size();

	size_t uiRefChunkCapicity = CHUNK_SIZE; // depends on the memory available
	SW_GPU_Processor<SCORE_TP4, SCORE_TP2, SCORE_TP, CHECKSUM_TP, TILE_WIDTH> xSW_GPU_Processor
	(	uiDeviceId, // device id 
		uiRefChunkCapicity ); // capacity for single reference chunk
	
	/* Set the segment size according to the freshly loaded reference.
	 * (At this point we know that the reference-size is a multiple of the number of segments.)
	 * TO DO: Set this as part of the GPU initialization.
	 */
	xSW_GPU_Processor.uiSegmentSize = CHUNK_SIZE / xSW_GPU_Processor.uiNumberOfSegments; 

	/* Create the segmented reference.
	 * For queries against the same reference this has to be done only once!
	 */
	ChunkedTransposedReference xChunkedReference( rvRefSeq, // reference sequence
												  uiRefChunkCapicity, // primary segment size (should be equal to strip capacity)
												  xSW_GPU_Processor.uiNumberOfSegments ); // secondary segment size
	
    std::vector<GPUReturn> vRet; 
    //std::vector<std::pair<SCORE_TP, std::vector<size_t>>> vRet;

	for( auto &rvQuerySeq : rvQuerySeqs )
    {
        SCORE_TP iOverallMaxScore = 0;
        std::vector<size_t> vMaxScorePositions;

	    /* Do the core alignment */
	    metaMeasureAndLogDuration<false>
	    (	"GPU time",	// logging text
	    	[&] () // lambda by reference
	    	{
	    		xSW_GPU_Processor.doQueryForChunkedReference
	    		( rvQuerySeq, // query
	    		  xChunkedReference, // reference in segmented form
	    		  vMaxScorePositions, // vector that receives maximum positions
	    		  iOverallMaxScore ); // overall maximum score
	    	} // lambda
	    ); // function call

        //// std::cout << "iOverallMaxScore: " << iOverallMaxScore << std::endl;
        vRet.emplace_back( iOverallMaxScore, vMaxScorePositions );
    } // for (iterating queries)

	return vRet;
} // function

size_t getNumberOfCards()
{
	int iDeviceCount; 

	cudaGetDeviceCount( &iDeviceCount );
	CUERR

	return iDeviceCount;
} // function


std::vector<GPUReturn> cudaAlign
(
    std::vector<char> &rvRefSeq, // reference sequence
	std::vector<std::vector<char>> &rvQuerySeqs, // vector of query sequences
    unsigned int uiDeviceId
)
{
    return cudaAlignTmpl<int4, int2, int, int>(rvRefSeq, rvQuerySeqs, uiDeviceId);
};

#if 0
std::vector<GPUReturn> cudaAlignPool
(
	std::vector<char> &rvRefSeq, // reference sequence
	std::vector<std::vector<char>> &rvQuerySeqs // vector of query sequences
)
{
	auto uiNumberOfCards = getNumberOfCards();
	std::cout << "Number of cards: " << uiNumberOfCards << std::endl;

	ThreadPool xThreadPool(uiNumberOfCards);

	std::vector< std::future<GPUReturn> > results;

	for ( size_t i = 0; i < rvQuerySeqs.size(); i++ )
	{	/* i must be passed by value into the lambda */
		results.push_back( xThreadPool.enqueue(
			[i] ( size_t uiThreadId )
			{
				cudaAlign()
			//// std::cout << "Thread " << uiThreadId << " starts max extraction" << std::endl;
				this->extractMaxima( rvMaxScorePositions,
									 iOverallMaxScore,
									 rxChunkedRef.numOfSymUsedInChunk( uiChunkIdLoc ),
									 uiChunkIdLoc * rxChunkedRef.uiChunkSize );
				//// std::cout << "Thread " << uiThreadId << " ends max extraction" << std::endl;
				return true;
			} // lambda
		)); 
		
			for ( size_t i = 0; i < results.size(); ++i )
				BOOST_LOG_TRIVIAL( trace ) << "*** See result for " << i << " as " << results[i].get();
		// std::cout << results[i].get() << ' ';



		xThreadPool.enqueue( [i] ( size_t id, int j ) { itemWorker( i, id, j ); }, 6 );
	} // for

	this->xMaxExtractSuccess
		

									 uiChunkId
		); // enqueue

	std::future<GPUReturn> xMaxExtractSuccess;


	
	// return cudaAlignTmpl<int4, int2, int, int>(rvRefSeq, rvQuerySeqs);
	return std::vector<GPUReturn>();
};
#endif
