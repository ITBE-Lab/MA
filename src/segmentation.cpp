#include "segmentation.h"

#define INCLUDE_SELF_TESTS (  1 )

//#include "BWTmem.h"
#include <vector>
//#include "analysis.h"
#include <memory>
#include <atomic>
#include <chrono>
//#include "assembly.h"
 

/**
* Delivers 4 intervals for a single input interval.
* Here we use only two fields in the BWT_Interval.
*/
SA_IndexInterval bwt_extend_backward(
									// current interval
									const SA_IndexInterval &ik,
									// the character to extend with
									const uint8_t c,
									// the FM index by reference
									std::shared_ptr<FM_Index> pxFM_Index
									)
{
	bwt64bitCounter cntk[4]; // Number of A, C, G, T in BWT until start of interval ik
	bwt64bitCounter cntl[4]; // Number of A, C, G, T in BWT until end of interval ik

	pxFM_Index->bwt_2occ4(
		// until start of SA index interval (-1 due to $ in BWT)
		ik.start() - 1,
		// until end of SA index interval (end is inclusive => -1)
		ik.end() - 1,
		cntk,						// output: Number of A, C, G, T until start of interval
		cntl						// output: Number of A, C, G, T until end of interval
	);
	//pxFM_Index->L2[c] start of nuc c in BWT
	//cntk[c] offset of new interval
	//cntl[c] end of new interval
	return SA_IndexInterval(pxFM_Index->L2[c] + 1 + cntk[c], cntl[c] - cntk[c]);
} // method

bool Segmentation::canExtendFurther(std::shared_ptr<SegmentTreeInterval> pxNode, nucSeqIndex uiCurrIndex, bool bBackwards, nucSeqIndex uiQueryLength)
{
	if (!bBackwards)
		return uiCurrIndex < uiQueryLength;
	return true;
}//function

SaSegment Segmentation::extend(
								 std::shared_ptr<SegmentTreeInterval> pxNode, 
								 nucSeqIndex uiStartIndex, bool bBackwards
								)
{
	nucSeqIndex i = uiStartIndex;
	assert(i >= 0);
	assert(i < pxQuerySeq->length());

	const uint8_t *q = pxQuerySeq->pGetSequenceRef(); // query sequence itself

	//select the correct fm_index based on the direction of the extension
	std::shared_ptr<FM_Index> pxUsedFmIndex;
	if (bBackwards)
		pxUsedFmIndex = pxFM_index;
	else
		pxUsedFmIndex = pxRev_FM_Index;
	
	/* Initialize ik on the foundation of the single base q[x].
	 * In order to understand this initialization you should have a look 
	 *to the corresponding PowerPoint slide.
	 */
	// start I(q[x]) in T (start in BWT used for backward search) + 1, 
	// because very first string in SA-array starts with $
	// size in T and T' is equal due to symmetry
	SA_IndexInterval ik(
						pxUsedFmIndex->L2[(int)q[i]] + 1, 
						pxUsedFmIndex->L2[(int)q[i] + 1] - pxUsedFmIndex->L2[(int)q[i]]
					);

	if (!bBackwards)
		i++;
	else if (i == 0)// unsigned int
		return SaSegment(0,0,ik,true);
	else
		i--;

	// while we can extend further
	while ( canExtendFurther(pxNode ,i ,bBackwards ,pxQuerySeq->length()) )
	{
		assert(i >= 0 && i < pxQuerySeq->length());
		if (q[i] >= 4 && bBreakOnAmbiguousBase) // An ambiguous base
			break; // break if the parameter is set.

		const uint8_t c = q[i]; // character at position i in the query

		// perform one step of the extension
		SA_IndexInterval ok = bwt_extend_backward(ik, c, pxUsedFmIndex);

		std::cout << i << " " << ik.size() << " -> " << ok.size();
		std::cout << " (" << ok.start() << ", " << ok.end() << ")" << std::endl;
			
		/*
		* In fact, if ok.getSize is zero, then there are no matches any more.
		*/
		if (ok.size() == 0)
			break; // the SA-index interval size is too small to be extended further
			
		
		if (!bBackwards)
			i++;
		else if (i-- == 0)// unsigned int
			break;

		/* Set input interval for the next iteration.
		*/
		ik = ok;

	}//while

	if (bBackwards)
		return SaSegment(i + 1, uiStartIndex - (i + 1), ik, false);
	else
		return SaSegment(uiStartIndex, i - uiStartIndex - 1, ik, true);
}//function

/* this function implements the segmentation of the query
 *
 * the process is synchronized using a thread pool
 *
 * to avoid unnecessary locking and unlocking we make sure that only one thread can process one interval at any given time.
 * ( there are always several intervals ready and waiting to be processed. )
 * this way it is only necessary to lock once we touch the structure of the list which stores the individual intervals.
 *
 * segmentation technique:
 *		start in the middle of the interval try to extend in both directions
 *		split the interval in 3 parts: prev:perfectMatch:post
 *		queue the prev and post intervals into the thread pool
 *		save the perfect match for later clustering
*/
void Segmentation::procesInterval(size_t uiThreadId, SegTreeItt pxNode, ThreadPoolAllowingRecursiveEnqueues *pxPool)
{
	DEBUG(
		std::cout << "interval (" << pxNode->start() << "," << pxNode->end() << ")" << std::endl;
	)

	//performs forward extension and records any perfect matches
	SaSegment xForw = extend(*pxNode, pxNode->getCenter(), false);
	/* we could already extend our matches from the center of the node to uiForwardExtensionReachedUntil
	* therefore we know that the next extension will reach at least as far as the center of the node.
	* we might be able to extend our matches further though.
	*/
	SaSegment xForwBack = extend(*pxNode, xForw.end(), true);
	assert(xForwBack.start() <= pxNode->getCenter());
	
	//performs backwards extension and records any perfect matches
	SaSegment xBack = extend(*pxNode, pxNode->getCenter(), true);
	/* we could already extend our matches from the center of the node to  uiBackwardsExtensionReachedUntil
	* therefore we know that the next extension will reach at least as far as the center of the node.
	* we might be able to extend our matches further though.
	*/
	SaSegment xBackForw = extend(*pxNode, xBack.start(), false);
	std::cout << xBackForw.end() << std::endl;
	std::cout << pxNode->getCenter() << std::endl;
	assert(xBackForw.end() >= pxNode->getCenter());

	nucSeqIndex uiFrom, uiTo;
	if (xBackForw.size() > xForwBack.size())
	{
		uiFrom = xBackForw.start();
		uiTo = xBackForw.end();
		pxNode->push_back(xBackForw);
	}//if
	else
	{
		uiFrom = xForwBack.start();
		uiTo = xForwBack.end();
		pxNode->push_back(xForwBack);
	}//else

	if (uiFrom != 0 && pxNode->start() + 1 < uiFrom)
	{
		//create a new list element and insert it before the current node
		auto pxPrevNode = pSegmentTree->insertBefore(std::shared_ptr<SegmentTreeInterval>(
			new SegmentTreeInterval(pxNode->start(), uiFrom - pxNode->start() - 1)), pxNode);
		//enqueue procesInterval() for the new interval
		pxPool->enqueue(
			[/*WARNING: do not catch anything here: the lambda function is enqueued to into a thread pool, 
			 local variables might not exist anymore during it's execution*/]
			(size_t uiThreadId, SegTreeItt pxPrevNode, Segmentation *pxAligner, ThreadPoolAllowingRecursiveEnqueues *pxPool)
			{
				pxAligner->procesInterval(uiThreadId, pxPrevNode, pxPool);
			},//lambda
			pxPrevNode, this, pxPool
		);//enqueue
	}//if
	if (pxNode->end() > uiTo + 1)
	{
		//create a new list element and insert it after the current node
		auto pxNextNode = pSegmentTree->insertAfter(std::shared_ptr<SegmentTreeInterval>(
			new SegmentTreeInterval(uiTo + 1, pxNode->end() - uiTo - 1)), pxNode);
		//enqueue procesInterval() for the new interval
		pxPool->enqueue(
			[/*WARNING: do not catch anything here: the lambda function is enqueued to into a thread pool,
			 local variables might not exist anymore during it's execution*/]
			(size_t uiThreadId, SegTreeItt pxNextNode, Segmentation *pxAligner, ThreadPoolAllowingRecursiveEnqueues *pxPool)
			{
				pxAligner->procesInterval(uiThreadId, pxNextNode, pxPool);
			},//lambda
			pxNextNode, this, pxPool
		);//enqueue
	}//if


	DEBUG(
		std::cout << "splitting interval (" << pxNode->start() << "," << pxNode->end() << ") at (" << uiFrom << "," << uiTo << ")" << std::endl;
	)
	pxNode->start(uiFrom);
	pxNode->end(uiTo);
}//function

void Segmentation::segment()
{
	assert(*pSegmentTree->begin() != nullptr);

	{//scope for xPool
		ThreadPoolAllowingRecursiveEnqueues xPool( 1 );//NUM_THREADS_ALIGNER);

		//enqueue the root interval for processing
		xPool.enqueue(
			[/*WARNING: do not catch anything here: the lambda function is enqueued to into a thread pool,
			 local variables might not exist anymore during it's execution*/]
			(size_t uiThreadId, SegTreeItt pxRoot, Segmentation *pxAligner, ThreadPoolAllowingRecursiveEnqueues *pxPool)
			{
				pxAligner->procesInterval(uiThreadId, pxRoot, pxPool);
			},//lambda
			pSegmentTree->begin(), this, &xPool
		);//enqueue

	}//end of scope xPool
}//function


std::vector<ContainerType> SegmentationContainer::getInputType()
{
	return std::vector<ContainerType>{
			//the forward fm_index
			ContainerType::fM_index,
			//the reversed fm_index
			ContainerType::fM_index,
			//the query sequence
			ContainerType::nucSeq,
			//the reference sequence (packed since it could be really long)
			ContainerType::packedNucSeq,
		};
}
ContainerType SegmentationContainer::getOutputType()
{
	return ContainerType::segmentList;
}


std::shared_ptr<Container> SegmentationContainer::execute(
		std::vector<std::shared_ptr<Container>> vpInput
	)
{
	std::shared_ptr<FM_Index> pFM_index = std::static_pointer_cast<FM_Index>(vpInput[0]);
	std::shared_ptr<FM_Index> pFM_indexReversed = std::static_pointer_cast<FM_Index>(vpInput[1]);
	std::shared_ptr<NucleotideSequence> pQuerySeq = 
		std::static_pointer_cast<NucleotideSequence>(vpInput[2]);
	std::shared_ptr<BWACompatiblePackedNucleotideSequencesCollection> pRefSeq = 
		std::static_pointer_cast<BWACompatiblePackedNucleotideSequencesCollection>(vpInput[3]);


	Segmentation xS(
			pFM_index, 
			pFM_indexReversed, 
			pQuerySeq, 
			bBreakOnAmbiguousBase,
			pRefSeq
		);
	xS.segment();

	return xS.pSegmentTree;
}//function

/**
 * @brief for analyzing a CRISPER experiment
 * @details
 * quick and dirty code that is used to analyze sequences for a CRISPER experiment.
 */
void analyseCRISPER()
{
	std::shared_ptr<FM_Index> xIndex(new FM_Index());

	/*BWACompatiblePackedNucleotideSequencesCollection xPack;
	std::string pref("/mnt/ssd0/chrom/human/");
	for(unsigned int i = 1; i <= 22; i++)
		xPack.vAppendFASTA(std::string(pref).append("chr").append(std::to_string(i)).append(".fna"));
	xPack.vAppendFASTA(std::string(pref).append("chrX.fna"));
	xPack.vAppendFASTA(std::string(pref).append("chrY.fna"));
	xIndex = std::shared_ptr<FM_Index>(new FM_Index(xPack));
	xIndex->vStoreFM_Index("/mnt/ssd0/chrom/human/index");*/
#if 0
	xIndex->vLoadFM_Index("/mnt/ssd0/chrom/mouse/index");
	std::string searchFor[] = 
	{
		/*
		"CGTCAGACTTACTTGGTAGG",
		"ACAGAATTTGCAACACAGGA",
		"ACCCGTCAGACTTACTTGGT",
		"GTCAGACTTACTTGGTAGGA",
		"GTTTACAGAATTTGCAACAC"
		*/

		"TACCAAGTAAGTCTGACGGG",
		"GTGGCCCTTTAGCTGGACAG",
		"ATTCTGCACCTCTTTCAACG",
		"AACACAGGATGGATACGCAT",
		"TTCTGCACCTCTTTCAACGT",
		"CGTCAGACTTACTTGGTAGG",
		"ATAGCGTCCCACGTTGAAAG",
		"GAGGTGCAGAATTTTGAAAG",
		"ACCCGTCAGACTTACTTGGT",
		"CTCCTACCAAGTAAGTCTGA",
		"AGCTGGACAGTGGATGGCAT",
		"TAGCTGGACAGTGGATGGCA",
		"TCTCAGCTCTAGGGGAATGA",
		"GCCATCCACTGTCCAGCTAA",
		"TCCTACCAAGTAAGTCTGAC",
		"GTCAGACTTACTTGGTAGGA",
		"ACAGAATTTGCAACACAGGA",
		"AAACAAGTTCTCAGCTCTAG",
		"GGCCACCCGTCAGACTTACT",
		"GTTTACAGAATTTGCAACAC",
		"CCCTTTAGCTGGACAGTGGA",
		"GTCTTTGAAAGTACAGTTGT",
		"GTAAACAAGTTCTCAGCTCT",
		"CCATCCACTGTCCAGCTAAA",
		"GAGGGATGTCTGTAGAGAAA",
		"AAGAGGTGCAGAATTTTGAA",
		"AGAGGTGCAGAATTTTGAAA",
		"TAAACAAGTTCTCAGCTCTA",
		"CTGACGGGTGGCCCTTTAGC",
		"GGGAATGAAGGCTGTTTTGC"
	};
#else
	xIndex->vLoadFM_Index("/mnt/ssd0/chrom/human/index");
	std::string searchFor[] = 
	{
		/*
		"ACAGAATTTGCAACACAGGA",
		"GTCCTTGAAACTACAAATGT",
		"GTTTACAGAATTTGCAACAC"
		*/

		/*
		"GGTGTAAACGTCCTACTTTG",
		"AAGTAGGACGTTTACACCTG",
		"TGTGTTAAATATCTCAGCAG",
		"ATTAACAGACTAGACTCCAC",
		"AATCTGCAAATGTAAACTGA",
		"GTTGTTGCAGGAAAGAAGAC",
		"GAAAGAAGACAGGAACTACA",
		"ACTACAAATGTTGGTAGACA",
		"ACAGAATTTGCAACACAGGA",
		"GCTCAGCTCTAGGGGAATGA",
		"GACATCGCATGTTTAAAATG",
		"CTGGCAAGTTGCCCCAAAGT",
		"TGTCTTCTTTCCTGCAACAA",
		"TGTTAATATCATAAACTGCC",
		"GTAAACAAGTGCTCAGCTCT",
		"CTTCTGGATGTTTCACAGCC",
		"ACATTTGTAGTTTCAAGGAC",
		"GTCCTTGAAACTACAAATGT",
		"GTTTACAGAATTTGCAACAC",
		"CAGGTGTAAACGTCCTACTT",
		"GTTGGTAGACATGGCGTCGC",
		"GAGGGATGTCTGTAGAGAAA"
		*/

		/*
		"TCTGTACTCTCTTCTATCTC",
		"AGTCCCCTCAGGTTTCTTCC",
		"CGGACCTTCCTCCGCTCTGC",
		"CGGGCGCGCGGTGGTGTTGG",
		"GCGCGAGCACGCTGAGGGCC",
		"TTTCCCTTTCTGCTGTCTAA",
		"CAGAGCTGAAGTCGCTTTAG",
		"GCTGAGGGCCGGGGGTTGCC",
		"AGTAGGGTAGTGTTGGATTC",
		"GGGCGCGCGGTGGTGTTGGT",
		"CGCGAGCACGCTGAGGGCCG",
		"GCCCCTCCCCCGGAGCGCTC",
		"AATGTAAAATTAATGGTATC",
		"TGGGCATAACAAAAAATTAA",
		"CCAGAAGTGTCCGTTGTTGC",
		"GAAGGGAGTGGAGTTTGAAT",
		"TAAACAAGTGCTCAGCTCTA",
		"TTAGAGGCTCAGTCGGCGCT",
		"TGGGTACTAGGCCTTTATTC",
		"GCAGGGAACCGGAGCGCTCC",
		"GGAGAAAAGGGGCGGGACGG",
		"CTTTTTCTGCCAGGCGATCC",
		"AGCAGGGAACCGGAGCGCTC",
		"CGCAGAGCGGAGGAAGGTCC",
		"CTGGTTCTCATAGGTTCTAT",
		"GAGGAAGGTCCGGGAGAAAA",
		"AGGAGGCAGGTGGAGTTTGA",
		"CACTCAGTTACCTCCCAATC",
		"GGCTGAGGGAATTGCTATGA",
		"AGTTCAACCGAGTTTTAATA",
		"AGGGGACTGCGGGGGGCACG",
		"CTCCCTTCTCTTTTTCTGCC",
		"TTAAAATGTGGGAACATTGT",
		"AGCGCGAGCACGCTGAGGGC",
		"CCGCAGAGCGGAGGAAGGTC",
		"AGAGGAGTGTATTAGATAAA",
		"CTCTTTAATATCTCAGTTGT",
		"GGGCCGGGGGTTGCCAGGAG",
		"CCTACTTTCCTATTCCTTTG",
		"ACCGGAGCGCTCCGGGGGAG",
		"CCTGCAACAACGGACACTTC",
		"GGTTCCCTGCTGTCGGATCT",
		"AGTTTTTTTTTCTGTTTTCC",
		"GGGTCCGCGGGCCCCTCTCC",
		"AGGTGTAAACGTCCTACTTT",
		"CACCAACACCACCGCGCGCC",
		"GCGCCCGGGAACAGCCGCTC",
		"GGGAATGAAGGCTGTTTTGC",
		"AGTATTTATTGAGTAGGAAC",
		"TCCCTTCTCCAAAAGTTAAA",
		"GCTTTTCTAAGTCCTAATAT",
		"TTTAGAGGCTCAGTCGGCGC",
		"TCTTTTTTGGAAAACAAAAA",
		"GGGCGGGACGGAGGAGAATC",
		"ACATTCTTCTAAAGCCATTA",
		"GGAGGAAGGTCCGGGAGAAA",
		"CTTCCCGGAGCGGCTGTTCC",
		"TACCAACATTTGTAGTTTCA",
		"TTTTTTTCTATTTAAGATTT",
		"GAACCGGAGCGCTCCGGGGG",
		"CTTTTCTAAGTCCTAATATT",
		"CTTTTTTGGAAAACAAAAAT",
		"TTGGCTTTCCTTTTAACTTT",
		"CAGTTTACATTTGCAGATTT",
		"CCGGGGGAGGGGCTTAATGC",
		"AAGTTCAACCGAGTTTTAAT",
		"CTAATTTCCTTCCCAATATT",
		"ATTCAGGAATGTAAAATTAA",
		"GGGAACAGGATAACTTATTT",
		"TCTATTTAAAAATATATTAA",
		"TAGTTGCCAGTGATCTTTTT"
		*/

		/*
		"CACTGCGCTGGAGAGAGAAA",
		"CCAGTGAAGTCCACTGCGCT",
		"CCTTTCTCTCTCCAGCGCAG",
		"TCCAGCGCAGTGGACTTCAC",
		"AGTTTACATTTGCAGATTTT",
		"AATCTGCAAATGTAAACTGA",
		"AGGACATCGCATGTTTAAAA",
		"ATGTTTAAAATGGTAATTTG",



		"ACATTTACTTTCTTTTCTAG",
		"CATTTACTTTCTTTTCTAGT",
		"TTTTCTAGTGGGAACATTGT",
		"AAAGAAGACAGGAACTACAT",
		"TTGTTGCAGGAAAGAAGACA",
		"CAGAAGTGTCCGTTGTTGCA",
		"TGTCTTCTTTCCTGCAACAA",
		"CCTGCAACAACGGACACTTC",
		"GTTAATATCATAAACTGCCT",
		"CTTCTGGATGTTTCACAGCC",
		"GTCAATTGTTTCTTACCTGT",
		"ATTAACAGACTAGACTCCAC",
		"CAGGTAAGAAACAATTGACT",



		"TTTTATGACTTTTTATTTTC",
		"TGGCAAGTTGCCCCAAAGTA",
		"CAGGTGTAAACGTCCTACTT",
		"AGGTGTAAACGTCCTACTTT",
		"GGTGTAAACGTCCTACTTTG",
		"CTGTACTCTCTTCTATCTCT",
		"TGTGTTAAATATCTCAGCAG",
		"AGAGGAGTGTATTAGATAAA"
		*/

				
		"TACCAAGTAAGTCATGCTAG",
		"AGCATGACTTACTTGGTAGG",
		"CAAGTAAGTCATGCTAGTGG",
		"ACTAGCATGACTTACTTGGT",
		"GTCGCAGGAAGGATGAGGTG",
		"GTAGACATGGCGTCGCAGGA",
		"TCGCAGGAAGGATGAGGTGT",
		"ACTACAAATGTTGGTAGACA",
		"ATGGCGTCGCAGGAAGGATG",
		"AAGTAAGTCATGCTAGTGGC",
		"ACAGAATTTGCAACACAGGA",
		"GCTCAGCTCTAGGGGAATGA",
		"GCATGACTTACTTGGTAGGA",
		"AAACAAGTGCTCAGCTCTAG",
		"GTAAACAAGTGCTCAGCTCT",
		"ACATTTGTAGTTTCAAGGAC",
		"GTCCTTGAAACTACAAATGT",
		"GTTTACAGAATTTGCAACAC",
		"GTTGGTAGACATGGCGTCGC",
		"GAGGGATGTCTGTAGAGAAA",
		"GAGGTGTGGGATTTTGAAAA",
		"CGCCACTAGCATGACTTACT",
		"TAAACAAGTGCTCAGCTCTA",
		"GGGAATGAAGGCTGTTTTGC",
		"TACCAACATTTGTAGTTTCA"
	};
#endif
	std::vector<std::vector<std::tuple<std::string, unsigned int, std::list<int64_t>>>> results;
    {
        ThreadPool xPool(48);		

		std::mutex m;
		unsigned int pos = 0;

        for(std::string& sequence : searchFor)
        {
				results.push_back(std::vector<std::tuple<std::string, unsigned int, std::list<int64_t>>>());
                for(unsigned int i = 0; i < sequence.size(); i++)
                {
                for(unsigned int i_mut = (i == 0 ? 0 : 1); i_mut < 4; i_mut++)
                {
                        for(unsigned int j = i+1; j < sequence.size(); j++)
                        {
                        for(unsigned int j_mut = (j == i+1 && i_mut == 0 ? 0 : 1); j_mut < 4; j_mut++)
                        {
                                for(unsigned int k = j+1; k < sequence.size(); k++)
                                {
                                for(unsigned int k_mut = (k == j+1 && j_mut == 0  ? 0 : 1); k_mut < 4; k_mut++)
                                {
									for(unsigned int l_mut = 0; l_mut < 1; l_mut++)
									{
                                        xPool.enqueue(
												[&]
												(
													size_t tID,
													unsigned int i,
													unsigned int j,
													unsigned int k,
													unsigned int i_mut,
													unsigned int j_mut,
													unsigned int k_mut,
													unsigned int l_mut,
													unsigned int pos,
													std::string sequence
												)
                                                {
		{
			uint8_t q[sequence.size()];
			for(int x = sequence.size()-1; x >=0; x--)
			{
					if(sequence.at(x) == 'A')
							q[x] = 0;
					if(sequence.at(x) == 'C')
							q[x] = 1;
					if(sequence.at(x) == 'G')
							q[x] = 2;
					if(sequence.at(x) == 'T')
							q[x] = 3;
			}//for
			q[k] = (q[k] + k_mut) % 4;
			q[j] = (q[j] + j_mut) % 4;
			q[i] = (q[i] + i_mut) % 4;
			q[sequence.size()-3] = (q[sequence.size()-3] + l_mut) % 4;
			SA_IndexInterval ik(
					xIndex->L2[(int)q[sequence.size()-1]] + 1, 
					xIndex->L2[(int)q[sequence.size()-1] + 1] - xIndex->L2[(int)q[sequence.size()-1]]
			);
			
			
			for(int x = sequence.size()-2; x >=0; x--)
			{
				const uint8_t c = q[x]; // character at position i in the query
				
				// perform one step of the extension
				ik = bwt_extend_backward(ik, c, xIndex);

				if(ik.size() == 0)
					break;
			}//for
			
			if(ik.size() != 0)
			{
			
				std::string seq("");
				for(unsigned int x = 0; x < sequence.size(); x++)
				{
					if( (i == x && i_mut != 0) || (j == x && j_mut != 0) || (k == x && k_mut != 0))
					{
						if(q[x] == 0)
							seq.append("A");
						if(q[x] == 1)
							seq.append("C");
						if(q[x] == 2)
							seq.append("G");
						if(q[x] == 3)
							seq.append("T");
					}
					else
					{
						if(q[x] == 0)
							seq.append("a");
						if(q[x] == 1)
							seq.append("c");
						if(q[x] == 2)
							seq.append("g");
						if(q[x] == 3)
							seq.append("t");
					}
				}//for
				int missmatchAmount = 0;
				if(i_mut != 0)
					missmatchAmount++;
				if(j_mut != 0)
					missmatchAmount++;
				if(k_mut != 0)
					missmatchAmount++;
				std::list<int64_t> list = std::list<int64_t>();
				for(auto p = ik.start(); p < ik.end(); p++)
				{
					auto ulIndexOnRefSeq = xIndex->bwt_sa(p);
					ulIndexOnRefSeq = xIndex->getRefSeqLength() - (ulIndexOnRefSeq + sequence.size()) - 1;
					list.push_back(ulIndexOnRefSeq);
				}//for
				m.lock();
				//std::tuple<std::string, unsigned int, std::list<bwtint_t>>
				results[pos].push_back(std::make_tuple(seq, missmatchAmount, list));
				m.unlock();
			}
		}

		return;
		//there is no need to look for that since the fm index contains the reverse complement..
		//thus should be changed in order to half the ram requirements
		//===========reversed complement=============
		{
			uint8_t q[sequence.size()];
			for(int x = sequence.size()-1; x >=0; x--)
			{
					if(sequence.at(x) == 'T')
							q[x] = 0;
					if(sequence.at(x) == 'G')
							q[x] = 1;
					if(sequence.at(x) == 'C')
							q[x] = 2;
					if(sequence.at(x) == 'A')
							q[x] = 3;
			}//for
			q[k] = (q[k] + k_mut) % 4;
			q[j] = (q[j] + j_mut) % 4;
			q[i] = (q[i] + i_mut) % 4;
			SA_IndexInterval ik(
				xIndex->L2[(int)q[0]] + 1, 
				xIndex->L2[(int)q[0] + 1] - xIndex->L2[(int)q[0]]
			);
			
			
			for(unsigned int x = 1; x < sequence.size(); x++)
			{
				const uint8_t c = q[x]; // character at position i in the query
				
				// perform one step of the extension
				ik = bwt_extend_backward(ik, c, xIndex);

				if(ik.size() == 0)
					break;
			}//for
			
			if(ik.size() == 0)
				return;
			
			std::string seq("");
			for(int x = sequence.size()-1; x >=0; x--)
			{
				if( ((int)i == x && i_mut != 0) || ((int)j == x && j_mut != 0) || ((int)k == x && k_mut != 0))
				{
					if(q[x] == 0)
						seq.append("A");
					if(q[x] == 1)
						seq.append("C");
					if(q[x] == 2)
						seq.append("G");
					if(q[x] == 3)
						seq.append("T");
				}
				else
				{
					if(q[x] == 0)
						seq.append("a");
					if(q[x] == 1)
						seq.append("c");
					if(q[x] == 2)
						seq.append("g");
					if(q[x] == 3)
						seq.append("t");
				}
			}//for
			int missmatchAmount = 0;
			if(i_mut != 0)
				missmatchAmount++;
			if(j_mut != 0)
				missmatchAmount++;
			if(k_mut != 0)
				missmatchAmount++;
			std::list<int64_t> list = std::list<int64_t>();
			for(auto p = ik.start(); p < ik.end(); p++)
			{
				auto ulIndexOnRefSeq = xIndex->bwt_sa(p);
				ulIndexOnRefSeq = xIndex->getRefSeqLength() - (ulIndexOnRefSeq + sequence.size()) - 1;
				list.push_back(ulIndexOnRefSeq);
			}//for
			m.lock();
			//std::tuple<std::string, unsigned int, std::list<bwtint_t>>
			results[pos].push_back(std::make_tuple(seq, missmatchAmount, list));
			m.unlock();
		}
												},//lambda
												i, j, k, i_mut, j_mut, k_mut, l_mut, pos, sequence
											);
									}//for
								}//for
                                }//for
						}//for
                        }//for
                }//for
				}//for
			pos++;
        }//for
	}//scope for xPool
	
	for( auto& result : results)
		std::sort(result.begin(), result.end(), 
			[](std::tuple<std::string, unsigned int, std::list<int64_t>> a,
				std::tuple<std::string, unsigned int, std::list<int64_t>> b)
			{
				return std::get<1>(a) < std::get<1>(b);
			}
		);
		
	
	std::sort(results.begin(), results.end(), 
		[](std::vector<std::tuple<std::string, unsigned int, std::list<int64_t>>> a,
			std::vector<std::tuple<std::string, unsigned int, std::list<int64_t>>> b)
		{
			unsigned int tota[] = {0,0,0,0};
			unsigned int totb[] = {0,0,0,0};
			for(auto entry : a)
			{
				tota[std::get<1>(entry)] += std::get<2>(entry).size();
			}
			for(auto entry : b)
			{
				totb[std::get<1>(entry)] += std::get<2>(entry).size();
			}
			if(tota[2] == totb[2] && tota[1] == totb[1])
				return tota[3] < totb[3];
			if(tota[1] == totb[1])
				return tota[2] < totb[2];
			return tota[1] < totb[1];
		}
	);

	/*std::cout << "num mismatches\tsequence\tlocation" << std::endl;
	for( auto result : results)
	{
		for( auto entry : result)
		{
			if(std::get<1>(entry) < 3)
			{
				std::cout << std::get<1>(entry) << ": " << std::get<0>(entry) << ": ";
				for( auto blub : std::get<2>(entry))
				{
					std::cout << blub << " ";
				}
				std::cout << std::endl;
			}
		}
		std::cout << "---------------------------------------------------------------" << std::endl;
	}//for*/
	std::cout << "0\t1\t2\t3\tmismatch hits (W\\O NGG)" << std::endl;
	std::cout << "===================================================================" << std::endl;
	for( auto result : results)
	{
		unsigned int tot[] = {0,0,0,0};
		for(auto entry : result)
		{
			tot[std::get<1>(entry)] += std::get<2>(entry).size();
		}
		std::cout << tot[0] << "\t" << tot[1] << "\t" << tot[2] << "\t" << tot[3]
		<< "\t";
		for(auto entry : result)
		{
			if(std::get<1>(entry) == 0)
			{
				std::cout << std::get<0>(entry);
				break;
			}
		}
		std::cout << std::endl;
	}//for

	results.clear();
    {
        ThreadPool xPool(48);		

		std::mutex m;
		unsigned int pos = 0;

        for(std::string& sequence_short : searchFor)
        {
				std::string sequence = std::string(sequence_short).append("AGG");
				results.push_back(std::vector<std::tuple<std::string, unsigned int, std::list<int64_t>>>());
                for(unsigned int i = 0; i < sequence.size()-3; i++)
                {
                for(unsigned int i_mut = (i == 0 ? 0 : 1); i_mut < 4; i_mut++)
                {
                        for(unsigned int j = i+1; j < sequence.size()-3; j++)
                        {
                        for(unsigned int j_mut = (j == i+1 && i_mut == 0 ? 0 : 1); j_mut < 4; j_mut++)
                        {
                                for(unsigned int k = j+1; k < sequence.size()-3; k++)
                                {
                                for(unsigned int k_mut = (k == j+1 && j_mut == 0  ? 0 : 1); k_mut < 4; k_mut++)
                                {
									for(unsigned int l_mut = 0; l_mut < 4; l_mut++)
									{
                                        xPool.enqueue(
												[&]
												(
													size_t tID,
													unsigned int i,
													unsigned int j,
													unsigned int k,
													unsigned int i_mut,
													unsigned int j_mut,
													unsigned int k_mut,
													unsigned int l_mut,
													unsigned int pos,
													std::string sequence
												)
                                                {
		{
			uint8_t q[sequence.size()];
			for(int x = sequence.size()-1; x >=0; x--)
			{
					if(sequence.at(x) == 'A')
							q[x] = 0;
					if(sequence.at(x) == 'C')
							q[x] = 1;
					if(sequence.at(x) == 'G')
							q[x] = 2;
					if(sequence.at(x) == 'T')
							q[x] = 3;
			}//for
			q[k] = (q[k] + k_mut) % 4;
			q[j] = (q[j] + j_mut) % 4;
			q[i] = (q[i] + i_mut) % 4;
			q[sequence.size()-3] = (q[sequence.size()-3] + l_mut) % 4;
			SA_IndexInterval ik(
					xIndex->L2[(int)q[sequence.size()-1]] + 1, 
					xIndex->L2[(int)q[sequence.size()-1] + 1] - xIndex->L2[(int)q[sequence.size()-1]]
			);
			
			
			for(int x = sequence.size()-2; x >=0; x--)
			{
				const uint8_t c = q[x]; // character at position i in the query
				
				// perform one step of the extension
				ik = bwt_extend_backward(ik, c, xIndex);

				if(ik.size() == 0)
					break;
			}//for
			
			if(ik.size() != 0)
			{
			
				std::string seq("");
				for(unsigned int x = 0; x < sequence.size(); x++)
				{
					if( (i == x && i_mut != 0) || (j == x && j_mut != 0) || (k == x && k_mut != 0))
					{
						if(q[x] == 0)
							seq.append("A");
						if(q[x] == 1)
							seq.append("C");
						if(q[x] == 2)
							seq.append("G");
						if(q[x] == 3)
							seq.append("T");
					}
					else
					{
						if(q[x] == 0)
							seq.append("a");
						if(q[x] == 1)
							seq.append("c");
						if(q[x] == 2)
							seq.append("g");
						if(q[x] == 3)
							seq.append("t");
					}
				}//for
				int missmatchAmount = 0;
				if(i_mut != 0)
					missmatchAmount++;
				if(j_mut != 0)
					missmatchAmount++;
				if(k_mut != 0)
					missmatchAmount++;
				std::list<int64_t> list = std::list<int64_t>();
				for(auto p = ik.start(); p < ik.end(); p++)
				{
					auto ulIndexOnRefSeq = xIndex->bwt_sa(p);
					ulIndexOnRefSeq = xIndex->getRefSeqLength() - (ulIndexOnRefSeq + sequence.size()) - 1;
					list.push_back(ulIndexOnRefSeq);
				}//for
				m.lock();
				//std::tuple<std::string, unsigned int, std::list<bwtint_t>>
				results[pos].push_back(std::make_tuple(seq, missmatchAmount, list));
				m.unlock();
			}
		}

		return;
		//there is no need to look for that since the fm index contains the reverse complement..
		//thus should be changed in order to half the ram requirements
		//===========reversed complement=============
		{
			uint8_t q[sequence.size()];
			for(int x = sequence.size()-1; x >=0; x--)
			{
					if(sequence.at(x) == 'T')
							q[x] = 0;
					if(sequence.at(x) == 'G')
							q[x] = 1;
					if(sequence.at(x) == 'C')
							q[x] = 2;
					if(sequence.at(x) == 'A')
							q[x] = 3;
			}//for
			q[k] = (q[k] + k_mut) % 4;
			q[j] = (q[j] + j_mut) % 4;
			q[i] = (q[i] + i_mut) % 4;
			SA_IndexInterval ik(
				xIndex->L2[(int)q[0]] + 1, 
				xIndex->L2[(int)q[0] + 1] - xIndex->L2[(int)q[0]]
			);
			
			
			for(unsigned int x = 1; x < sequence.size(); x++)
			{
				const uint8_t c = q[x]; // character at position i in the query
				
				// perform one step of the extension
				ik = bwt_extend_backward(ik, c, xIndex);

				if(ik.size() == 0)
					break;
			}//for
			
			if(ik.size() == 0)
				return;
			
			std::string seq("");
			for(int x = sequence.size()-1; x >=0; x--)
			{
				if( ((int)i == x && i_mut != 0) || ((int)j == x && j_mut != 0) || ((int)k == x && k_mut != 0))
				{
					if(q[x] == 0)
						seq.append("A");
					if(q[x] == 1)
						seq.append("C");
					if(q[x] == 2)
						seq.append("G");
					if(q[x] == 3)
						seq.append("T");
				}
				else
				{
					if(q[x] == 0)
						seq.append("a");
					if(q[x] == 1)
						seq.append("c");
					if(q[x] == 2)
						seq.append("g");
					if(q[x] == 3)
						seq.append("t");
				}
			}//for
			int missmatchAmount = 0;
			if(i_mut != 0)
				missmatchAmount++;
			if(j_mut != 0)
				missmatchAmount++;
			if(k_mut != 0)
				missmatchAmount++;
			std::list<int64_t> list = std::list<int64_t>();
			for(auto p = ik.start(); p < ik.end(); p++)
			{
				auto ulIndexOnRefSeq = xIndex->bwt_sa(p);
				ulIndexOnRefSeq = xIndex->getRefSeqLength() - (ulIndexOnRefSeq + sequence.size()) - 1;
				list.push_back(ulIndexOnRefSeq);
			}//for
			m.lock();
			//std::tuple<std::string, unsigned int, std::list<bwtint_t>>
			results[pos].push_back(std::make_tuple(seq, missmatchAmount, list));
			m.unlock();
		}
												},//lambda
												i, j, k, i_mut, j_mut, k_mut, l_mut, pos, sequence
											);
									}//for
								}//for
                                }//for
						}//for
                        }//for
                }//for
				}//for
			pos++;
        }//for
	}//scope for xPool
	
	for( auto& result : results)
		std::sort(result.begin(), result.end(), 
			[](std::tuple<std::string, unsigned int, std::list<int64_t>> a,
				std::tuple<std::string, unsigned int, std::list<int64_t>> b)
			{
				return std::get<1>(a) < std::get<1>(b);
			}
		);
		
	
	std::sort(results.begin(), results.end(), 
		[](std::vector<std::tuple<std::string, unsigned int, std::list<int64_t>>> a,
			std::vector<std::tuple<std::string, unsigned int, std::list<int64_t>>> b)
		{
			unsigned int tota[] = {0,0,0,0};
			unsigned int totb[] = {0,0,0,0};
			for(auto entry : a)
			{
				tota[std::get<1>(entry)] += std::get<2>(entry).size();
			}
			for(auto entry : b)
			{
				totb[std::get<1>(entry)] += std::get<2>(entry).size();
			}
			if(tota[2] == totb[2] && tota[1] == totb[1])
				return tota[3] < totb[3];
			if(tota[1] == totb[1])
				return tota[2] < totb[2];
			return tota[1] < totb[1];
		}
	);

	/*std::cout << "num mismatches\tsequence\tlocation" << std::endl;
	for( auto result : results)
	{
		for( auto entry : result)
		{
			if(std::get<1>(entry) < 3)
			{
				std::cout << std::get<1>(entry) << ": " << std::get<0>(entry) << ": ";
				for( auto blub : std::get<2>(entry))
				{
					std::cout << blub << " ";
				}
				std::cout << std::endl;
			}
		}
		std::cout << "---------------------------------------------------------------" << std::endl;
	}//for*/
	std::cout << "0\t1\t2\t3\tmismatch hits (W NGG)" << std::endl;
	std::cout << "===================================================================" << std::endl;
	for( auto result : results)
	{
		unsigned int tot[] = {0,0,0,0};
		for(auto entry : result)
		{
			tot[std::get<1>(entry)] += std::get<2>(entry).size();
		}
		std::cout << tot[0] << "\t" << tot[1] << "\t" << tot[2] << "\t" << tot[3]
		<< "\t";
		for(auto entry : result)
		{
			if(std::get<1>(entry) == 0)
			{
				std::cout << std::get<0>(entry);
				break;
			}
		}
		std::cout << std::endl;
	}//for
}//main

void exportSegmentation()
{
	boost::python::def("analyse_crisper", analyseCRISPER);
	//export the segmentation class
	boost::python::class_<SegmentationContainer, boost::python::bases<Module>>(
			"Segmentation",
			"bBreakOnAmbiguousBase: weather the extension of "
			"intervals shall be stopped at N's\n",
			boost::python::init<boost::python::optional<bool>>(
				"arg1: self\n"
				"arg2: weather the extension of "
				"intervals shall be stopped at N's\n"
			)
		)
		.def_readwrite("bBreakOnAmbiguousBase", &SegmentationContainer::bBreakOnAmbiguousBase)
		;
	
}//function