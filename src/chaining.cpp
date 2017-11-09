#include "chaining.h"

std::vector<ContainerType> Chaining::getInputType()
{
	return std::vector<ContainerType>{
			//the strip of consideration
			ContainerType::seeds,
		};
}//function

ContainerType Chaining::getOutputType()
{
	return ContainerType::seeds;
}//function


std::shared_ptr<Container> Chaining::execute(
        std::vector<std::shared_ptr<Container>> vpInput
    )
{
    std::shared_ptr<Seeds> pSeeds = std::static_pointer_cast<Seeds>(vpInput[0]);

    std::vector<RMQ<int64_t>::RMQData> data1;
    std::vector<RMQ<int64_t>::RMQData> data2;
    data1.push_back(RMQ<int64_t>::RMQData(-100000,0,nullptr, -1000000));
    data2.push_back(RMQ<int64_t>::RMQData(0,-100000,nullptr, -1000000));

    std::vector<std::shared_ptr<Chain>> chains;
    for(Seed seed : *pSeeds)
    {
        std::shared_ptr<Chain> chain = std::shared_ptr<Chain>(new Chain(seed));
        chains.push_back(chain);
        data1.push_back(t1(seed, chain));
        data2.push_back(t2(seed, chain));
    }//for

    //first octant
    RMQ<int64_t> d1 = RMQ<int64_t>(data1);
    //2nd octant
    RMQ<int64_t> d2 = RMQ<int64_t>(data2);

    //it's importent to do this after the initialization of the trees
    //sonce the trees sort the vector datastructure...
    for(unsigned int i = 0; i < data1.size(); i++)
        if(data1[i].pChain != nullptr)
            data1[i].pChain->t1 = &data1[i];
    for(unsigned int i = 0; i < data2.size(); i++)
        if(data2[i].pChain != nullptr)
            data2[i].pChain->t2 = &data2[i];
    std::shared_ptr<Chain> bestChain = chains[0];



    std::sort(
        chains.begin(), chains.end(),
        []
        (const std::shared_ptr<Chain> a, const std::shared_ptr<Chain> b)
        {
            if(a->s.end_ref() == b->s.end_ref())
                return a->s.end() < b->s.end();
            return a->s.end_ref() < b->s.end_ref();
        }//lambda
    );//function call

    //the actual chaining
    for(std::shared_ptr<Chain> chain : chains)
    {
        DEBUG_2(
            std::cout << "computing chain for (ref/query/size): " << chain->s.end_ref() << " "
                    << chain->s.end() << " "
                    << chain->s.size() << std::endl;
        )
    /*
    * SWITCH between allowing overlaps of chains our not
    * true = no overlaps
    * WARNING: chaining code is not 100% correct when allowing overlaps
    * (basically each certain match will be scored as possible match only in that case)
    */ 
    #define STARTS true
    #if STARTS
        RMQ<int64_t>::RMQData& a = d1.rmq(
                -1000000,-1,
                (int64_t)chain->s.start_ref()-(int64_t)chain->s.start() - 1, chain->s.start() - 1//TODO: replace with function
            );
        RMQ<int64_t>::RMQData& b = d2.rmq(
                -1,(int64_t)chain->s.start()-(int64_t)chain->s.start_ref() - 1,
                chain->s.start_ref() - 1,-1000000 //TODO: replace with function
            );
    #else
        RMQ<int64_t>::RMQData& a = d1.rmq(
                -1000000,-1,
                (int64_t)chain->s.end_ref()-(int64_t)chain->s.end() - 1, chain->s.end() - 1//TODO: replace with function
            );
        RMQ<int64_t>::RMQData& b = d2.rmq(
                -1,(int64_t)chain->s.end()-(int64_t)chain->s.end_ref() - 1,
                chain->s.end_ref() - 1,-1000000 //TODO: replace with function
            );
    #endif
        DEBUG_2(
            std::cout << "score adjustments (second/first octant): -"
                    << gc1_start(chain->s) << " -"
                    << gc2_start(chain->s) << std::endl;
        )
        //using the RMQdata's to check for smaller insted of the chains
        //since we do not have to check for nullptrs this way
    #if STARTS
        if(a.score - gc1_start(chain->s) < b.score - gc2_start(chain->s))
    #else
        if(a.score - gc1_end(chain->s) < b.score - gc2_end(chain->s))
    #endif
            chain->pred = b.pChain;
        else
            chain->pred = a.pChain;

        if(chain->pred == nullptr)
            continue;
        assert(chain->pred != chain);
        DEBUG_2(
            std::cout << "candidate: " << chain->pred->s.end_ref() << " "
                    << chain->pred->s.end() << " "
                    << chain->pred->s.size() << std::endl;
        )
        assert(chain->pred->s.end_ref() < chain->s.end_ref());
        assert(chain->pred->s.end() < chain->s.end());
        //update score
        int64_t x = chain->s.end_ref() - chain->pred->s.end_ref();
        int64_t y = chain->s.end() - chain->pred->s.end();

        int64_t possibleMatches = std::min(x,y);
        int64_t certainMatches = std::min(possibleMatches, (int64_t)chain->s.size());
        possibleMatches -= certainMatches;
        int64_t insOrDls = std::max(x,y) - std::min(x,y);

        //set the score for the new chain
        int64_t addScore = (certainMatches * SCORE_MATCH
                - possibleMatches * COST_POSS_MATCH)
                - insOrDls * COST_INS_DEL;

        DEBUG_2(
            std::cout << "best candidate would require " << insOrDls << " insertions or deletions "
            << possibleMatches << " possible matches and " << certainMatches << " certain matches for a score of " << chain->pred->score << "; add score is: " << addScore << std::endl;
        )

        //if chaining this seed results in a worse chain then don't...
        //this makes our global chaining into a local chaining
        if(chain->score > addScore + chain->pred->score)
        {
            DEBUG(
                std::cout << "best chain not worth"<< std::endl;
            )
            chain->pred = nullptr;
        }//if
        else{
            chain->score = addScore + chain->pred->score;
            chain->t1->score = chain->score + gc1_end(chain->s);
            chain->t2->score = chain->score + gc2_end(chain->s);
        }//else


        DEBUG(
            std::cout << "current score: "<< chain->score << std::endl;
        )


        //remember best chain
        if(*bestChain < *chain)
            bestChain = chain;
    }//for

    std::shared_ptr<Seeds> pRet = std::shared_ptr<Seeds>(new Seeds());
//toggle for generating debug output
#if 1
    for(std::shared_ptr<Chain> chain : chains)
    {
        if(chain->pred != nullptr)
        {
            pRet->push_front(chain->s);
            pRet->push_front(chain->pred->s);
            pRet->push_front(Seed(0,0,0));
        }
    }//for
#endif
    while(bestChain != nullptr)
    {
        DEBUG_2(
            std::cout << "solution:" << std::endl;
            std::cout << bestChain->s.start() << " "
                    << bestChain->s.start_ref() << " "
                    << bestChain->s.size() << std::endl;
        )
        pRet->push_front(bestChain->s);
        bestChain = bestChain->pred;
    }//while
    DEBUG_2(
        std::cout << "done" << std::endl;
    )
    return pRet;
}//function

void exportChaining()
{
    //export the LineSweepContainer class
	boost::python::class_<Chaining, boost::python::bases<CppModule>>(
        "Chaining",
        "Uses chaining to remove contradicting "
        "matches within one strip of consideration.\n"
        "\n"
        "Execution:\n"
        "	Expects query, ref, strip_vec as input.\n"
        "		query: the query as NucleotideSequence\n"
        "		ref: the reference sequence as Pack\n"
        "		strip_vec: the areas that shall be evaluated as StripOfConsiderationVector\n"
        "	returns strip_vec.\n"
        "		strip_vec: the evaluated areas\n"
    );
}//function
