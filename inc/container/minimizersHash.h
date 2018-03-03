/** 
 * @file alignment.h
 * @brief Implements a Container that holds a finished alignment.
 * @author Markus Schmidt
 */

#ifndef MINIMIZER_HASH_H
#define MINIMIZER_HASH_H

#include "container/segment.h"
#include "container/seed.h"
#include "util/support.h"
#include <map>
#include <fstream>
#include <queue>
#include <vector>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#define complement(x) (uint8_t)NucSeq::nucleotideComplement(x)

namespace libMA
{

    /*std::function<int(uint8_t a, uint8_t b, unsigned int iPos)>*/ 
    template<size_t k>
    class Minimizer
    {
    public:
        /**
         * @brief The minimizer sequence
         */
        uint8_t xSeq[k];

        int ordering (uint8_t a, uint8_t b, unsigned int iPos) const
        {
            /*
            * ordering in nucSeq is 'A', 'C', 'G', 'T', N
            * we want to have C, A, T, G, N however therefore translation table:
            */
            static const uint8_t translate[5] = {1, 0, 3, 2, 4};

            /*
            * In DNA sequences, the letters C and G often occur less frequently than
            * A and T.We assign the values 0, 1, 2, 3 to C, A, T, G, respectively,
            * for the odd numbered bases of k-mers, and reverse the
            * ordering for even numbered bases. This tends to start minimizers
            * with the valuable (in the sense of the significance
            * of a match) letters C and G, and makes the minimum k-mer
            * CGCGCG.... There are many other possibilities.
            */
            if(a == b)
                return 0;
            if(iPos % 2 == 0)
                return translate[b] < translate[a] ? 1 : -1;
            else
                return translate[a] < translate[b] ? 1 : -1;
        }//function

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        int compare(const uint8_t* pOtherSeq) const
        {
            unsigned int uiPos = 0;
            while(uiPos < k)
            {
                int res = ordering(xSeq[uiPos], pOtherSeq[uiPos], uiPos);
                if(res != 0)
                    return res;
                uiPos++;
            }//while
            return 0;
        }//function

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        bool operator<(const Minimizer& rOther) const
        {
            return compare(rOther.xSeq) < 0;
        }//operator

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        bool operator<=(const Minimizer& rOther) const
        {
            return compare(rOther.xSeq) <= 0;
        }//operator

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        bool operator>(const Minimizer& rOther) const
        {
            return compare(rOther.xSeq) > 0;
        }//operator

        void operator=(const uint8_t* pOtherSeq)
        {
            for(unsigned int i = 0; i < k; i++)
                xSeq[i] = pOtherSeq[i];
        }//operator

        void operator=(const Minimizer& rOther)
        {
            operator=(rOther.xSeq);
        }//operator

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        bool operator==(const Minimizer& rOther) const
        {
            return compare(rOther.xSeq) == 0;
        }//operator

        /**
         * @brief Compares two minimizers lexicographically according to the alphabet order given by
         * ordering()
         */
        bool operator!=(const Minimizer& rOther) const
        {
            return compare(rOther.xSeq) != 0;
        }//operator

        Minimizer(std::shared_ptr<NucSeq> pInit, nucSeqIndex uiStart)
        {
            assert(pInit->length() > k + uiStart);
            //this makes sure that a sequence and it's reverse complement have the same minimizer
            uint8_t xSeq2[k];
            for(int i = 0; i < (int)k; i++)
            {
                xSeq[i] = (*pInit)[i + uiStart];
                assert( ((int)k)-(i+1) >= 0);
                assert( ((int)k)-(i+1) < (int)k);
                xSeq2[k-(i+1)] = (uint8_t)complement((*pInit)[i + uiStart]);
            }//for

            //if the reverse complement minimizer is smaller switch to that
            if(compare(xSeq2) > 0)
                *this = xSeq2;
        }//constructor

        Minimizer(std::shared_ptr<NucSeq> pInit)
                :
            Minimizer(pInit, 0)
        {}//constructor

        Minimizer()
        {}//default constructor

    private:
        friend class boost::serialization::access;
        /**
         * This makes the class serializable
         */
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            for(size_t i = 0; i < k; i++)
                ar & xSeq[i];
        }
    };//class

    template<size_t w, size_t k>
    class MinimizersHash;

    template<size_t w, size_t k>
    class MinimizersVector 
            :
        public std::vector<std::pair<Minimizer<k>, nucSeqIndex>>,
        public Container
    {
    public:

        MinimizersVector()
        {}//constructor

        /**
         * creates a hash table from the vector
         */
        std::shared_ptr<MinimizersHash<w,k>> EXPORTED toHash();

        //overload
        bool canCast(std::shared_ptr<Container> c) const
        {
            return std::dynamic_pointer_cast<MinimizersVector>(c) != nullptr;
        }//function

        //overload
        std::string getTypeName() const
        {
            return "MinimizersVector";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new MinimizersVector());
        }//function
    };//class

    template<size_t w, size_t k>
    class MinimizersHash : public Container
    {
    public:
        //it might be better to actually use a hash table here...
        std::map<Minimizer<k>, std::vector<nucSeqIndex>> xHashMap;


        MinimizersHash()
        {}//constructor


        //overload
        bool canCast(std::shared_ptr<Container> c) const
        {
            return std::dynamic_pointer_cast<MinimizersHash>(c) != nullptr;
        }//function

        //overload
        std::string getTypeName() const
        {
            return "MinimizersHash";
        }//function

        //overload
        std::shared_ptr<Container> getType() const
        {
            return std::shared_ptr<Container>(new MinimizersHash());
        }//function
    private:
        friend class boost::serialization::access;
        /**
         * This makes the class serializable
         */
        template<class Archive>
        void save(Archive & ar, const unsigned int version) const
        {
            ar & k;
            ar & w;
            ar & xHashMap.size();
            for(auto& rContent : xHashMap)
            {
                ar & rContent.first;
                ar & rContent.second.size();
                for(auto& rElement : rContent.second)
                    ar & rElement;
            }//for
        }//function
        template<class Archive>
        void load(Archive & ar, const unsigned int version)
        {
            size_t check;
            ar & check;
            if(k != check)
                throw "@todo proper exception";
            ar & check;
            if(w != check)
                throw "@todo proper exception";
            size_t size = 0;
            ar & size;
            //variables to load into
            Minimizer<k> a;
            size_t s;
            nucSeqIndex b;
            while(size-- > 0)
            {
                ar & a;
                ar & s;
                auto pVec = &xHashMap[a];
                while(s-- > 0)
                {
                    ar & b;
                    pVec->push_back(b);
                }//while
            }//while
        }//function
        BOOST_SERIALIZATION_SPLIT_MEMBER()
    public:
        static std::shared_ptr<MinimizersHash> fromFile(std::string sFileName)
        {
            auto pRet = std::make_shared<MinimizersHash>();
            // create and open an archive for input
            std::ifstream ifs(sFileName.c_str());
            boost::archive::text_iarchive ia(ifs);
            // read class state from archive
            ia >> *pRet;
            // archive and stream closed when destructors are called
            return pRet;
        }

        void toFile(std::string sFileName)
        {
            std::ofstream ofs(sFileName.c_str());
            boost::archive::text_oarchive oa(ofs);
            // write class instance to archive
            oa << *this;
        }

        std::vector<nucSeqIndex>& operator[](Minimizer<k>& rKey)
        {
            return xHashMap[rKey];
        }//function

        class iterator : public Container
        {
        public:
            typedef std::tuple< //content
                    Minimizer<k>, //minimizer sequence
                    std::vector<nucSeqIndex>::iterator, // reference position
                    nucSeqIndex, // query position
                    std::vector<nucSeqIndex>::iterator // end iterator
                > queueContent;
            typedef std::function<bool(queueContent,queueContent)> queueCompFunc;
            std::priority_queue<queueContent, std::vector<queueContent>, queueCompFunc> xQueue;
            iterator(
                    std::shared_ptr<MinimizersVector<w,k>> pQueryMinimizers, 
                    MinimizersHash& rHash
                )
                    :
                xQueue(
                    []
                    (
                        std::tuple<
                            Minimizer<k>, 
                            std::vector<nucSeqIndex>::iterator, 
                            nucSeqIndex, 
                            std::vector<nucSeqIndex>::iterator
                        > a,
                        std::tuple<
                            Minimizer<k>, 
                            std::vector<nucSeqIndex>::iterator, 
                            nucSeqIndex, 
                            std::vector<nucSeqIndex>::iterator
                        > b
                    )
                    {
                        if(*std::get<1>(a) == *std::get<1>(b))
                            //sort so that larger query positions come first
                            return std::get<2>(a) > std::get<2>(b);
                        //sort so that smaller reference positions come first (with priority)
                        return *std::get<1>(a) < *std::get<1>(b);
                    }//lambda
                )//constructor for priority queue
            {
                auto xIter = pQueryMinimizers->begin();
                while(xIter != pQueryMinimizers->end())
                {
                    xQueue.push(std::make_tuple(
                        xIter->first, //minimizer sequence
                        rHash[xIter->first].begin(), //vector of reference positions
                        xIter->second, // query position
                        rHash[xIter->first].end() // end iterator
                    ));
                }//while
            }//constructor

            iterator()
                    :
                xQueue()
            {}//default constructor
            
            //overload
            bool canCast(std::shared_ptr<Container> c) const
            {
                return std::dynamic_pointer_cast<iterator>(c) != nullptr;
            }//function

            //overload
            std::string getTypeName() const
            {
                return "MinimizersHash::iterator";
            }//function

            //overload
            std::shared_ptr<Container> getType() const
            {
                return std::shared_ptr<Container>(new iterator());
            }//function

            bool empty()
            {
                return xQueue.empty();
            }//function

            Seed operator*() const
            {
                assert(!xQueue.empty());
                return Seed(std::get<2>(xQueue.top()), k, *std::get<1>(xQueue.top()));
            }//operator

            iterator& operator++()
            {
                //copy the top element
                queueContent xTuple;
                xTuple = xQueue.top();
                //increment the vector iterator of the top element
                std::get<1>(xTuple)++;
                //remove the top element from the priority queue
                xQueue.pop();
                //if there are still elements left in the reference vector
                //readd the previous top element
                if(std::get<1>(xTuple) != std::get<3>(xTuple))
                    xQueue.push(xTuple);
                return *this;
            }//operator
        };//class

        iterator begin(std::shared_ptr<MinimizersVector<w,k>> pQueryMinimizers)
        {
            return iterator(pQueryMinimizers, *this);
        }//function

        std::shared_ptr<Seeds> toSeeds(std::shared_ptr<MinimizersVector<w,k>> pQueryMinimizers)
        {
            auto it = begin(pQueryMinimizers);
            auto pRet = std::make_shared<Seeds>();
            while(!it.empty())
            {
                pRet->push_back(*it);
                ++it;
            }//while
            return pRet;
        }//function

    };//inner class

    template<size_t w, size_t k>
    std::shared_ptr<MinimizersHash<w,k>> MinimizersVector<w,k>::toHash()
    {
        assert(!this->empty());
        //sort now so that there is no need to sort for the SOC
        //sideeffect: lets us fill the hashtable much faster and easier
        std::cout << "sorting..." << std::endl;
        std::sort(
            this->begin(), 
            this->end(),
            []
            (std::pair<Minimizer<k>, nucSeqIndex>& a, 
                std::pair<Minimizer<k>, nucSeqIndex>& b)
            {
                if(a.first == b.first)
                    //the soc order
                    return a.second < b.second;
                //the fill order (this takes priority)
                return a.first < b.first;
            }//lambda
        );//sort function call
        std::cout << "done sorting" << std::endl;
        auto pRet = std::make_shared<MinimizersHash<w,k>>();
        //remember the vector in the hash table that we are currently filling
        auto pCurrent = this->begin();
        //remember the last element we inserted so that we know when to change to the next vec
        auto xLast = pCurrent->first;
        //this is in order to check for duplicates
        nucSeqIndex uiLast = (nucSeqIndex)-1;// set to max value...
        //fill in the first element
        std::vector<nucSeqIndex>* pAppend = &(*pRet)[xLast];
        unsigned int i = 0;
        //fill in all other elements
        while(++pCurrent != this->end())
        {
            if(i++ % 1000000 == 0)
                std::cout << i << "/" << this->size() << std::endl;
            //check that there are no duplicates
            assert(uiLast != pCurrent->second);
            //in this case we need to create a new vector in the hash table
            if(pCurrent->first != xLast)
            {
                xLast = pCurrent->first;
                pAppend = &(*pRet)[xLast];
            }//if
            //save the current element
            pAppend->push_back(pCurrent->second);
            uiLast = pCurrent->second;
        }//while
        return pRet;
    }//function

}//namespace libMA

/**
 * @brief Exposes the MinimizersHash container to boost python.
 * @ingroup export
 */
void exportMinimizersHash();

#endif