#ifndef DEBRUIJNGRAPH_H
#define DEBRUIJNGRAPH_H

#include<unordered_map>
#include<vector>
#include<set>
#include<string>
#include<bitset>
#include<queue>
#include<tuple>
#include"BitCalc.cpp"
#include"bloomfilter.cpp"


template<typename BITLENGTH>
class DeBruijnGraph{
    public:
      int bitset_length;
      std::vector<BITLENGTH> end_bases;
      std::vector< std::tuple<BITLENGTH,BITLENGTH> > junction_edges;  
      std::vector< std::tuple<BITLENGTH,std::string,std::string> > straight_edges;  
      std::unordered_map<std::string,std::string> straight_nodes;  
      std::unordered_map<BITLENGTH,std::string> junction_nodes;
      int straight_nodes_id;
      int junction_nodes_id;
      DeBruijnGraph(int k);
      void MakeDBG(KmerSet &seedkmer,BF<BITLENGTH> &all_kmers,uint64_t filtersize ,uint8_t numhashes);
      void SearchNode(BITLENGTH target_kmer,BF<BITLENGTH> *visited_kmers,std::queue<BITLENGTH> *visiting,BF<BITLENGTH> &all_kmers);
      BITLENGTH ExtendLeft(BITLENGTH &target_kmer,BITLENGTH &previous_kmer ,std::vector<char> *extend_bases,int previous_base,BF<BITLENGTH> &all_kmers,BF<BITLENGTH> *visited_kmers,std::queue<BITLENGTH> *visiting);
      BITLENGTH ExtendRight(BITLENGTH target_kmer,BITLENGTH &previous_kmer,std::vector<char> *extend_bases,int previous_base,BF<BITLENGTH> &all_kmers,BF<BITLENGTH> *visited_kmers,std::queue<BITLENGTH> *visiting);
      void AddJunctionNode(BITLENGTH &added_junction_node);
      void AddStraightNode(std::string &added_straight_node);
};
#endif