#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<vector>
#include<set>
#include<string>
#include<bitset>
#include<tuple>
#include<queue>
#include"Load.cpp"
#include"bloomfilter.cpp"
#include"MakeBloomFilter.cpp"
#include"DeBruijnGraph.cpp"

int main(int argc,char **argv){

    //input
    std::string Input_FastaFile=argv[1];
    FastaFile inputreads(Input_FastaFile);
    std::cerr<<"fasta file loaded"<<std::endl;

    //parameter
    uint64_t FilterSize=1000000000;
    uint8_t NumHashes=5;
    uint32_t  kmer_length=501;
    const uint32_t bitset_length=1002; // kmer_length*2
    //seed k-mer
    KmerSet seed_kmer=inputreads.GetSeedKmer(kmer_length);
    
    //make bloomfilter
    BF<std::bitset<bitset_length> > first_bloom_filter=MakeBF<std::bitset<bitset_length> >(&inputreads.Reads,FilterSize,NumHashes,kmer_length);
    //make debruijngraph
    DeBruijnGraph<std::bitset<bitset_length> > first_dbg(kmer_length);
    first_dbg.MakeDBG(seed_kmer,first_bloom_filter,FilterSize,NumHashes);

    return 0;
}



