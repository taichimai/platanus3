#include"common.h"
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
    uint64_t FilterSize=10000000;
    uint8_t NumHashes=5;
    const uint32_t bitset_length=10; // kmer_length*2
    uint32_t  kmer_length=bitset_length/2;
    //seed k-mer
    KmerSet seed_kmer=inputreads.GetSeedKmer(kmer_length);
    //make bloomfilter
    BF<std::bitset<bitset_length> > first_bloom_filter=MakeBF<std::bitset<bitset_length> >(&inputreads.Reads,FilterSize,NumHashes,kmer_length);
    //make debruijngraph
    DeBruijnGraph<std::bitset<bitset_length> > first_dbg(kmer_length);
    first_dbg.MakeDBG(seed_kmer,first_bloom_filter,FilterSize,NumHashes);
    first_dbg.PrintGraph();

    return 0;
}



