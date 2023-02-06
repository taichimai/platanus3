#ifndef ASSEMBLE_CPP
#define ASSEMBLE_CPP
#include"common.h"

template<typename BITLENGTH>
int Assemble (ReadFile &input_reads,Options &parameters){
    KmerSet seed_kmer=input_reads.GetSeedKmer(parameters.kmer_length);
    BF<BITLENGTH> first_bloom_filter=MakeBF<BITLENGTH>(input_reads.reads,parameters.filter_size,parameters.num_hashes,parameters.kmer_length); //make bloomfilter
    std::cerr<<"bloom filter loaded"<<std::endl;
    DeBruijnGraph<BITLENGTH > first_dbg(parameters.kmer_length,first_bloom_filter);   //make debruijngraph
    std::cerr<<"make de bruijn graph"<<std::endl;
    first_dbg.MakeDBG(seed_kmer,parameters.filter_size,parameters.num_hashes);
    std::cerr<<"de bruijn graph loaded"<<std::endl;
    first_dbg.CountNodeCoverage(input_reads.reads);
    first_dbg.PrintGraph();

    return 0;
}

int Assemble_k(ReadFile &input_reads,Options &parameters){
    switch (parameters.kmer_length){
        case 5:
            return Assemble<std::bitset<10>>(input_reads,parameters);
        case 21:
            return Assemble<std::bitset<42>>(input_reads,parameters);
        case 63:
            return Assemble<std::bitset<126>>(input_reads,parameters);
        case 501:
            return Assemble<std::bitset<1002>>(input_reads,parameters);
        case 1001:
            return Assemble<std::bitset<2002>>(input_reads,parameters);
        
        default:
            return 1;
    }
}
#endif
