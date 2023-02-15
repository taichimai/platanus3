#ifndef ASSEMBLE_CPP
#define ASSEMBLE_CPP
#include"common.h"

template<typename LARGE_BITSET>
int Assemble (ReadFile &input_reads,Options &parameters){
    input_reads.CountShortKmer(parameters.shortk_length);
    std::cerr<<"count short kmer"<<"\n";
    std::set<LARGE_BITSET> seedkmer;
    BF<LARGE_BITSET> first_bloom_filter=MakeBF<LARGE_BITSET>(input_reads.reads,input_reads.shortk_database,parameters.filter_size,parameters.num_hashes,parameters.kmer_length,&seedkmer); //make bloomfilter
    std::cerr<<"bloom filter loaded"<<"\n";
    std::cerr<<"get seed kmer"<<"\n";
    //make debruijngraph
    DeBruijnGraph<LARGE_BITSET> first_dbg(parameters.kmer_length,first_bloom_filter);
    std::cerr<<"make de bruijn graph"<<"\n";
    first_dbg.MakeDBG(seedkmer,parameters.filter_size,parameters.num_hashes);
    std::cerr<<"de bruijn graph loaded"<<"\n";
    first_dbg.CountNodeCoverage(input_reads.reads);
    std::cerr<<"count node coverage"<<"\n";
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