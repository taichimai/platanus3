#ifndef ASSEMBLE_CPP
#define ASSEMBLE_CPP
#include"common.h"


template<typename LARGE_BITSET>
int Assemble (ReadFile &input_reads,Options &options,Logging &logging){
    logging.WriteLog("Assemble");
    input_reads.CountShortKmer(options.shortk_length);
    logging.WriteLog("count short kmer");
    std::set<std::string> seedkmer;
    BF<LARGE_BITSET> first_bloom_filter=MakeBF<LARGE_BITSET>(input_reads.reads,input_reads.shortk_database,options.filter_size,options.num_hashes,options.kmer_length,&seedkmer); //make bloomfilter

    logging.WriteLog("bloom filter loaded");
    logging.WriteLog("get seed kmer");
    logging.WriteLog("seed kmer num= "+std::to_string(seedkmer.size()));

    //make debruijngraph
    DeBruijnGraph<LARGE_BITSET> first_dbg(options.kmer_length,first_bloom_filter);
    logging.WriteLog("start graph extention");
    first_dbg.MakeDBG(seedkmer,options.filter_size,options.num_hashes,options.threads_num,logging);

    logging.WriteLog("de bruijn graph loaded");
    first_dbg.CountNodeCoverage(input_reads.reads);
    logging.WriteLog("count node coverage");
    first_dbg.PrintGraph();
    return 0;
}

int Assemble_k(ReadFile &input_reads,Options &options,Logging &logging){
    switch (options.kmer_length){
        case 5:
            return Assemble<std::bitset<10>>(input_reads,options,logging);
        case 21:
            return Assemble<std::bitset<42>>(input_reads,options,logging);
        case 25:
            return Assemble<std::bitset<50>>(input_reads,options,logging);
        case 63:
            return Assemble<std::bitset<126>>(input_reads,options,logging);
        case 101:
            return Assemble<std::bitset<202>>(input_reads,options,logging);
        case 501:
            return Assemble<std::bitset<1002>>(input_reads,options,logging);
        case 1001:
            return Assemble<std::bitset<2002>>(input_reads,options,logging);
        case 2001:
            return Assemble<std::bitset<4002>>(input_reads,options,logging);
        case 3001:
            return Assemble<std::bitset<6002>>(input_reads,options,logging);
        
        default:
            return 1;
    }
}
#endif
