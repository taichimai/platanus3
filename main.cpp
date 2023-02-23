#include"common.h"
#include"Options.cpp"
#include"Load.cpp"
#include"bloomfilter.cpp"
#include"MakeBloomFilter.cpp"
#include"DeBruijnGraph.cpp"
#include"Assemble.cpp"


int main(int argc,char **argv){
    //get parameter and data
    Options parameters;
    parameters.Parse(argc,argv);
    ReadFile input_reads(parameters);
    input_reads.LoadFile();
    parameters.EstimateBloomfilter(input_reads.all_bases);
    std::cerr<<"read file loaded"<<"\n";
    parameters.PrintParameters();

    //assemble
    Assemble_k(input_reads,parameters);
    std::cerr<<"finish"<<"\n";
    return 0;
}



