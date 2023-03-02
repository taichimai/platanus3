#include"common.h"
#include"ShowInfo.cpp"
#include"Options.cpp"
#include"Load.cpp"
#include"bloomfilter.cpp"
#include"MakeBloomFilter.cpp"
#include"DeBruijnGraph.cpp"
#include"Assemble.cpp"
#include"Logging.cpp"

int main(int argc,char **argv){
    //get parameter and data
    Logging logging;
    Options options;
    options.Parse(argc,argv,logging);
    if (options.readfile_name==""){
        ShowUsage();
        return 0;
    }

    ReadFile input_reads(options);
    input_reads.LoadFile();
    options.EstimateBloomfilter(input_reads.all_bases,logging);
    logging.WriteLog("read file loaded");
    options.PrintParameters(logging);

    //assemble
    Assemble_k(input_reads,options,logging);
    logging.WriteLog("finish");
    return 0;
}



