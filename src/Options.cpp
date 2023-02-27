#ifndef OPTIONS_CPP
#define OPTIONS_CPP
#include"common.h"
#include"Logging.cpp"


class Options{
    public:
        std::string readfile_name;
        uint64_t filter_size=0;
        uint8_t num_hashes=10;
        int threads_num=8;
        uint32_t  kmer_length=5;
        uint32_t  shortk_length=21;
        double  error_rate=0.0005;
        

        void EstimateBloomfilter(uint64_t all_bases);
        bool Parse(int argc,char **argv,Logging &logging);
        void PrintParameters(Logging &logging);
};

bool Options::Parse(int argc,char **argv,Logging &logging){
    int opt;
    const char* options="i:m:k:";
    while((opt = getopt(argc, argv, options)) != -1){
        switch(opt)
        {
            case 'i':
                readfile_name = optarg;
                break;
            case 'm':
                filter_size = (uint64_t)std::stoll(optarg);
                break;
            case 'k':
                kmer_length =std::stoi(optarg);
                break;
            case 't':
                threads_num =std::stoi(optarg);
                break;
            default: 
                logging.WriteLog("Invalid option");
                return false;
                break;
        }
    }

    return true;
}

void Options::EstimateBloomfilter(uint64_t all_bases){
    if (filter_size!=0) return;
    double false_positive_rate=1.0e-6;
    uint64_t item_number=all_bases*error_rate*kmer_length;
    filter_size=(item_number*(-(std::log(false_positive_rate))))/std::pow(std::log(2),2);
    num_hashes=(std::log(2) * filter_size)/item_number;
    std::cerr<<"all_bases"<<" : "<<all_bases<<"\n";
    std::cerr<<"item_number"<<" : "<<item_number<<"\n";
    std::cerr<<"get filter_size"<<" : "<<filter_size<<"\n";
}
void Options::PrintParameters(Logging &logging){
    logging.WriteLog("readfile_name : "+readfile_name);
    logging.WriteLog("filter_size : "+std::to_string(filter_size));
    logging.WriteLog("num_hashes : "+std::to_string((int)num_hashes));
    logging.WriteLog("kmer_length : "+std::to_string(kmer_length));
    logging.WriteLog("error_rate : "+std::to_string(error_rate));
}


#endif

