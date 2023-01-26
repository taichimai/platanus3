#ifndef OPTIONS_CPP
#define OPTIONS_CPP
#include<common.h>


class Options{
    public:
        std::string readfile_name;
        uint64_t filter_size=0;
        uint8_t num_hashes=10;
        uint32_t  kmer_length=5;
        double  error_rate=0.0005;
        

        void EstimateBloomfilter(uint64_t all_bases);
        bool Parse(int argc,char **argv);
        void PrintParameters();

};

bool Options::Parse(int argc,char **argv){
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
            default: 
                std::cerr << "Invalid option" << std::endl;
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
    filter_size=(-item_number*(std::log(false_positive_rate)))/std::pow(std::log(2),2);
    num_hashes=(std::log(2) * filter_size)/item_number;
    std::cerr<<"all_bases"<<" : "<<all_bases<<std::endl;
    std::cerr<<"item_number"<<" : "<<item_number<<std::endl;
    std::cerr<<"get filter_size"<<" : "<<filter_size<<std::endl;
}

void Options::PrintParameters(){
    std::cerr<<"readfile_name"<<" : "<<readfile_name<<std::endl;
    std::cerr<<"filter_size"<<" : "<<filter_size<<std::endl;
    std::cerr<<"num_hashes"<<" : "<<(int)num_hashes<<std::endl;
    std::cerr<<"kmer_length"<<" : "<<kmer_length<<std::endl;
    std::cerr<<"error_rate"<<" : "<<error_rate<<std::endl;
    return;
}


#endif

