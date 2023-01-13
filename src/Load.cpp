#ifndef LOAD_CPP
#define LOAD_CPP
#include"common.h"

class FastaFile{
    public:
        ReadSet Reads;
        FastaFile(std::string FileName);
        KmerSet GetSeedKmer(int kmer_length);
};

FastaFile::FastaFile(std::string FileName){
    ReadSet  Input_Reads;
    std::ifstream lines(FileName);
    std::string line;
    std::string seq="";
    std::string readname="";
    while (getline(lines,line)) {
        if (line[0]=='>'){
            if (readname!=""){
                Input_Reads[readname]=seq;
                seq="";
            }
            readname=line;
        }
        else{
            seq=seq+line;
        }
    }
    Input_Reads[readname]=seq;
    this->Reads=Input_Reads;
}

KmerSet FastaFile::GetSeedKmer(int kmer_length){
    KmerSet EdgeKmers;
    for (auto itr=Reads.begin();itr!=Reads.end();++itr){
        EdgeKmers.insert((itr->second).substr(0,kmer_length));
    }
    return EdgeKmers;
}
#endif 