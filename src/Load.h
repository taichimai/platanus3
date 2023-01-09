#ifndef LOAD_H
#define LOAD_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<vector>
#include<set>
#include<string>


typedef std::unordered_map<std::string,std::string> ReadSet;
typedef std::set<std::string>  KmerSet;

class FastaFile{
    public:
        ReadSet Reads;
        FastaFile(std::string FileName);
        KmerSet GetSeedKmer(int kmer_length);
};
#endif 