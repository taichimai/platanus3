#ifndef MAKEBLOOMFILTER_H
#define MAKEBLOOMFILTER_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<vector>
#include<set>
#include<string>
#include<array>
#include<bitset>
#include<omp.h>
#include"BitCalc.cpp"
#include"Load.h"
#include"bloomfilter.cpp"

template<typename BITLENGTH>
BF<BITLENGTH> MakeBF(ReadSet *RS,uint64_t filtersize ,uint8_t numhashes,int kmer_length);

#endif 
