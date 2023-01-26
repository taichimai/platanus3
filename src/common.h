#ifndef COMMON_H
#define COMMON_H

#include<iostream>
#include<unistd.h>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<array>
#include<vector>
#include<set>
#include<string>
#include<bitset>
#include<tuple>
#include<queue>
#include<cmath>
#include<omp.h>

typedef int Error;
typedef std::unordered_map<std::string,std::string> ReadSet;
typedef std::set<std::string>  KmerSet;
typedef int id_counter;

char bit_to_base[4]={'A','C','G','T'};
std::unordered_map<char,int> base_to_bit={{'A',0},{'C',1},{'G',2},{'T',3}};
std::unordered_map<char,char> trans_base={{'A','T'},{'C','G'},{'G','C'},{'T','A'}};

#endif