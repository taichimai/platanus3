#ifndef BITCALC_CPP
#define BITCALC_CPP

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include<vector>
#include<set>
#include<string>
#include<array>
#include<bitset>

char bit_to_base[4]={'A','C','G','T'};

template<typename BITLENGTH>
BITLENGTH GetFirstKmerForward(std::string kmer){
    BITLENGTH A(0),C(1),G(2),T(3);
    BITLENGTH left_kmer(0);
    for (int i=0;i<(kmer.size()*2);i++){
        left_kmer<<=2;
        if (kmer[i]=='A') left_kmer=(left_kmer | A);
        else if (kmer[i]=='C') left_kmer=(left_kmer | C);
        else if (kmer[i]=='G') left_kmer=(left_kmer | G);
        else if  (kmer[i]=='T') left_kmer=(left_kmer | T);
    }
    return left_kmer;
}

template<typename BITLENGTH>
BITLENGTH GetFirstKmerBackward(std::string kmer){
    BITLENGTH A(0),C(1),G(2),T(3);
    BITLENGTH left_kmer(0);
    for (int i=(kmer.size()*2)-1;i>=0;i--){
        left_kmer<<=2;
        if (kmer[i]=='A') left_kmer=(left_kmer | T);
        else if (kmer[i]=='C') left_kmer=(left_kmer | G);
        else if (kmer[i]=='G') left_kmer=(left_kmer | C);
        else if (kmer[i]=='T') left_kmer=(left_kmer | A);
    }
    return left_kmer;
} 

template<typename BITLENGTH>
BITLENGTH GetComplementKmer(BITLENGTH forward_kmer){
    BITLENGTH complement_kmer(0);
    for (int i=0;i<(complement_kmer.size()/2);i++){
        complement_kmer<<=1;
        if (forward_kmer[i*2+1]==0) complement_kmer.set(0);
        complement_kmer<<=1;
        if (forward_kmer[i*2]==0) complement_kmer.set(0);
    }
    return complement_kmer;
}

template<typename BITLENGTH>
BITLENGTH CompareBit(BITLENGTH Bit_Fw,BITLENGTH Bit_Bw,int bitset_length){
    for (int i=0;i<bitset_length;i++){
        if (Bit_Fw[bitset_length-1-i]<Bit_Bw[bitset_length-1-i]) return Bit_Fw;
        else return Bit_Bw;
    }
    return Bit_Fw;
}


template<typename BITLENGTH>
std::string GetStringKmer(BITLENGTH &BitKmer){
    std::string restored_kmer;
    int base_id;
    for (int i=(BitKmer.size()/2)-1;i>=0;i--){
        base_id=2*BitKmer[i*2+1]+BitKmer[i*2];
        restored_kmer+=bit_to_base[base_id];
    }
    return restored_kmer;
}
#endif