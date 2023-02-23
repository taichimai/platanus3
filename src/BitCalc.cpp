#ifndef BITCALC_CPP
#define BITCALC_CPP
#include"common.h"


template<typename LARGE_BITSET>
LARGE_BITSET GetFirstKmerForward(std::string kmer){
    LARGE_BITSET A(0),C(1),G(2),T(3);
    LARGE_BITSET left_kmer(0);
    for (int i=0;i<kmer.size();i++){
        left_kmer<<=2;
        if (kmer[i]=='A') left_kmer=(left_kmer | A);
        else if (kmer[i]=='C') left_kmer=(left_kmer | C);
        else if (kmer[i]=='G') left_kmer=(left_kmer | G);
        else if  (kmer[i]=='T') left_kmer=(left_kmer | T);
    }
    return left_kmer;
}

template<typename LARGE_BITSET>
LARGE_BITSET GetFirstKmerBackward(std::string kmer){
    LARGE_BITSET A(0),C(1),G(2),T(3);
    LARGE_BITSET left_kmer(0);
    for (int i=kmer.size()-1;i>=0;i--){
        left_kmer<<=2;
        if (kmer[i]=='A') left_kmer=(left_kmer | T);
        else if (kmer[i]=='C') left_kmer=(left_kmer | G);
        else if (kmer[i]=='G') left_kmer=(left_kmer | C);
        else if (kmer[i]=='T') left_kmer=(left_kmer | A);
    }
    return left_kmer;
} 

template<typename LARGE_BITSET>
LARGE_BITSET GetComplementKmer(LARGE_BITSET forward_kmer){
    LARGE_BITSET complement_kmer(0);
    for (int i=0;i<(complement_kmer.size()/2);i++){
        complement_kmer<<=1;
        if (forward_kmer[i*2+1]==0) complement_kmer.set(0);
        complement_kmer<<=1;
        if (forward_kmer[i*2]==0) complement_kmer.set(0);
    }
    return complement_kmer;
}

template<typename LARGE_BITSET>
LARGE_BITSET CompareBit(LARGE_BITSET Bit_Fw,LARGE_BITSET Bit_Bw,int bitset_length){
    for (int i=0;i<bitset_length;i++){
        if (Bit_Fw[bitset_length-1-i]<Bit_Bw[bitset_length-1-i]) return Bit_Fw;
        else if (Bit_Fw[bitset_length-1-i]>Bit_Bw[bitset_length-1-i]) return Bit_Bw;
    }
    return Bit_Fw;
}

template<typename LARGE_BITSET>
std::string GetStringKmer(LARGE_BITSET &BitKmer){
    std::string restored_kmer;
    int base_id;
    for (int i=(BitKmer.size()/2)-1;i>=0;i--){
        base_id=2*BitKmer[i*2+1]+BitKmer[i*2];
        restored_kmer+=bit_to_base[base_id];
    }
    return restored_kmer;
}
#endif