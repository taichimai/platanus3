#ifndef BITCALC_CPP
#define BITCALC_CPP
#include"common.h"



template<typename BITSET>
BITSET GetFirstKmerForward(std::string kmer){
    BITSET A(0),C(1),G(2),T(3);
    BITSET left_kmer(0);
    for (int i=0;i<kmer.size();i++){
        left_kmer<<=2;
        if (kmer[i]=='A') left_kmer=(left_kmer | A);
        else if (kmer[i]=='C') left_kmer=(left_kmer | C);
        else if (kmer[i]=='G') left_kmer=(left_kmer | G);
        else if  (kmer[i]=='T') left_kmer=(left_kmer | T);
    }
    return left_kmer;
}

template<typename BITSET>
BITSET GetFirstKmerBackward(std::string kmer){
    BITSET A(0),C(1),G(2),T(3);
    BITSET left_kmer(0);
    for (int i=kmer.size()-1;i>=0;i--){
        left_kmer<<=2;
        if (kmer[i]=='A') left_kmer=(left_kmer | T);
        else if (kmer[i]=='C') left_kmer=(left_kmer | G);
        else if (kmer[i]=='G') left_kmer=(left_kmer | C);
        else if (kmer[i]=='T') left_kmer=(left_kmer | A);
    }
    return left_kmer;
} 

template<typename BITSET>
BITSET GetComplementKmer(BITSET forward_kmer){
    BITSET complement_kmer(0);
    for (int i=0;i<(complement_kmer.size()/2);i++){
        complement_kmer<<=1;
        if (forward_kmer[i*2+1]==0) complement_kmer.set(0);
        complement_kmer<<=1;
        if (forward_kmer[i*2]==0) complement_kmer.set(0);
    }
    return complement_kmer;
}

template<typename BITSET>
BITSET CompareBit(BITSET Bit_Fw,BITSET Bit_Bw,int bitset_length){
    for (int i=0;i<bitset_length;i++){
        if (Bit_Fw[bitset_length-1-i]<Bit_Bw[bitset_length-1-i]) return Bit_Fw;
        else if (Bit_Fw[bitset_length-1-i]>Bit_Bw[bitset_length-1-i]) return Bit_Bw;
    }
    return Bit_Fw;
}

template<typename BITSET>
std::string GetStringKmer(BITSET &BitKmer){
    std::string restored_kmer;
    int base_id;
    for (int i=(BitKmer.size()/2)-1;i>=0;i--){
        base_id=2*BitKmer[i*2+1]+BitKmer[i*2];
        restored_kmer+=bit_to_base[base_id];
    }
    return restored_kmer;
}
#endif