#ifndef MAKEBLOOMFILTER_CPP
#define MAKEBLOOMFILTER_CPP

#include"common.h"
#include"BitCalc.cpp"
#include"bloomfilter.cpp"

template<typename BITLENGTH>
BF<BITLENGTH> MakeBF(ReadSet &RS,uint64_t filtersize ,uint8_t numhashes,int kmer_length){
    BF<BITLENGTH> Kmer_BF(filtersize,numhashes);

    BITLENGTH A_right(0); BITLENGTH A_left=(A_right<<(kmer_length*2-2)); 
    BITLENGTH C_right(1); BITLENGTH C_left=(C_right<<(kmer_length*2-2));
    BITLENGTH G_right(2); BITLENGTH G_left=(G_right<<(kmer_length*2-2));
    BITLENGTH T_right(3); BITLENGTH T_left=(T_right<<(kmer_length*2-2));
    std::vector<BITLENGTH> end_bases={A_left,C_left,G_left,T_left,A_right,C_right,G_right,T_right};
    
    #pragma omp parallel for  num_threads(4) 
    for(size_t b=0;b<RS.bucket_count();b++)
    for(auto bi=RS.begin(b);bi!=RS.end(b);bi++){
        std::string target_read = (bi->second);  
        BITLENGTH kmer_Fw=GetFirstKmerForward<BITLENGTH>(target_read.substr(0,kmer_length));
        BITLENGTH kmer_Bw=GetFirstKmerBackward<BITLENGTH>(target_read.substr(0,kmer_length));
        BITLENGTH KmerItem=CompareBit(kmer_Fw,kmer_Bw,kmer_length*2);       
        Kmer_BF.add(&KmerItem,kmer_length*2);
        for (int i=kmer_length;i<target_read.size();i++){
            kmer_Fw=((kmer_Fw<<2)| end_bases[ base_to_bit[ target_read[i] ] + 4 ] );
            kmer_Bw=((kmer_Bw>>2)| end_bases[ base_to_bit[ trans_base[target_read[i]] ] ] );
            KmerItem=CompareBit(kmer_Fw,kmer_Bw,kmer_length*2);
            Kmer_BF.add(&KmerItem,kmer_length*2);
        }
    }
    return Kmer_BF;
}
#endif

