#include"MakeBloomFilter.h"

template<typename BITLENGTH>
BF<BITLENGTH> MakeBF(ReadSet *RS,uint64_t filtersize ,uint8_t numhashes,int kmer_length){
    BF<BITLENGTH> Kmer_BF(filtersize,numhashes);

    BITLENGTH mp_Fw[100];
    BITLENGTH A_Fw(0);
    BITLENGTH C_Fw(1); 
    BITLENGTH G_Fw(2);
    BITLENGTH T_Fw(3);
    mp_Fw['A']=A_Fw; mp_Fw['C']=C_Fw; mp_Fw['G']=G_Fw; mp_Fw['T']=T_Fw;

    BITLENGTH mp_Bw[100];
    BITLENGTH A_Bw(3); A_Bw<<=(kmer_length*2-2);
    BITLENGTH C_Bw(2); C_Bw<<=(kmer_length*2-2);
    BITLENGTH G_Bw(1); G_Bw<<=(kmer_length*2-2);
    BITLENGTH T_Bw(0); T_Bw<<=(kmer_length*2-2);
    mp_Bw['A']=A_Bw; mp_Bw['C']=C_Bw; mp_Bw['G']=G_Bw; mp_Bw['T']=T_Bw;
    
    #pragma omp parallel for  num_threads(4) 
    for(size_t b=0;b<(*RS).bucket_count();b++)
    for(auto bi=(*RS).begin(b);bi!=(*RS).end(b);bi++){
        std::string target_read = (bi->second);  
        BITLENGTH kmer_Fw=GetFirstKmerForward<BITLENGTH>(target_read.substr(0,kmer_length));
        BITLENGTH kmer_Bw=GetFirstKmerBackward<BITLENGTH>(target_read.substr(0,kmer_length));
        BITLENGTH KmerItem=CompareBit(kmer_Fw,kmer_Bw,kmer_length*2);       
       
        Kmer_BF.add(&KmerItem,kmer_length*2);
        for (int i=kmer_length;i<target_read.size();i++){
            kmer_Fw=((kmer_Fw<<2)| mp_Fw[target_read[i]]);
            kmer_Bw=((kmer_Bw>>2)| mp_Bw[target_read[i]]);

            KmerItem=CompareBit(kmer_Fw,kmer_Bw,kmer_length*2);
            Kmer_BF.add(&KmerItem,kmer_length*2);
        }
    }
    return Kmer_BF;
}
