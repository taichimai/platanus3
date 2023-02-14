#ifndef MAKEBLOOMFILTER_CPP
#define MAKEBLOOMFILTER_CPP

#include"common.h"
#include"BitCalc.cpp"
#include"bloomfilter.cpp"


struct SegmentTree {
private:
    uint64_t n;
    std::vector<uint64_t> node;

public:
    SegmentTree(std::vector<uint64_t> v) {
        int sz = v.size();
        n = 1; while(n < sz) n *= 2;
        node.resize(2*n-1, INF);
        for(int i=0; i<sz; i++) node[i+n-1] = v[i];
        for(int i=n-2; i>=0; i--) node[i] = std::min(node[2*i+1], node[2*i+2]);
    }
    int getmin(int a, int b, int k=0, int l=0, int r=-1) {
        if(r < 0) r = n;
        if(r <= a || b <= l) return INF;
        if(a <= l && r <= b) return node[k];

        int vl = getmin(a, b, 2*k+1, l, (l+r)/2);
        int vr = getmin(a, b, 2*k+2, (l+r)/2, r);
        return std::min(vl, vr);
    }
};

template<typename LARGE_BITSET>
BF<LARGE_BITSET> MakeBF(ReadSet &RS,KmerCount &KC,uint64_t filtersize ,uint8_t numhashes,int kmer_length){
    BF<LARGE_BITSET> Kmer_BF(filtersize,numhashes);

    LARGE_BITSET A_right(0); LARGE_BITSET A_left=(A_right<<(kmer_length*2-2)); 
    LARGE_BITSET C_right(1); LARGE_BITSET C_left=(C_right<<(kmer_length*2-2));
    LARGE_BITSET G_right(2); LARGE_BITSET G_left=(G_right<<(kmer_length*2-2));
    LARGE_BITSET T_right(3); LARGE_BITSET T_left=(T_right<<(kmer_length*2-2));
    std::vector<LARGE_BITSET> end_bases={A_left,C_left,G_left,T_left,A_right,C_right,G_right,T_right};

    const uint32_t shortk_length=21;
    const uint64_t cov_threshold=3;

    #pragma omp parallel for  num_threads(20) 
    for(size_t b=0;b<RS.bucket_count();b++)
    for(auto bi=RS.begin(b);bi!=RS.end(b);bi++){
        std::string target_read = (bi->second); 
        //calculate 21-mer coverage
        std::vector<uint64_t> shortk_cov;
        shortk_cov.resize(target_read.size()+1-shortk_length,0);
        for (int i = 0; i <target_read.size()+1-shortk_length; i++){   
            std::string shortk_for=target_read.substr(i,shortk_length);
            std::string shortk_rev=shortk_for;
            std::reverse(shortk_rev.begin(), shortk_rev.end()); // reverse
		    for (int j = 0; j < shortk_length; ++j) { // complement
			    if (shortk_rev[j] == 'A') shortk_rev[j] = 'T';
			    else if (shortk_rev[j] == 'T') shortk_rev[j] = 'A';
			    else if (shortk_rev[j] == 'C') shortk_rev[j] = 'G';
			    else if (shortk_rev[j] == 'G') shortk_rev[j] = 'C';
		    }
            if (shortk_for<shortk_rev){
                shortk_cov[i]=KC[shortk_for];
            }
            else{
                shortk_cov[i]=KC[shortk_rev];
            }
        }
        //use segment tree for estimating large kmer coverage
        SegmentTree shortk_tree(shortk_cov);

        LARGE_BITSET kmer_Fw=GetFirstKmerForward<LARGE_BITSET>(target_read.substr(0,kmer_length));
        LARGE_BITSET kmer_Bw=GetFirstKmerBackward<LARGE_BITSET>(target_read.substr(0,kmer_length));
        LARGE_BITSET KmerItem;
        if (shortk_tree.getmin(0,kmer_length-shortk_length+1)>=cov_threshold){
            KmerItem=CompareBit(kmer_Fw,kmer_Bw,kmer_length*2);
            Kmer_BF.add(&KmerItem,kmer_length*2);
        }
        for (int i=kmer_length;i<target_read.size();i++){
            kmer_Fw=((kmer_Fw<<2)| end_bases[ base_to_bit[ target_read[i] ] + 4 ] );
            kmer_Bw=((kmer_Bw>>2)| end_bases[ base_to_bit[ trans_base[target_read[i]] ] ] );
            KmerItem=CompareBit(kmer_Fw,kmer_Bw,kmer_length*2);
            if (shortk_tree.getmin((i-kmer_length)+1,(i-kmer_length+1)+kmer_length-shortk_length+1)>=cov_threshold){
                KmerItem=CompareBit(kmer_Fw,kmer_Bw,kmer_length*2);
                Kmer_BF.add(&KmerItem,kmer_length*2);
            } 
            Kmer_BF.add(&KmerItem,kmer_length*2);
            std::cerr<<"test"<<i<<"\n";
        }
    }
    return Kmer_BF;
}
#endif

