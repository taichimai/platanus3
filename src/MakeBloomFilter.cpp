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

  uint64_t getmin(int a, int b, int k=0, int l=0, int r=-1) {
    if(r < 0) r = n;
    if(r <= a || b <= l) return INF;
    if(a <= l && r <= b) return node[k];

    uint64_t vl = getmin(a, b, 2*k+1, l, (l+r)/2);
    uint64_t vr = getmin(a, b, 2*k+2, (l+r)/2, r);
    return std::min(vl, vr);
  }

};

template<typename LARGE_BITSET>
BF<LARGE_BITSET> MakeBF(ReadSet &RS,KmerCount &KC,uint64_t filtersize ,uint8_t numhashes,int kmer_length,std::set<std::string> *seed_kmer,Logging &logging){
  BF<LARGE_BITSET> Kmer_BF(filtersize,numhashes);
  const uint32_t shortk_length=21;
  const uint64_t cov_threshold=2;

  LARGE_BITSET A_right(0); LARGE_BITSET A_left=(A_right<<(kmer_length*2-2)); 
  LARGE_BITSET C_right(1); LARGE_BITSET C_left=(C_right<<(kmer_length*2-2));
  LARGE_BITSET G_right(2); LARGE_BITSET G_left=(G_right<<(kmer_length*2-2));
  LARGE_BITSET T_right(3); LARGE_BITSET T_left=(T_right<<(kmer_length*2-2));
  std::vector<LARGE_BITSET> end_bases={A_left,C_left,G_left,T_left,A_right,C_right,G_right,T_right};

  std::bitset<42>  A_right_short(0); std::bitset<42> A_left_short=(A_right_short<<(shortk_length*2-2)); 
  std::bitset<42>  C_right_short(1); std::bitset<42> C_left_short=(C_right_short<<(shortk_length*2-2));
  std::bitset<42>  G_right_short(2); std::bitset<42> G_left_short=(G_right_short<<(shortk_length*2-2));
  std::bitset<42>  T_right_short(3); std::bitset<42> T_left_short=(T_right_short<<(shortk_length*2-2));
  std::vector<std::bitset<42>> end_bases_short={A_left_short,C_left_short,G_left_short,T_left_short,A_right_short,C_right_short,G_right_short,T_right_short};


  //for(size_t b=0;b<RS.bucket_count();b++)
  //for(auto bi=RS.begin(b);bi!=RS.end(b);bi++){

  
  logging.WriteLog(std::to_string(RS.size()));
  for(auto bi=RS.begin();bi!=RS.end();++bi){
    logging.WriteLog((bi->second)+"_S"); 
    //calculate 21-mer coverage
    std::vector<uint64_t> shortk_cov;
    shortk_cov.resize((bi->second).size()-shortk_length+1,0);

    std::bitset<42> shortk_for=GetFirstKmerForward<std::bitset<42>>((bi->second).substr(0,shortk_length));
    std::bitset<42> shortk_rev=GetFirstKmerBackward<std::bitset<42>>((bi->second).substr(0,shortk_length));


    for (int i=shortk_length-1;i<(bi->second).size();i++){
      logging.WriteLog(std::to_string(i));
      if (i!=shortk_length-1){
        shortk_for=( (shortk_for<<2) | end_bases_short[ base_to_bit[ (bi->second)[i] ] + 4 ]   );
        shortk_rev=( (shortk_rev>>2) | end_bases_short[ base_to_bit[ trans_base[(bi->second)[i]] ] ] );
      }
      shortk_cov[i-shortk_length+1]=KC[CompareBit(shortk_for,shortk_rev,42)];
    }
    
    logging.WriteLog((bi->second)); 

    //use segment tree for estimating large kmer coverage
    SegmentTree shortk_tree(shortk_cov);
    bool is_seedkmer_recorded=false;
    
    LARGE_BITSET kmer_Fw=GetFirstKmerForward<LARGE_BITSET>((bi->second).substr(0,kmer_length));
    LARGE_BITSET kmer_Bw=GetFirstKmerBackward<LARGE_BITSET>((bi->second).substr(0,kmer_length));
    LARGE_BITSET KmerItem;

    for (int i=kmer_length-1;i<(bi->second).size();i++){
      if (i!=kmer_length-1){
        kmer_Fw=((kmer_Fw<<2)| end_bases[ base_to_bit[ (bi->second)[i] ] + 4 ] );
        kmer_Bw=((kmer_Bw>>2)| end_bases[ base_to_bit[ trans_base[(bi->second)[i]] ] ] );
      }
      if (shortk_tree.getmin((i-kmer_length)+1,(i-kmer_length+1)+kmer_length-shortk_length+1)>=cov_threshold){
        KmerItem=CompareBit(kmer_Fw,kmer_Bw,kmer_length*2);
        Kmer_BF.add(&KmerItem,kmer_length*2);

        if (!is_seedkmer_recorded){
          std::string recording_kmer=GetStringKmer<LARGE_BITSET>(kmer_Fw);
          (*seed_kmer).insert(recording_kmer);
          is_seedkmer_recorded=true;
        }
      }
    }
    logging.WriteLog((bi->second)+"_F"); 
  }

  return Kmer_BF;
}
#endif

