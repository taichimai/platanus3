#ifndef DEBRUIJNGRAPH_CPP
#define DEBRUIJNGRAPH_CPP
#include"common.h"

template<typename BITLENGTH>
class DeBruijnGraph{
    public:
      //variables
      int kmer_length;
      int bitset_length;
      std::vector<BITLENGTH> end_bases;
      
      std::vector< std::tuple<BITLENGTH,char,BITLENGTH,char> > junction_edges;  
      std::vector< std::tuple<BITLENGTH,std::string,std::string> > straight_edges;  

      std::unordered_map<std::string,std::string> straight_nodes;  
      std::unordered_map<BITLENGTH,std::tuple<std::string,int> > junction_nodes;
      std::unordered_map<BITLENGTH,std::tuple<std::string,int> > joint_nodes;  

      id_counter straight_nodes_id=0;
      id_counter junction_nodes_id=0;
      id_counter joint_nodes_id=0;

      BF<BITLENGTH> *all_kmers;
      BF<BITLENGTH> visited_kmers;     //for making dbg
      std::queue<BITLENGTH> visiting;  //for making dbg

      //constructor
      DeBruijnGraph(int k,BF<BITLENGTH> &first_bloom_filter);

      //methods
      void MakeDBG(KmerSet &seedkmer,uint64_t filtersize ,uint8_t numhashes);
      void SearchNode(BITLENGTH target_kmer);

      BITLENGTH ExtendLeft(BITLENGTH &target_kmer,BITLENGTH &previous_kmer ,std::vector<char> *extend_bases,int previous_base);
      BITLENGTH ExtendRight(BITLENGTH target_kmer,BITLENGTH &previous_kmer,std::vector<char> *extend_bases,int previous_base);

      bool IsRecorded(BF<BITLENGTH> &bloomfilter,BITLENGTH &seqrching_kmer);
      void RecordKmer(BF<BITLENGTH> &bloomfilter,BITLENGTH &seqrching_kmer);
      void CheckDirections(std::vector<BITLENGTH> *stock_left,std::vector<BITLENGTH> *stock_right,BITLENGTH &target_kmer,int ignored_direction);

      void AddJunctionEdge(BITLENGTH &junction_node,BITLENGTH &nearby_node,std::string direction);
      void AddStraightEdge(BITLENGTH &joint_node,std::string &straight_node,std::string direction);

      void AddJunctionNode(BITLENGTH &added_junction_node);
      void AddJointNode(BITLENGTH &added_joint_node);
      void AddStraightNode(std::string &added_straight_node);

      void CountNodeCoverage(ReadSet &RS);
      void AddNodeCoverage(BITLENGTH &target_kmer);
      void PrintGraph();
      
};

template<typename BITLENGTH>
DeBruijnGraph<BITLENGTH>::DeBruijnGraph(int k,BF<BITLENGTH> &first_bloom_filter){
    all_kmers=&first_bloom_filter;
    kmer_length=k;
    bitset_length=k*2;
    BITLENGTH A_right(0); BITLENGTH A_left=(A_right<<(bitset_length-2)); 
    BITLENGTH C_right(1); BITLENGTH C_left=(C_right<<(bitset_length-2));
    BITLENGTH G_right(2); BITLENGTH G_left=(G_right<<(bitset_length-2));
    BITLENGTH T_right(3); BITLENGTH T_left=(T_right<<(bitset_length-2));
    end_bases={A_left,C_left,G_left,T_left,A_right,C_right,G_right,T_right};
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::MakeDBG(KmerSet &seedkmer,uint64_t filtersize ,uint8_t numhashes){

    visited_kmers.Set_BF(filtersize ,numhashes);

    for (auto itr=seedkmer.begin();itr!=seedkmer.end();++itr){
        BITLENGTH fw_kmer=GetFirstKmerForward<BITLENGTH>(*itr);
        if (IsRecorded(visited_kmers,fw_kmer)) continue;
        visiting.push(fw_kmer);

        while(!visiting.empty()){
            BITLENGTH visiting_kmer=visiting.front();
            visiting.pop();
            SearchNode(visiting_kmer);
        } 
    }
    visited_kmers.Delete_Filter();
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::SearchNode(BITLENGTH target_kmer){

    if (IsRecorded(visited_kmers,target_kmer)) return;

    //search eight directions
    std::vector<BITLENGTH> stock_left;
    std::vector<BITLENGTH> stock_right;
    CheckDirections(&stock_left,&stock_right,target_kmer,-1);

    //junction    
    if (stock_left.size()!=1 or stock_right.size()!=1) {
        for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
            visiting.push(*itr);
            AddJunctionEdge(target_kmer,*itr,"left");
        }
        for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
            visiting.push(*itr);
            AddJunctionEdge(target_kmer,*itr,"right");
        }
        AddJunctionNode(target_kmer);
        RecordKmer(visited_kmers,target_kmer);
        return;
    }

    std::string left_part;
    BITLENGTH left_end_kmer;
    std::vector<char> extend_bases_left;
    int right_base=2*target_kmer[1]+target_kmer[0];
    left_end_kmer=ExtendLeft(stock_left[0],target_kmer,&extend_bases_left,right_base);
    for (int i=extend_bases_left.size()-1;i>=0;i--){
        left_part+=extend_bases_left[i];
    }

    std::string right_part;
    BITLENGTH right_end_kmer;
    std::vector<char> extend_bases_right;
    int left_base=2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2];
    right_end_kmer=ExtendRight(stock_right[0],target_kmer,&extend_bases_right,left_base);
    for (int i=0;i<extend_bases_right.size();i++){
        right_part+=extend_bases_right[i]; 
    }
    AddJointNode(left_end_kmer);
    AddJointNode(right_end_kmer);
    RecordKmer(visited_kmers,left_end_kmer);
    RecordKmer(visited_kmers,right_end_kmer);

    if ((left_part.size()+right_part.size())>=1){
        std::string straightnode=left_part+GetStringKmer(target_kmer)+right_part;
        AddStraightNode(straightnode);
        AddStraightEdge(left_end_kmer,straightnode,"right");
        AddStraightEdge(right_end_kmer,straightnode,"left");


        
    }


    return;
}

template<typename BITLENGTH>
BITLENGTH DeBruijnGraph<BITLENGTH>::ExtendLeft(BITLENGTH &target_kmer,BITLENGTH &previous_kmer ,std::vector<char> *extend_bases,int previous_base){

    //search eight directions(except Node before transition)
    std::vector<BITLENGTH> stock_left;
    std::vector<BITLENGTH> stock_right;
    CheckDirections(&stock_left,&stock_right,target_kmer,4+previous_base);

    if (stock_left.size()==1 and stock_right.size()==0){

        if (IsRecorded(visited_kmers,target_kmer)){
            (*extend_bases).clear();
            return target_kmer;
        }

        (*extend_bases).push_back(bit_to_base[2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2]]);
        int right_base=2*target_kmer[1]+target_kmer[0];
        return ExtendLeft(stock_left[0],target_kmer,extend_bases,right_base);
    }
    //junction
    else{
        if (!IsRecorded(visited_kmers,target_kmer)){
            for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
                visiting.push(*itr);
                AddJunctionEdge(target_kmer,*itr,"left");
            }
            for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
                visiting.push(*itr);
                AddJunctionEdge(target_kmer,*itr,"right");
            }
            AddJunctionNode(target_kmer);
            AddJunctionEdge(target_kmer,previous_kmer,"right");
            RecordKmer(visited_kmers,target_kmer);
        }
        return previous_kmer;
    }
}

template<typename BITLENGTH>
BITLENGTH DeBruijnGraph<BITLENGTH>::ExtendRight(BITLENGTH target_kmer,BITLENGTH &previous_kmer,std::vector<char> *extend_bases,int previous_base){

    //search eight directions(except Node before transition)
    std::vector<BITLENGTH> stock_left;
    std::vector<BITLENGTH> stock_right;
    CheckDirections(&stock_left,&stock_right,target_kmer,previous_base);

    if (stock_left.size()==0 and stock_right.size()==1){
        if (IsRecorded(visited_kmers,target_kmer)){
            (*extend_bases).clear();
            return target_kmer;
        }
        (*extend_bases).push_back(bit_to_base[2*target_kmer[1]+target_kmer[0]]);
        int left_base=2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2];
        return ExtendRight(stock_right[0],target_kmer,extend_bases,left_base);
    }

    //junction
    else{
        if (!IsRecorded(visited_kmers,target_kmer)){
            for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
                visiting.push(*itr);
                AddJunctionEdge(target_kmer,*itr,"left");
            }
            for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
                visiting.push(*itr);
                AddJunctionEdge(target_kmer,*itr,"right");
            }
            AddJunctionNode(target_kmer);
            AddJunctionEdge(target_kmer,previous_kmer,"left");

            RecordKmer(visited_kmers,target_kmer);
        }
        
        return previous_kmer;
    }
}

template<typename BITLENGTH>
bool DeBruijnGraph<BITLENGTH>::IsRecorded(BF<BITLENGTH> &bloomfilter,BITLENGTH &seqrching_kmer){
    BITLENGTH seqrching_kmer_bw=GetComplementKmer(seqrching_kmer);
    BITLENGTH query_kmer=CompareBit(seqrching_kmer,seqrching_kmer_bw,bitset_length);
    return bloomfilter.possiblyContains(&query_kmer,bitset_length);
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::RecordKmer(BF<BITLENGTH> &bloomfilter,BITLENGTH &seqrching_kmer){
    BITLENGTH seqrching_kmer_bw=GetComplementKmer(seqrching_kmer);
    BITLENGTH query_kmer=CompareBit(seqrching_kmer,seqrching_kmer_bw,bitset_length);
    bloomfilter.add(&query_kmer,bitset_length);
}



template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::CheckDirections(std::vector<BITLENGTH> *stock_left,std::vector<BITLENGTH> *stock_right,BITLENGTH &target_kmer,int ignored_direction){ 
    BITLENGTH back_shifted_kmer=(target_kmer>>2);
    BITLENGTH front_shifted_kmer=(target_kmer<<2);
    BITLENGTH adjacent_kmer;
    for (int i=0;i<end_bases.size();i++){
      if (i==ignored_direction) continue; 
      if (i<4){
          adjacent_kmer=(back_shifted_kmer  | end_bases[i] );  
          if (IsRecorded(*all_kmers,adjacent_kmer)) (*stock_left).push_back(adjacent_kmer);
      }
      else{
          adjacent_kmer=(front_shifted_kmer | end_bases[i] ); 
          if (IsRecorded(*all_kmers,adjacent_kmer)) (*stock_right).push_back(adjacent_kmer);
      }
    }
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddJunctionEdge(BITLENGTH &junction_node,BITLENGTH &nearby_node,std::string direction){
    BITLENGTH nearby_node_bw=GetComplementKmer(nearby_node);
    if (junction_node==nearby_node_bw){
        if (direction=="right"){
            junction_edges.push_back({junction_node,'+',junction_node,'-'}); 
        }
        else{
            junction_edges.push_back({junction_node,'-',junction_node,'+'}); 
        }
    }
    else{
        if (direction=="right"){
            junction_edges.push_back({junction_node,'+',nearby_node,'+'}); 
        }
        else{
            junction_edges.push_back({nearby_node,'+',junction_node,'+'}); 
        }
    }
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddStraightEdge(BITLENGTH &joint_node,std::string &straight_node,std::string direction){
    straight_edges.push_back({joint_node,straight_node,direction});
}



template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddJunctionNode(BITLENGTH &added_junction_node){
    junction_nodes_id++;
    std::get<0>(junction_nodes[added_junction_node])="Junction_"+std::to_string(junction_nodes_id);

    //junction_nodes[added_junction_node]=GetStringKmer(added_junction_node); //debug用

}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddJointNode(BITLENGTH &added_joint_node){
    if ( joint_nodes.find(added_joint_node)!=joint_nodes.end() ) return;
    joint_nodes_id++;
    std::get<0>(joint_nodes[added_joint_node])="Joint_"+std::to_string(joint_nodes_id);

    //joint_nodes[added_joint_node]=GetStringKmer(added_joint_node); //debug用

}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddStraightNode(std::string &added_straight_node){
    straight_nodes_id++;
    straight_nodes[added_straight_node]="Straight_"+std::to_string(straight_nodes_id);

    //straight_nodes[added_straight_node]=added_straight_node; //debug用
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::CountNodeCoverage(ReadSet &RS){
    #pragma omp parallel for  num_threads(4) 
    for(size_t b=0;b<RS.bucket_count();b++)
    for(auto bi=RS.begin(b);bi!=RS.end(b);bi++){
        std::string target_read = (bi->second);
        BITLENGTH kmer_Fw=GetFirstKmerForward<BITLENGTH>(target_read.substr(0,kmer_length));
        BITLENGTH kmer_Bw=GetFirstKmerBackward<BITLENGTH>(target_read.substr(0,kmer_length));
        AddNodeCoverage(kmer_Fw);
        AddNodeCoverage(kmer_Bw);
        
        for (int i=kmer_length;i<target_read.size();i++){
            kmer_Fw=((kmer_Fw<<2)| end_bases[ base_to_bit[ target_read[i] ] + 4 ] );
            kmer_Bw=((kmer_Bw>>2)| end_bases[ base_to_bit[ trans_base[target_read[i]] ] ] );
            AddNodeCoverage(kmer_Fw);
            AddNodeCoverage(kmer_Bw);
        }
    }
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddNodeCoverage(BITLENGTH &target_kmer){
    if (junction_nodes.find(target_kmer)!=junction_nodes.end() ){
        std::get<1>(junction_nodes[target_kmer])++;
    }
    else if (joint_nodes.find(target_kmer)!=joint_nodes.end()){
        std::get<1>(joint_nodes[target_kmer])++;
    }
}


template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::PrintGraph(){
    //std::vector< std::tuple<BITLENGTH,BITLENGTH> > junction_edges;  
    //std::vector< std::tuple<BITLENGTH,std::string,std::string> > straight_edges;  
    //std::unordered_map<std::string,std::string> straight_nodes;  
    //std::unordered_map<BITLENGTH,std::string> junction_nodes;
    //std::unordered_map<BITLENGTH,std::string> joint_nodes; 

    std::ofstream writing_gfa;
    std::string junction_file ="./de_bruijn_graph.gfa";
    writing_gfa.open(junction_file, std::ios::out);
    //header
    writing_gfa<<"H\tVN:Z:1.0"<<std::endl;
    //nodes(straight)
    for (auto itr=straight_nodes.begin();itr!=straight_nodes.end();++itr){
        writing_gfa<<"S"<<"\t"<<(itr->second)<<"\t"<<(itr->first)<<"\t"<<"KC:i:"<<1<<std::endl;
    }
    //nodes(junction)
    for (auto itr=junction_nodes.begin();itr!=junction_nodes.end();++itr){
        writing_gfa<<"S"<<"\t"<<std::get<0>(itr->second)<<"\t"<<GetStringKmer((itr->first))<<"\t"<<"KC:i:"<<std::get<1>(itr->second)*kmer_length<<std::endl;
    }
    //nodes(joint)
    for (auto itr=joint_nodes.begin();itr!=joint_nodes.end();++itr){
        writing_gfa<<"S"<<"\t"<<std::get<0>(itr->second)<<"\t"<<GetStringKmer((itr->first))<<"\t"<<"KC:i:"<<std::get<1>(itr->second)*kmer_length<<std::endl;
    }

    
    //edges(straight) 
    for (auto itr=straight_edges.begin();itr!=straight_edges.end();++itr){
        writing_gfa<<"L"<<"\t";
        if (std::get<2>(*itr)=="right"){
            writing_gfa<<std::get<0>(joint_nodes[std::get<0>(*itr)])<<"\t"<<"+"<<"\t";
            writing_gfa<<straight_nodes[std::get<1>(*itr)]<<"\t"<<"+"<<"\t";
            
        }
        else {
            writing_gfa<<straight_nodes[std::get<1>(*itr)]<<"\t"<<"+"<<"\t";
            writing_gfa<<std::get<0>(joint_nodes[std::get<0>(*itr)])<<"\t"<<"+"<<"\t";
        }
        writing_gfa<<kmer_length<<"M"<<std::endl;
    }

    //edges(junction)
    for (auto itr=junction_edges.begin();itr!=junction_edges.end();++itr){
        writing_gfa<<"L"<<"\t";

        if (junction_nodes.find(std::get<0>(*itr))!=junction_nodes.end()){
            writing_gfa<<std::get<0>(junction_nodes[std::get<0>(*itr)])<<"\t"<<std::get<1>(*itr)<<"\t";
        }
        else {
            writing_gfa<<std::get<0>(joint_nodes[std::get<0>(*itr)])<<"\t"<<std::get<1>(*itr)<<"\t";
        }
        if (junction_nodes.find(std::get<2>(*itr))!=junction_nodes.end()){
            writing_gfa<<std::get<0>(junction_nodes[std::get<2>(*itr)])<<"\t"<<std::get<3>(*itr)<<"\t";
        }
        else {
            writing_gfa<<std::get<0>(joint_nodes[std::get<2>(*itr)])<<"\t"<<std::get<3>(*itr)<<"\t";
        }

        writing_gfa<<kmer_length-1<<"M"<<std::endl;

    }
}
#endif