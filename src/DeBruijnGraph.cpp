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
      std::vector< std::tuple<BITLENGTH,BITLENGTH> > junction_edges;  
      std::vector< std::tuple<BITLENGTH,std::string,std::string> > straight_edges;  

      std::unordered_map<std::string,std::string> straight_nodes;  
      std::unordered_map<BITLENGTH,std::string> junction_nodes;
      std::unordered_map<BITLENGTH,std::string> joint_nodes;  

      id_counter straight_nodes_id;
      id_counter junction_nodes_id;
      id_counter joint_nodes_id;

      BF<BITLENGTH> visited_kmers;     //for making dbg
      std::queue<BITLENGTH> visiting;  //for making dbg

      //constructor
      DeBruijnGraph(int k);

      //methods
      void MakeDBG(KmerSet &seedkmer,BF<BITLENGTH> &all_kmers,uint64_t filtersize ,uint8_t numhashes);
      void SearchNode(BITLENGTH target_kmer,BF<BITLENGTH> &all_kmers);
      BITLENGTH ExtendLeft(BITLENGTH &target_kmer,BITLENGTH &previous_kmer ,std::vector<char> *extend_bases,int previous_base,BF<BITLENGTH> &all_kmers);
      BITLENGTH ExtendRight(BITLENGTH target_kmer,BITLENGTH &previous_kmer,std::vector<char> *extend_bases,int previous_base,BF<BITLENGTH> &all_kmers);
      void AddJunctionNode(BITLENGTH &added_junction_node);
      void AddJointNode(BITLENGTH &added_joint_node);
      void AddStraightNode(std::string &added_straight_node);
      void PrintGraph();
      
};

template<typename BITLENGTH>
DeBruijnGraph<BITLENGTH>::DeBruijnGraph(int k){
    kmer_length=k;
    bitset_length=k*2;
    BITLENGTH A_right(0); BITLENGTH A_left=(A_right<<(bitset_length-2)); 
    BITLENGTH C_right(1); BITLENGTH C_left=(C_right<<(bitset_length-2));
    BITLENGTH G_right(2); BITLENGTH G_left=(G_right<<(bitset_length-2));
    BITLENGTH T_right(3); BITLENGTH T_left=(T_right<<(bitset_length-2));
    end_bases={A_left,C_left,G_left,T_left,A_right,C_right,G_right,T_right};
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::MakeDBG(KmerSet &seedkmer,BF<BITLENGTH> &all_kmers,uint64_t filtersize ,uint8_t numhashes){

    visited_kmers.Set_BF(filtersize ,numhashes);

    for (auto itr=seedkmer.begin();itr!=seedkmer.end();++itr){
        BITLENGTH fw_kmer=GetFirstKmerForward<BITLENGTH>(*itr);
        BITLENGTH bw_kmer=GetFirstKmerBackward<BITLENGTH>(*itr);
        BITLENGTH target_kmer=CompareBit(fw_kmer,bw_kmer,bitset_length);

        if (visited_kmers.possiblyContains(&target_kmer,bitset_length)) continue;
        visiting.push(target_kmer);

        while(!visiting.empty()){
            BITLENGTH visiting_kmer=visiting.front();
            std::cout<<GetStringKmer(visiting_kmer)<<std::endl;
            visiting.pop();
            SearchNode(visiting_kmer,all_kmers);
        } 
    }

    visited_kmers.Delete_Filter();
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::SearchNode(BITLENGTH target_kmer,BF<BITLENGTH> &all_kmers){

    if (visited_kmers.possiblyContains(&target_kmer,bitset_length)) return;

    //search eight directions
    std::vector<BITLENGTH> stock_left;
    BITLENGTH back_shifted_kmer=(target_kmer>>2);
    BITLENGTH left_fw_kmer;
    BITLENGTH left_bw_kmer;
    BITLENGTH left_search_kmer;
    std::vector<BITLENGTH> stock_right;
    BITLENGTH front_shifted_kmer=(target_kmer<<2);
    BITLENGTH right_fw_kmer;
    BITLENGTH right_bw_kmer;
    BITLENGTH right_search_kmer;
    
    for (int i=0;i<end_bases.size();i++){
        if (i<4){
            left_fw_kmer=(back_shifted_kmer | end_bases[i] );  
            left_bw_kmer=GetComplementKmer(left_fw_kmer);
            left_search_kmer=CompareBit(left_fw_kmer,left_bw_kmer,bitset_length);
            if (all_kmers.possiblyContains(&left_search_kmer,bitset_length)) {
                stock_left.push_back(left_fw_kmer);
            }
        }
        else{
            right_fw_kmer=(front_shifted_kmer | end_bases[i] ); 
            right_bw_kmer=GetComplementKmer(right_fw_kmer);
            right_search_kmer=CompareBit(right_fw_kmer,right_bw_kmer,bitset_length);
            if (all_kmers.possiblyContains(&right_search_kmer,bitset_length)) {
                stock_right.push_back(right_fw_kmer);
            }
        }
    }

    //junction    
    if (stock_left.size()!=1 or stock_right.size()!=1) {
        for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
            visiting.push(*itr);
           
            junction_edges.push_back({*itr,target_kmer}); 
        }
        for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
            visiting.push(*itr);
            junction_edges.push_back({target_kmer,*itr}); 
        }
        AddJunctionNode(target_kmer);
        visited_kmers.add(&target_kmer,bitset_length);
        return;
    }

    std::string left_part;
    BITLENGTH left_end_kmer;
    std::vector<char> extend_bases_left;
    int right_base=2*target_kmer[1]+target_kmer[0];
    left_end_kmer=ExtendLeft(stock_left[0],target_kmer,&extend_bases_left,right_base,all_kmers);
    for (int i=extend_bases_left.size()-1;i>=0;i--){
        left_part+=extend_bases_left[i];
    }

    std::string right_part;
    BITLENGTH right_end_kmer;
    std::vector<char> extend_bases_right;
    int left_base=2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2];
    right_end_kmer=ExtendRight(stock_right[0],target_kmer,&extend_bases_right,left_base,all_kmers);
    for (int i=0;i<extend_bases_right.size();i++){
        right_part+=extend_bases_right[i]; 
    }
    //std::cout<<"left "<<"end "<<GetStringKmer(left_end_kmer)<<std::endl;
    //std::cout<<"right "<<"end "<<GetStringKmer(right_end_kmer)<<std::endl;
    AddJointNode(left_end_kmer);
    AddJointNode(right_end_kmer);

    if ((left_part.size()+right_part.size())>=1){
        std::string straightnode=left_part+GetStringKmer(target_kmer)+right_part;
        AddStraightNode(straightnode);
        straight_edges.push_back({left_end_kmer,straightnode,"right"}); 
        straight_edges.push_back({right_end_kmer,straightnode,"left"}); 
    }


    return;
}

template<typename BITLENGTH>
BITLENGTH DeBruijnGraph<BITLENGTH>::ExtendLeft(BITLENGTH &target_kmer,BITLENGTH &previous_kmer ,std::vector<char> *extend_bases,int previous_base,BF<BITLENGTH> &all_kmers){

    //search eight directions(except Node before transition)
    std::vector<BITLENGTH> stock_left;
    BITLENGTH back_shifted_kmer=(target_kmer>>2);
    BITLENGTH left_fw_kmer;
    BITLENGTH left_bw_kmer;
    BITLENGTH left_search_kmer;

    std::vector<BITLENGTH> stock_right;
    BITLENGTH front_shifted_kmer=(target_kmer<<2);
    BITLENGTH right_fw_kmer;
    BITLENGTH right_bw_kmer;
    BITLENGTH right_search_kmer;

    for (int i=0;i<end_bases.size();i++){
        if (i<4){
            left_fw_kmer=(back_shifted_kmer | end_bases[i] );  
            left_bw_kmer=GetComplementKmer(left_fw_kmer);
            left_search_kmer=CompareBit(left_fw_kmer,left_bw_kmer,bitset_length);
            if (all_kmers.possiblyContains(&left_search_kmer,bitset_length)) stock_left.push_back(left_fw_kmer);
        }
        else{
            if (i==4+previous_base) continue; 
            right_fw_kmer=(front_shifted_kmer | end_bases[i] ); 
            right_bw_kmer=GetComplementKmer(right_fw_kmer);
            right_search_kmer=CompareBit(right_fw_kmer,right_bw_kmer,bitset_length);
            if (all_kmers.possiblyContains(&right_search_kmer,bitset_length)) stock_right.push_back(right_fw_kmer);
        }
    }

    if (stock_left.size()==1 and stock_right.size()==0){
        BITLENGTH target_kmer_bw=GetComplementKmer(target_kmer);
        BITLENGTH query_kmer=CompareBit(target_kmer,target_kmer_bw,bitset_length);
        if (visited_kmers.possiblyContains(&query_kmer,bitset_length)){
            (*extend_bases).clear();
            return target_kmer;
        }

        (*extend_bases).push_back(bit_to_base[2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2]]);
        int right_base=2*target_kmer[1]+target_kmer[0];
        return ExtendLeft(stock_left[0],target_kmer,extend_bases,right_base,all_kmers);
    }
    //junction
    else{
        BITLENGTH target_kmer_bw=GetComplementKmer(target_kmer);
        BITLENGTH query_kmer=CompareBit(target_kmer,target_kmer_bw,bitset_length);
        if (!visited_kmers.possiblyContains(&query_kmer,bitset_length)){
            for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
                visiting.push(*itr);
                junction_edges.push_back({*itr,target_kmer}); 
            }
            for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
                visiting.push(*itr);
                junction_edges.push_back({target_kmer,*itr});
            }
            AddJunctionNode(target_kmer);
            junction_edges.push_back({target_kmer,previous_kmer});
            visited_kmers.add(&target_kmer,bitset_length);
        }
        visited_kmers.add(&previous_kmer,bitset_length);
        return previous_kmer;
    }
}

template<typename BITLENGTH>
BITLENGTH DeBruijnGraph<BITLENGTH>::ExtendRight(BITLENGTH target_kmer,BITLENGTH &previous_kmer,std::vector<char> *extend_bases,int previous_base,BF<BITLENGTH> &all_kmers){

    //search eight directions(except Node before transition)
    std::vector<BITLENGTH> stock_left;
    BITLENGTH back_shifted_kmer=(target_kmer>>2);
    BITLENGTH left_fw_kmer;
    BITLENGTH left_bw_kmer;
    BITLENGTH left_search_kmer;

    std::vector<BITLENGTH> stock_right;
    BITLENGTH front_shifted_kmer=(target_kmer<<2);
    BITLENGTH right_fw_kmer;
    BITLENGTH right_bw_kmer;
    BITLENGTH right_search_kmer;

    for (int i=0;i<end_bases.size();i++){
        if (i<4){
            if (i==previous_base) continue; 
            left_fw_kmer=(back_shifted_kmer | end_bases[i] );  
            left_bw_kmer=GetComplementKmer(left_fw_kmer);
            left_search_kmer=CompareBit(left_fw_kmer,left_bw_kmer,bitset_length);
            if (all_kmers.possiblyContains(&left_search_kmer,bitset_length)) stock_left.push_back(left_fw_kmer);
        }
        else{
            right_fw_kmer=(front_shifted_kmer | end_bases[i] ); 
            right_bw_kmer=GetComplementKmer(right_fw_kmer);
            right_search_kmer=CompareBit(right_fw_kmer,right_bw_kmer,bitset_length);
            if (all_kmers.possiblyContains(&right_search_kmer,bitset_length)) stock_right.push_back(right_fw_kmer);
        }
    }

    if (stock_left.size()==0 and stock_right.size()==1){
        BITLENGTH target_kmer_bw=GetComplementKmer(target_kmer);
        BITLENGTH query_kmer=CompareBit(target_kmer,target_kmer_bw,bitset_length);
        if (visited_kmers.possiblyContains(&query_kmer,bitset_length)){
            (*extend_bases).clear();
            return target_kmer;
        }
        (*extend_bases).push_back(bit_to_base[2*target_kmer[1]+target_kmer[0]]);
        int left_base=2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2];
        return ExtendRight(stock_right[0],target_kmer,extend_bases,left_base,all_kmers);
    }

    //junction
    else{
        BITLENGTH target_kmer_bw=GetComplementKmer(target_kmer);
        BITLENGTH query_kmer=CompareBit(target_kmer,target_kmer_bw,bitset_length);
        if (!visited_kmers.possiblyContains(&query_kmer,bitset_length)){
            for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
                visiting.push(*itr);
                junction_edges.push_back({*itr,target_kmer});
            }
            for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
                visiting.push(*itr);
                junction_edges.push_back({target_kmer,*itr}); 
            }
            AddJunctionNode(target_kmer);
            junction_edges.push_back({previous_kmer,target_kmer});
            visited_kmers.add(&target_kmer,bitset_length);
        }
        visited_kmers.add(&previous_kmer,bitset_length);
        return previous_kmer;
    }
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddJunctionNode(BITLENGTH &added_junction_node){
    junction_nodes_id++;
    junction_nodes[added_junction_node]="Junction_"+std::to_string(junction_nodes_id);

    //junction_nodes[added_junction_node]=GetStringKmer(added_junction_node); //debug用

}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddJointNode(BITLENGTH &added_joint_node){
    if ( joint_nodes.find(added_joint_node)!=joint_nodes.end() ) return;
    joint_nodes_id++;
    joint_nodes[added_joint_node]="Joint_"+std::to_string(joint_nodes_id);

    //joint_nodes[added_joint_node]=GetStringKmer(added_joint_node); //debug用

}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddStraightNode(std::string &added_straight_node){
    straight_nodes_id++;
    straight_nodes[added_straight_node]="Straight_"+std::to_string(straight_nodes_id);

    //straight_nodes[added_straight_node]=added_straight_node; //debug用
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
        writing_gfa<<"S"<<"\t"<<(itr->second)<<"\t"<<(itr->first)<<std::endl;
    }
    //nodes(junction)
    for (auto itr=junction_nodes.begin();itr!=junction_nodes.end();++itr){
        writing_gfa<<"S"<<"\t"<<(itr->second)<<"\t"<<GetStringKmer((itr->first))<<std::endl;
    }
    //nodes(joint)
    for (auto itr=joint_nodes.begin();itr!=joint_nodes.end();++itr){
        writing_gfa<<"S"<<"\t"<<(itr->second)<<"\t"<<GetStringKmer((itr->first))<<std::endl;
    }

    
    //edges(straight) 
    for (auto itr=straight_edges.begin();itr!=straight_edges.end();++itr){
        writing_gfa<<"L"<<"\t";
        if (std::get<2>(*itr)=="right"){
            //writing_gfa<<GetStringKmer(std::get<0>(*itr))<<"\t"<<"+"<<"\t";
            //writing_gfa<<std::get<1>(*itr)<<"\t"<<"+"<<"\t";
            writing_gfa<<joint_nodes[std::get<0>(*itr)]<<"\t"<<"+"<<"\t";
            writing_gfa<<straight_nodes[std::get<1>(*itr)]<<"\t"<<"+"<<"\t";
            
        }
        else {
            //writing_gfa<<std::get<1>(*itr)<<"\t"<<"+"<<"\t";
            //writing_gfa<<GetStringKmer(std::get<0>(*itr))<<"\t"<<"+"<<"\t";
            writing_gfa<<straight_nodes[std::get<1>(*itr)]<<"\t"<<"+"<<"\t";
            writing_gfa<<joint_nodes[std::get<0>(*itr)]<<"\t"<<"+"<<"\t";
        }
        writing_gfa<<kmer_length<<"M"<<std::endl;
    }

    //edges(junction)
    for (auto itr=junction_edges.begin();itr!=junction_edges.end();++itr){
        writing_gfa<<"L"<<"\t";
        //writing_gfa<<GetStringKmer(std::get<0>(*itr))<<"\t"<<"+"<<"\t";
        //writing_gfa<<GetStringKmer(std::get<1>(*itr))<<"\t"<<"+"<<"\t";

        if (junction_nodes.find(std::get<0>(*itr))!=junction_nodes.end()){
            writing_gfa<<junction_nodes[std::get<0>(*itr)]<<"\t"<<"+"<<"\t";
        }
        else {
            writing_gfa<<joint_nodes[std::get<0>(*itr)]<<"\t"<<"+"<<"\t";
        }
        if (junction_nodes.find(std::get<1>(*itr))!=junction_nodes.end()){
            writing_gfa<<junction_nodes[std::get<1>(*itr)]<<"\t"<<"+"<<"\t";
        }
        else {
            writing_gfa<<joint_nodes[std::get<1>(*itr)]<<"\t"<<"+"<<"\t";
        }

        writing_gfa<<kmer_length-1<<"M"<<std::endl;

    }
}
#endif