#include"DeBruijnGraph.h"

template<typename BITLENGTH>
DeBruijnGraph<BITLENGTH>::DeBruijnGraph(int k){
    bitset_length=k*2;

    BITLENGTH A_right(0); BITLENGTH A_left=(A_right<<(bitset_length-2)); 
    BITLENGTH C_right(1); BITLENGTH C_left=(C_right<<(bitset_length-2));
    BITLENGTH G_right(2); BITLENGTH G_left=(G_right<<(bitset_length-2));
    BITLENGTH T_right(3); BITLENGTH T_left=(T_right<<(bitset_length-2));
    end_bases={A_left,C_left,G_left,T_left,A_right,C_right,G_right,T_right};
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::MakeDBG(KmerSet &seedkmer,BF<BITLENGTH> &all_kmers,uint64_t filtersize ,uint8_t numhashes){

    BF<BITLENGTH> visited_kmers(filtersize,numhashes);
    std::queue<BITLENGTH> visiting;

    for (auto itr=seedkmer.begin();itr!=seedkmer.end();++itr){
        BITLENGTH fw_kmer=GetFirstKmerForward<BITLENGTH>(*itr);
        BITLENGTH bw_kmer=GetFirstKmerBackward<BITLENGTH>(*itr);
        BITLENGTH target_kmer=CompareBit(fw_kmer,bw_kmer,bitset_length);
    
        if (!visited_kmers.possiblyContains(&target_kmer,bitset_length)) continue;
        visiting.push(target_kmer);

        while(!visiting.empty()){
            BITLENGTH visiting_kmer=visiting.front();
            visiting.pop();
            SearchNode(visiting_kmer,&visited_kmers,&visiting,all_kmers);
        }   
    }
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::SearchNode(BITLENGTH target_kmer,BF<BITLENGTH> *visited_kmers,std::queue<BITLENGTH> *visiting,BF<BITLENGTH> &all_kmers){

    if ((*visited_kmers).possiblyContains(&target_kmer,bitset_length)) return;

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
            right_fw_kmer=(front_shifted_kmer | end_bases[i] ); 
            right_bw_kmer=GetComplementKmer(right_fw_kmer);
            right_search_kmer=CompareBit(right_fw_kmer,right_bw_kmer,bitset_length);

            if (all_kmers.possiblyContains(&right_search_kmer,bitset_length)) stock_right.push_back(right_fw_kmer);
        }
    }

    if (stock_left.size()!=1 or stock_left.size()!=1) {
        for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
            (*visiting).push(*itr);
            junction_edges.push_back({*itr,target_kmer}); 
        }
        for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
            (*visiting).push(*itr); 
            junction_edges.push_back({target_kmer,*itr}); 
        }
        AddJunctionNode(target_kmer);
        (*visited_kmers).add(&target_kmer,bitset_length);
        return;
    }

    std::string left_part;
    BITLENGTH left_end_kmer;
    std::vector<char> extend_bases_left;
    int right_base=2*target_kmer[1]+target_kmer[0];
    left_end_kmer=ExtendLeft(stock_left[0],target_kmer,&extend_bases_left,right_base,all_kmers,visited_kmers,visiting);
    for (int i=extend_bases_left.size()-1;i>=0;i--){
        left_part+=extend_bases_left[i];
    }
    
    std::string right_part;
    BITLENGTH right_end_kmer;
    std::vector<char> extend_bases_right;
    int left_base=2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2];
    right_end_kmer=ExtendRight(stock_right[0],target_kmer,&extend_bases_right,left_base,all_kmers,visited_kmers,visiting);
    for (int i=0;i<extend_bases_right.size();i++){
        right_part+=extend_bases_right[i]; 
    }

    if ((left_part.size()+right_part.size())>=1){
        std::string straightnode=left_part+GetStringKmer(target_kmer)+right_part;
        AddStraightNode(straightnode);
        straight_edges.push_back({left_end_kmer,straightnode,"right"}); 
        straight_edges.push_back({right_end_kmer,straightnode,"left"}); 
    }
    return;
}

template<typename BITLENGTH>
BITLENGTH DeBruijnGraph<BITLENGTH>::ExtendLeft(BITLENGTH &target_kmer,BITLENGTH &previous_kmer ,std::vector<char> *extend_bases,int previous_base,BF<BITLENGTH> &all_kmers,BF<BITLENGTH> *visited_kmers,std::queue<BITLENGTH> *visiting){

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
        if ((*visited_kmers).possiblyContains(&target_kmer,bitset_length)){
            (*extend_bases).clear();
            return target_kmer;
        }

        (*extend_bases).push_back(bit_to_base[2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2]]);
        int right_base=2*target_kmer[1]+target_kmer[0];
        return ExtendLeft(stock_left[0],target_kmer,extend_bases,right_base,all_kmers,visited_kmers,visiting);
    }
    else{
        if (!(*visited_kmers).possiblyContains(&target_kmer,bitset_length)){
            for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
                (*visiting).push(*itr);
                junction_edges.push_back({*itr,target_kmer}); 
            }
            for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
                (*visiting).push(*itr);
                junction_edges.push_back({target_kmer,*itr});
            }
            junction_edges.push_back({target_kmer,previous_kmer});
            (*visited_kmers).add(&target_kmer,bitset_length);
        }
        (*visited_kmers).add(&previous_kmer,bitset_length);
        return previous_kmer;
    }
}

template<typename BITLENGTH>
BITLENGTH DeBruijnGraph<BITLENGTH>::ExtendRight(BITLENGTH target_kmer,BITLENGTH &previous_kmer,std::vector<char> *extend_bases,int previous_base,BF<BITLENGTH> &all_kmers,BF<BITLENGTH> *visited_kmers,std::queue<BITLENGTH> *visiting){

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
        if ((*visited_kmers).possiblyContains(&target_kmer,bitset_length)){
            (*extend_bases).clear();
            return target_kmer;
        }
        (*extend_bases).push_back(bit_to_base[2*target_kmer[1]+target_kmer[0]]);
        int left_base=2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2];
        return ExtendRight(stock_right[0],target_kmer,extend_bases,left_base,all_kmers,visited_kmers,visiting);
    }
    else{

        if (!(*visited_kmers).possiblyContains(&target_kmer,bitset_length)){
            for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
                (*visiting).push(*itr);
                junction_edges.push_back({*itr,target_kmer});
            }
            for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
                (*visiting).push(*itr);
                junction_edges.push_back({target_kmer,*itr}); 
            }
            junction_edges.push_back({previous_kmer,target_kmer});
            (*visited_kmers).add(&target_kmer,bitset_length);
        }
        (*visited_kmers).add(&previous_kmer,bitset_length);
        return previous_kmer;
    }
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddJunctionNode(BITLENGTH &added_junction_node){
    junction_nodes_id++;
    junction_nodes[added_junction_node]="Junction_"+std::to_string(junction_nodes_id);
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddStraightNode(std::string &added_straight_node){
    straight_nodes_id++;
    straight_nodes[added_straight_node]="Straight_"+std::to_string(straight_nodes_id);
}