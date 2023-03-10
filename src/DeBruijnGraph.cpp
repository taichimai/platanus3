#ifndef DEBRUIJNGRAPH_CPP
#define DEBRUIJNGRAPH_CPP
#include"common.h"
#include"bloomfilter.cpp"

struct JunctionInfo;
struct JointInfo;
struct StraightInfo;

struct JunctionInfo{
    int id;
    int coverage=0;
    int left_kmers_cov[4]={0,0,0,0};
    int right_kmers_cov[4]={0,0,0,0};
    JunctionInfo(){}
};

struct JointInfo{
    int id;
    int coverage=0;
    StraightInfo*  connected_straight;
    JointInfo(){}
};

struct StraightInfo{
    int id;
    std::string sequence;
    JointInfo* left_joint;
    JointInfo* right_joint;
    StraightInfo(){}
};

template<typename LARGE_BITSET>
class DeBruijnGraph{
    public:
        //variables
        int kmer_length;
        int bitset_length;
        std::vector<LARGE_BITSET> end_bases;
        BF<LARGE_BITSET> *all_kmers;
        Logging *logging;

        std::unordered_map<LARGE_BITSET,JunctionInfo> junctions;
        std::unordered_map<LARGE_BITSET,JointInfo> joints;
        std::unordered_map<int,StraightInfo> straights;
        int straight_nodes_id=0;
        int junction_nodes_id=0;
        int joint_nodes_id=0;

        std::queue<LARGE_BITSET> visiting;  //for making dbg
        
        //mutex
        std::mutex mtx_junctions;
        std::mutex mtx_joints;
        std::mutex mtx_straights;
        std::mutex mtx_visiting;

        //constructor
        DeBruijnGraph(int k,BF<LARGE_BITSET> &first_bloom_filter,Logging &input_logging);

        //methods
        void MakeDBG(std::set<std::string> &seedkmer,uint64_t filtersize ,uint8_t numhashes,int threads_num);
        void SearchNode(LARGE_BITSET target_kmer);
        LARGE_BITSET ExtendLeft(LARGE_BITSET target_kmer,LARGE_BITSET previous_kmer ,std::vector<char> *extend_bases,int previous_base);
        LARGE_BITSET ExtendRight(LARGE_BITSET target_kmer,LARGE_BITSET previous_kmer,std::vector<char> *extend_bases,int previous_base);
        bool IsVisited(LARGE_BITSET &seqrching_kmer);

        bool IsRecorded(BF<LARGE_BITSET> &bloomfilter,LARGE_BITSET seqrching_kmer);
        void CheckDirections(std::vector<LARGE_BITSET> *stock_left,std::vector<LARGE_BITSET> *stock_right,LARGE_BITSET &target_kmer,int ignored_direction);

        void AddJunctionNode(LARGE_BITSET &added_node);
        void AddJointNode(LARGE_BITSET &added_node);
        void AddStraightNode(std::string &added_straight_node,LARGE_BITSET &left_joint_node,LARGE_BITSET &right_joint_node);

        void CountNodeCoverage(ReadSet &RS);
        void AddNodeCoverage(LARGE_BITSET &target_kmer);
        void PrintGraph();
};

template<typename LARGE_BITSET>
DeBruijnGraph<LARGE_BITSET>::DeBruijnGraph(int k,BF<LARGE_BITSET> &first_bloom_filter,Logging &input_logging){
    all_kmers=&first_bloom_filter;
    logging=&input_logging;
    kmer_length=k;
    bitset_length=k*2;
    LARGE_BITSET A_right(0); LARGE_BITSET A_left=(A_right<<(bitset_length-2)); 
    LARGE_BITSET C_right(1); LARGE_BITSET C_left=(C_right<<(bitset_length-2));
    LARGE_BITSET G_right(2); LARGE_BITSET G_left=(G_right<<(bitset_length-2));
    LARGE_BITSET T_right(3); LARGE_BITSET T_left=(T_right<<(bitset_length-2));
    end_bases={A_left,C_left,G_left,T_left,A_right,C_right,G_right,T_right};
}

template<typename LARGE_BITSET>
void DeBruijnGraph<LARGE_BITSET>::MakeDBG(std::set<std::string> &seedkmer,uint64_t filtersize ,uint8_t numhashes,int threads_num){
    int seedk_cnt=0;
    for (auto itr=seedkmer.begin();itr!=seedkmer.end();++itr){
        seedk_cnt++;
      
        (*logging).WriteLog("search new read"); 
        LARGE_BITSET first_kmer=GetFirstKmerForward<LARGE_BITSET>(*itr);
        if (IsVisited(first_kmer)) {
            (*logging).WriteLog("this read is visited"); 
            continue;
        }
        visiting.push(first_kmer);

        if (seedk_cnt%20!=0) continue;
        while(!visiting.empty()){
            std::vector<std::thread> ths;
            std::queue<LARGE_BITSET> mid_point;
            
            for (int i = 0; i < threads_num  && !visiting.empty(); i++) {
                LARGE_BITSET next_kmer=visiting.front();
                visiting.pop();
                mid_point.push(next_kmer);
            }

            while(!mid_point.empty()){
                LARGE_BITSET visiting_kmer=mid_point.front();
                mid_point.pop();
                ths.emplace_back(
                    [this,visiting_kmer]() { 
                        this->SearchNode(visiting_kmer);
                    }
                );
            }

            for (auto& t : ths) {
                t.join();
            }
        } 
    }
    //??????
    while(!visiting.empty()){
            std::vector<std::thread> ths;
            std::queue<LARGE_BITSET> mid_point;
            for (int i = 0; i < threads_num  && !visiting.empty(); i++) {
                LARGE_BITSET next_kmer=visiting.front();
                visiting.pop();
                mid_point.push(next_kmer);
            }
            while(!mid_point.empty()){
                LARGE_BITSET visiting_kmer=mid_point.front();
                mid_point.pop();
                ths.emplace_back(
                    [this,visiting_kmer]() {
                        this->SearchNode(visiting_kmer);
                    }
                );
            }
            for (auto& t : ths) {
                t.join();
            }
    } 
}

template<typename LARGE_BITSET>
void DeBruijnGraph<LARGE_BITSET>::SearchNode(LARGE_BITSET target_kmer){
    if (IsVisited(target_kmer)) return;

    //search eight directions
    std::vector<LARGE_BITSET> stock_left;
    std::vector<LARGE_BITSET> stock_right;
    CheckDirections(&stock_left,&stock_right,target_kmer,-1);

    //junction    
    if (stock_left.size()!=1 or stock_right.size()!=1) {
        (*logging).WriteLog("junction recording...");
        mtx_visiting.lock();
        for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
            visiting.push(*itr);
        }
        for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
            visiting.push(*itr);
        }
        mtx_visiting.unlock();
        AddJunctionNode(target_kmer);
        return;
    }
    (*logging).WriteLog("check left part");
    ////joint or straight part
    //check left part
    std::string left_part;
    LARGE_BITSET left_end_kmer;
    std::vector<char> extend_bases_left;
    int right_base=2*target_kmer[1]+target_kmer[0];
    left_end_kmer=ExtendLeft(stock_left[0],target_kmer,&extend_bases_left,right_base);
    for (int i=extend_bases_left.size()-1;i>=0;i--){
        left_part+=extend_bases_left[i];
    }
    if (IsVisited(left_end_kmer)){
        (*logging).WriteLog("this straight left part has already visited");
        return;
    } 
    
    (*logging).WriteLog("check right part");
    //check right part
    std::string right_part;
    LARGE_BITSET right_end_kmer;
    std::vector<char> extend_bases_right;
    int left_base=2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2];
    right_end_kmer=ExtendRight(stock_right[0],target_kmer,&extend_bases_right,left_base);
    for (int i=0;i<extend_bases_right.size();i++){
        right_part+=extend_bases_right[i]; 
    }
    if (IsVisited(right_end_kmer)){
      (*logging).WriteLog("this straight right part has already visited");
      return;
    } 

    //if node cannot extend
    if (left_end_kmer==right_end_kmer){
        (*logging).WriteLog("this node cannot extend");
        AddJunctionNode(left_end_kmer);
        return;
    }
    //get straight node
    if ((left_part.size()+right_part.size())>=1){
        (*logging).WriteLog("record straight node");

        std::string straightnode=left_part+GetStringKmer(target_kmer)+right_part;
        AddStraightNode(straightnode,left_end_kmer,right_end_kmer);
    }
    return;
}


template<typename LARGE_BITSET>
LARGE_BITSET DeBruijnGraph<LARGE_BITSET>::ExtendLeft(LARGE_BITSET target_kmer,LARGE_BITSET previous_kmer ,std::vector<char> *extend_bases,int previous_base){
    std::vector<LARGE_BITSET> stock_left;
    std::vector<LARGE_BITSET> stock_right;
    CheckDirections(&stock_left,&stock_right,target_kmer,4+previous_base);
    while(stock_left.size()==1 and stock_right.size()==0){
        if (IsVisited(target_kmer)){
            (*extend_bases).clear();
            return target_kmer;
        }
        (*extend_bases).push_back(bit_to_base[2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2]]);

        previous_base=2*target_kmer[1]+target_kmer[0];
        previous_kmer=target_kmer;
        target_kmer=stock_left[0];
        stock_left.clear();
        stock_right.clear();
        CheckDirections(&stock_left,&stock_right,target_kmer,4+previous_base);
    }

    if (!IsVisited(target_kmer)){
        mtx_visiting.lock();
        for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
            visiting.push(*itr);
        }
        for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
            visiting.push(*itr);
        }
        mtx_visiting.unlock();
        AddJunctionNode(target_kmer);
    }
    return previous_kmer;
}


template<typename LARGE_BITSET>
LARGE_BITSET DeBruijnGraph<LARGE_BITSET>::ExtendRight(LARGE_BITSET target_kmer,LARGE_BITSET previous_kmer,std::vector<char> *extend_bases,int previous_base){

    std::vector<LARGE_BITSET> stock_left;
    std::vector<LARGE_BITSET> stock_right;

    CheckDirections(&stock_left,&stock_right,target_kmer,previous_base);
    while(stock_left.size()==0 and stock_right.size()==1){
        if (IsVisited(target_kmer)){
            (*extend_bases).clear();
            return target_kmer;
        }
        (*extend_bases).push_back(bit_to_base[2*target_kmer[1]+target_kmer[0]]);

        previous_base=2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2];
        previous_kmer=target_kmer;
        target_kmer=stock_right[0];
        stock_left.clear();
        stock_right.clear();
        CheckDirections(&stock_left,&stock_right,target_kmer,previous_base);
    }

    if (!IsVisited(target_kmer)){
        mtx_visiting.lock();
        for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
            visiting.push(*itr);
        }
        for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
            visiting.push(*itr);
        }
        mtx_visiting.unlock();
        AddJunctionNode(target_kmer);
    }
    return previous_kmer;   
}

template<typename LARGE_BITSET>
bool DeBruijnGraph<LARGE_BITSET>::IsVisited(LARGE_BITSET &seqrching_kmer){
    bool is_visited=false;
    LARGE_BITSET seqrching_kmer_bw=GetComplementKmer(seqrching_kmer);
    std::lock(mtx_junctions, mtx_joints);
    if (junctions.find(seqrching_kmer)!=junctions.end() || junctions.find(seqrching_kmer_bw)!=junctions.end()){
        is_visited=true;
    }
    else if (joints.find(seqrching_kmer)!=joints.end() || joints.find(seqrching_kmer_bw)!=joints.end()){
        is_visited=true;
    }
    else is_visited=false;

    mtx_joints.unlock();
    mtx_junctions.unlock();
    return is_visited;
}

template<typename LARGE_BITSET>
bool DeBruijnGraph<LARGE_BITSET>::IsRecorded(BF<LARGE_BITSET> &bloomfilter,LARGE_BITSET seqrching_kmer){

    LARGE_BITSET seqrching_kmer_bw=GetComplementKmer(seqrching_kmer);
    LARGE_BITSET query_kmer=CompareBit(seqrching_kmer,seqrching_kmer_bw,bitset_length);
    return bloomfilter.possiblyContains(&query_kmer,bitset_length);
}

template<typename LARGE_BITSET>
void DeBruijnGraph<LARGE_BITSET>::CheckDirections(std::vector<LARGE_BITSET> *stock_left,std::vector<LARGE_BITSET> *stock_right,LARGE_BITSET &target_kmer,int ignored_direction){ 
    LARGE_BITSET back_shifted_kmer=(target_kmer>>2);
    LARGE_BITSET front_shifted_kmer=(target_kmer<<2);
    LARGE_BITSET adjacent_kmer;
    for (int i=0;i<end_bases.size();i++){
      if (i==ignored_direction) continue; 
      if (i<4){
          adjacent_kmer=(back_shifted_kmer  | end_bases[i] );
          if (IsRecorded(*all_kmers,adjacent_kmer)) {
                (*stock_left).push_back(adjacent_kmer);
          };
      }
      else{
          adjacent_kmer=(front_shifted_kmer | end_bases[i] ); 
          if (IsRecorded(*all_kmers,adjacent_kmer)) {
                (*stock_right).push_back(adjacent_kmer);
          }
      }
    }
}

template<typename LARGE_BITSET>
void DeBruijnGraph<LARGE_BITSET>::AddJunctionNode(LARGE_BITSET &added_node){
    if (IsVisited(added_node)){
        (*logging).WriteLog("this is visited junction node");
        return;
    }
    std::lock_guard<std::mutex> lock(mtx_junctions);
    junction_nodes_id++;
    (*logging).WriteLog("add junction node "+std::to_string(junction_nodes_id));
    junctions[added_node].id=junction_nodes_id;
}

template<typename LARGE_BITSET>
void DeBruijnGraph<LARGE_BITSET>::AddJointNode(LARGE_BITSET &added_node){
    if (IsVisited(added_node)){
        (*logging).WriteLog("this is visited joint node");
        return;
    }
    std::lock_guard<std::mutex> lock(mtx_joints);
    joint_nodes_id++;
    (*logging).WriteLog("add joint node "+std::to_string(joint_nodes_id));
    joints[added_node].id=joint_nodes_id;
}



template<typename LARGE_BITSET>
void DeBruijnGraph<LARGE_BITSET>::AddStraightNode(std::string &added_straight_node,LARGE_BITSET &left_joint_node,LARGE_BITSET &right_joint_node){
    std::lock_guard<std::mutex> lock(mtx_straights);
    if (IsVisited(left_joint_node)) return;
    AddJointNode(left_joint_node);
    AddJointNode(right_joint_node);
    straight_nodes_id++;
    (*logging).WriteLog("add straight node "+std::to_string(straight_nodes_id));
    straights[straight_nodes_id].sequence=added_straight_node;
    straights[straight_nodes_id].id=straight_nodes_id;

    //add edge
    straights[straight_nodes_id].left_joint=&joints[left_joint_node];
    straights[straight_nodes_id].right_joint=&joints[right_joint_node];
    mtx_joints.lock();
    joints[left_joint_node].connected_straight=&straights[straight_nodes_id];
    joints[right_joint_node].connected_straight=&straights[straight_nodes_id];
    mtx_joints.unlock();
}

template<typename LARGE_BITSET>
void DeBruijnGraph<LARGE_BITSET>::CountNodeCoverage(ReadSet &RS){
    #pragma omp parallel for  num_threads(20) 
    for(size_t b=0;b<RS.bucket_count();b++)
    for(auto bi=RS.begin(b);bi!=RS.end(b);bi++){
        std::string target_read = (bi->second);
        LARGE_BITSET kmer_Fw=GetFirstKmerForward<LARGE_BITSET>(target_read.substr(0,kmer_length));
        LARGE_BITSET kmer_Bw=GetFirstKmerBackward<LARGE_BITSET>(target_read.substr(0,kmer_length));
        //add node coverage
        #pragma omp critical
        {
            AddNodeCoverage(kmer_Fw);
            AddNodeCoverage(kmer_Bw);
            //add path coverage
            if (junctions.find(kmer_Fw)!=junctions.end()){
                junctions[kmer_Fw].right_kmers_cov[base_to_bit[ target_read[kmer_length] ]]++;
                
            }
            else if (junctions.find(kmer_Bw)!=junctions.end()){
                junctions[kmer_Bw].left_kmers_cov[base_to_bit[ trans_base[target_read[kmer_length]]]]++;
            }
        }
        for (int i=kmer_length;i<target_read.size();i++){
            kmer_Fw=((kmer_Fw<<2)| end_bases[ base_to_bit[ target_read[i] ] + 4 ] );
            kmer_Bw=((kmer_Bw>>2)| end_bases[ base_to_bit[ trans_base[target_read[i]] ] ] );
            #pragma omp critical
            {
                //add node coverage
                AddNodeCoverage(kmer_Fw);
                AddNodeCoverage(kmer_Bw);
                //add path coverage
                if (junctions.find(kmer_Fw)!=junctions.end()){
                    junctions[kmer_Fw].left_kmers_cov[base_to_bit[ target_read[i-kmer_length] ]]++;
                    if (i<target_read.size()-1){
                        junctions[kmer_Fw].right_kmers_cov[base_to_bit[ target_read[i+1] ]]++;
                    }
                }
                else if (junctions.find(kmer_Bw)!=junctions.end()){
                    junctions[kmer_Bw].right_kmers_cov[base_to_bit[ trans_base[target_read[i-kmer_length]]]]++;
                    if (i<target_read.size()-1){
                        junctions[kmer_Bw].left_kmers_cov[base_to_bit[ trans_base[target_read[i+1]]]]++;
                    }
                }
            }
        }
    }
}

template<typename LARGE_BITSET>
void DeBruijnGraph<LARGE_BITSET>::AddNodeCoverage(LARGE_BITSET &target_kmer){
    if (junctions.find(target_kmer)!=junctions.end()){
        junctions[target_kmer].coverage++;
    }
    if (joints.find(target_kmer)!=joints.end()) {
        joints[target_kmer].coverage++;
    }
}

template<typename LARGE_BITSET>
void DeBruijnGraph<LARGE_BITSET>::PrintGraph(){
    std::ofstream writing_gfa;
    std::string output_file ="./de_bruijn_graph.gfa";
    writing_gfa.open(output_file, std::ios::out);
    //header
    writing_gfa<<"H\tVN:Z:1.0"<<"\n";
    //nodes(straight)
    for (auto itr=straights.begin();itr!=straights.end();++itr){
        writing_gfa<<"S"<<"\t"<<"Straight_"<<(itr->first)<<"\t"<<(itr->second).sequence<<"\t"<<"KC:i:"<<(itr->second).sequence.size()<<"\n";
    }
    //nodes(junction)
    for (auto itr=junctions.begin();itr!=junctions.end();++itr){
        writing_gfa<<"S"<<"\t"<<"Junction_"<<(itr->second).id<<"\t"<<GetStringKmer((itr->first))<<"\t"<<"KC:i:"<<(itr->second).coverage*kmer_length<<"\n";
    }

    //edges(junction and straight)
    for (auto itr = junctions.begin(); itr!=junctions.end(); ++itr){
        //left
        for (int i = 0; i < 4 ; i++){
            if ((itr->second).left_kmers_cov[i]==0) continue;
            LARGE_BITSET output_left_kmer=( ((itr->first)>>2) | end_bases[i] );
            if (!IsRecorded(*all_kmers,output_left_kmer)) continue;
            if (junctions.find(output_left_kmer)!=junctions.end()){
                writing_gfa<<"L"<<"\t";
                writing_gfa<<"Junction_"<<junctions[output_left_kmer].id<<"\t+\t";
                writing_gfa<<"Junction_"<<(itr->second).id<<"\t+\t";
                writing_gfa<<kmer_length-1<<"M"<<"\n";
            }
            else if (joints.find(output_left_kmer)!=joints.end()){
                writing_gfa<<"L"<<"\t";
                writing_gfa<<"Straight_"<<(*(joints[output_left_kmer].connected_straight)).id<<"\t+\t";
                writing_gfa<<"Junction_"<<(itr->second).id<<"\t+\t";
                writing_gfa<<kmer_length-1<<"M"<<"\n";
            }
            else{
                //complement
                LARGE_BITSET output_left_kmer_bw=GetComplementKmer(output_left_kmer);
                if (junctions.find(output_left_kmer_bw)!=junctions.end()){
                    writing_gfa<<"L"<<"\t";
            
                    writing_gfa<<"Junction_"<<junctions[output_left_kmer_bw].id<<"\t-\t";
                    writing_gfa<<"Junction_"<<(itr->second).id<<"\t+\t";
                    writing_gfa<<kmer_length-1<<"M"<<"\n";
                   
                }
                else if (joints.find(output_left_kmer_bw)!=joints.end()){
                    writing_gfa<<"L"<<"\t";
                  
                    writing_gfa<<"Straight_"<<(*(joints[output_left_kmer_bw].connected_straight)).id<<"\t-\t";
                    writing_gfa<<"Junction_"<<(itr->second).id<<"\t+\t";
                    writing_gfa<<kmer_length-1<<"M"<<"\n";
                  
                }
            }
        }
      
        //right
        for (int i = 0; i < 4 ; i++){
            if ((itr->second).right_kmers_cov[i]==0) continue;
            LARGE_BITSET output_right_kmer=( ((itr->first)<<2) | end_bases[i+4] );
            if (!IsRecorded(*all_kmers,output_right_kmer)) continue;
            if (junctions.find(output_right_kmer)!=junctions.end()){
                writing_gfa<<"L"<<"\t";
                writing_gfa<<"Junction_"<<(itr->second).id<<"\t+\t";
                writing_gfa<<"Junction_"<<junctions[output_right_kmer].id<<"\t+\t";
                writing_gfa<<kmer_length-1<<"M"<<"\n";
            }
            else if (joints.find(output_right_kmer)!=joints.end()) {
                //joint
                writing_gfa<<"L"<<"\t";
                writing_gfa<<"Junction_"<<(itr->second).id<<"\t+\t";
                writing_gfa<<"Straight_"<<(*(joints[output_right_kmer].connected_straight)).id<<"\t+\t";
                writing_gfa<<kmer_length-1<<"M"<<"\n";
            }
            else{
                //complement
                LARGE_BITSET output_right_kmer_bw=GetComplementKmer(output_right_kmer);
                if (junctions.find(output_right_kmer_bw)!=junctions.end()){
                    writing_gfa<<"L"<<"\t";
                    writing_gfa<<"Junction_"<<(itr->second).id<<"\t+\t";
                    writing_gfa<<"Junction_"<<junctions[output_right_kmer_bw].id<<"\t-\t";
                    writing_gfa<<kmer_length-1<<"M"<<"\n";
                }
                else if (joints.find(output_right_kmer_bw)!=joints.end()){
                    writing_gfa<<"L"<<"\t";
                    writing_gfa<<"Junction_"<<(itr->second).id<<"\t+\t";
                    writing_gfa<<"Straight_"<<(*(joints[output_right_kmer_bw].connected_straight)).id<<"\t-\t";
                    writing_gfa<<kmer_length-1<<"M"<<"\n";
                }
            }
        }
    }
}
#endif