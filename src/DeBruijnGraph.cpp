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
    int connect_direction; // direction straight is connected in
    JointInfo(){}
};

struct StraightInfo{
    int id;
    std::string sequence;
    JointInfo* left_joint;
    JointInfo* right_joint;
    StraightInfo(){}
};

template<typename BITLENGTH>
class DeBruijnGraph{
    public:
        //variables
        int kmer_length;
        int bitset_length;
        std::vector<BITLENGTH> end_bases;

        std::unordered_map<BITLENGTH,JunctionInfo> junctions;
        std::unordered_map<BITLENGTH,JointInfo> joints;
        std::unordered_map<int,StraightInfo> straights;

        int straight_nodes_id=0;
        int junction_nodes_id=0;
        int joint_nodes_id=0;

        BF<BITLENGTH> *all_kmers;
        BF<BITLENGTH> visited_kmers;     //for making dbg
        std::queue<BITLENGTH> visiting;  //for making dbg

        //constructor
        DeBruijnGraph(int k,BF<BITLENGTH> &first_bloom_filter);

        //methods
        void MakeDBG(KmerSet &seedkmer,uint64_t filtersize ,uint8_t numhashes);
        void SearchNode(BITLENGTH &target_kmer);
        BITLENGTH ExtendLeft(BITLENGTH &target_kmer,BITLENGTH &previous_kmer ,std::vector<char> *extend_bases,int previous_base);
        BITLENGTH ExtendRight(BITLENGTH &target_kmer,BITLENGTH &previous_kmer,std::vector<char> *extend_bases,int previous_base);
        bool IsVisited(BITLENGTH &seqrching_kmer);

        bool IsRecorded(BF<BITLENGTH> &bloomfilter,BITLENGTH &seqrching_kmer);
        void CheckDirections(std::vector<BITLENGTH> *stock_left,std::vector<BITLENGTH> *stock_right,BITLENGTH &target_kmer,int ignored_direction);

        void AddStraightEdge(BITLENGTH &left_joint_node,BITLENGTH &right_joint_node);
        void AddJunctionNode(BITLENGTH &added_node);
        void AddJointNode(BITLENGTH &added_node);
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

    for (auto itr=seedkmer.begin();itr!=seedkmer.end();++itr){
        std::cout<<"search new read"<<std::endl;
        BITLENGTH fw_kmer=GetFirstKmerForward<BITLENGTH>(*itr);
        if (IsVisited(fw_kmer)) {
            std::cout<<"this read is visited"<<std::endl;
            continue;
        }
        visiting.push(fw_kmer);
        while(!visiting.empty()){
            BITLENGTH visiting_kmer=visiting.front();
            visiting.pop();
            SearchNode(visiting_kmer);
        } 
    }
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::SearchNode(BITLENGTH &target_kmer){

    if (IsVisited(target_kmer)) return;

    //search eight directions
    std::vector<BITLENGTH> stock_left;
    std::vector<BITLENGTH> stock_right;
    CheckDirections(&stock_left,&stock_right,target_kmer,-1);

    //junction    
    if (stock_left.size()!=1 or stock_right.size()!=1) {
        for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
            visiting.push(*itr);
        }
        for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
            visiting.push(*itr);
        }
        AddJunctionNode(target_kmer);
        return;
    }

    ////joint or straight part
    //check left part
    std::string left_part;
    BITLENGTH left_end_kmer;
    std::vector<char> extend_bases_left;
    int right_base=2*target_kmer[1]+target_kmer[0];
    left_end_kmer=ExtendLeft(stock_left[0],target_kmer,&extend_bases_left,right_base);
    for (int i=extend_bases_left.size()-1;i>=0;i--){
        left_part+=extend_bases_left[i];
    }
    if (IsVisited(left_end_kmer)){
        std::cout<<"this straight left part has already visited"<<std::endl;
        return;
    } 

    //check right part
    std::string right_part;
    BITLENGTH right_end_kmer;
    std::vector<char> extend_bases_right;
    int left_base=2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2];
    right_end_kmer=ExtendRight(stock_right[0],target_kmer,&extend_bases_right,left_base);
    for (int i=0;i<extend_bases_right.size();i++){
        right_part+=extend_bases_right[i]; 
    }
    if (IsVisited(right_end_kmer)){
      std::cout<<"this straight right part has already visited"<<std::endl;
      return;
    } 

    //if node cannot extend
    if (left_end_kmer==right_end_kmer){
        std::cout<<"this node cannot extend"<<std::endl;
        AddJunctionNode(left_end_kmer);
        return;
    }
    //get straight node
    if ((left_part.size()+right_part.size())>=1){
        std::cout<<"extend node"<<std::endl;
        std::string straightnode=left_part+GetStringKmer(target_kmer)+right_part;
        AddJointNode(left_end_kmer);
        AddJointNode(right_end_kmer);
        AddStraightNode(straightnode);
        AddStraightEdge(left_end_kmer,right_end_kmer);
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
        //if visited
        if (IsVisited(target_kmer)){
            (*extend_bases).clear();
            return target_kmer;
        }

        (*extend_bases).push_back(bit_to_base[2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2]]);
        int right_base=2*target_kmer[1]+target_kmer[0];
        return ExtendLeft(stock_left[0],target_kmer,extend_bases,right_base);
    }
    //junction
    else{
        if (!IsVisited(target_kmer)){
            for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
                visiting.push(*itr);
            }
            for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
                visiting.push(*itr);
            }
            AddJunctionNode(target_kmer);
        }
        return previous_kmer;
    }
}

template<typename BITLENGTH>
BITLENGTH DeBruijnGraph<BITLENGTH>::ExtendRight(BITLENGTH &target_kmer,BITLENGTH &previous_kmer,std::vector<char> *extend_bases,int previous_base){

    //search eight directions(except Node before transition)
    std::vector<BITLENGTH> stock_left;
    std::vector<BITLENGTH> stock_right;
    CheckDirections(&stock_left,&stock_right,target_kmer,previous_base);

    if (stock_left.size()==0 and stock_right.size()==1){
        //if visited
        if (IsVisited(target_kmer)){
            (*extend_bases).clear();
            return target_kmer;
        }
        (*extend_bases).push_back(bit_to_base[2*target_kmer[1]+target_kmer[0]]);
        int left_base=2*target_kmer[bitset_length-1]+target_kmer[bitset_length-2];
        return ExtendRight(stock_right[0],target_kmer,extend_bases,left_base);
    }
    //junction
    else{
        if (!IsVisited(target_kmer)){
            for (auto itr=stock_left.begin();itr!=stock_left.end();++itr){
                visiting.push(*itr);
            }
            for (auto itr=stock_right.begin();itr!=stock_right.end();++itr){
                visiting.push(*itr);
            }
            AddJunctionNode(target_kmer);
        }
        return previous_kmer;
    }
}

template<typename BITLENGTH>
bool DeBruijnGraph<BITLENGTH>::IsVisited(BITLENGTH &seqrching_kmer){
    BITLENGTH seqrching_kmer_bw=GetComplementKmer(seqrching_kmer);
    if (junctions.find(seqrching_kmer)!=junctions.end() || junctions.find(seqrching_kmer_bw)!=junctions.end()){
        return true;
    }
    else if (joints.find(seqrching_kmer)!=joints.end() || joints.find(seqrching_kmer_bw)!=joints.end()){
        return true;
    }
    else return false;
}

template<typename BITLENGTH>
bool DeBruijnGraph<BITLENGTH>::IsRecorded(BF<BITLENGTH> &bloomfilter,BITLENGTH &seqrching_kmer){
    BITLENGTH seqrching_kmer_bw=GetComplementKmer(seqrching_kmer);
    BITLENGTH query_kmer=CompareBit(seqrching_kmer,seqrching_kmer_bw,bitset_length);
    return bloomfilter.possiblyContains(&query_kmer,bitset_length);
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
void DeBruijnGraph<BITLENGTH>::AddStraightEdge(BITLENGTH &left_joint_node,BITLENGTH &right_joint_node){
    straights[straight_nodes_id].left_joint=&joints[left_joint_node];
    straights[straight_nodes_id].right_joint=&joints[right_joint_node];
    joints[left_joint_node].connected_straight=&straights[straight_nodes_id];
    joints[left_joint_node].connect_direction=1;
    joints[right_joint_node].connected_straight=&straights[straight_nodes_id];
    joints[right_joint_node].connect_direction=0;
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddJunctionNode(BITLENGTH &added_node){
    if (IsVisited(added_node)){
        std::cout<<"this is visited junction node"<<std::endl;
        return;
    }
    junction_nodes_id++;
    std::cout<<"add junction node "<<junction_nodes_id<<std::endl;
    junctions[added_node].id=junction_nodes_id;
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddJointNode(BITLENGTH &added_node){
    if (IsVisited(added_node)){
        std::cout<<"this is visited joint node"<<std::endl;
        return;
    }
    joint_nodes_id++;
    std::cout<<"add joint node "<<joint_nodes_id<<std::endl;
    joints[added_node].id=joint_nodes_id;
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddStraightNode(std::string &added_straight_node){
    straight_nodes_id++;
    std::cout<<"add straight node"<<straight_nodes_id<<std::endl;
    straights[straight_nodes_id].sequence=added_straight_node;
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::CountNodeCoverage(ReadSet &RS){
    #pragma omp parallel for  num_threads(20) 
    for(size_t b=0;b<RS.bucket_count();b++)
    for(auto bi=RS.begin(b);bi!=RS.end(b);bi++){
        std::string target_read = (bi->second);
        BITLENGTH kmer_Fw=GetFirstKmerForward<BITLENGTH>(target_read.substr(0,kmer_length));
        BITLENGTH kmer_Bw=GetFirstKmerBackward<BITLENGTH>(target_read.substr(0,kmer_length));
        //add node coverage
        AddNodeCoverage(kmer_Fw);
        AddNodeCoverage(kmer_Bw);
        //add path coverage
        if (junctions.find(kmer_Fw)!=junctions.end()){
            junctions[kmer_Fw].right_kmers_cov[base_to_bit[ target_read[kmer_length] ]]++;
        }
        else if (junctions.find(kmer_Bw)!=junctions.end()){
            junctions[kmer_Bw].left_kmers_cov[base_to_bit[ trans_base[target_read[kmer_length]]]]++;
        }
        for (int i=kmer_length;i<target_read.size();i++){
            kmer_Fw=((kmer_Fw<<2)| end_bases[ base_to_bit[ target_read[i] ] + 4 ] );
            kmer_Bw=((kmer_Bw>>2)| end_bases[ base_to_bit[ trans_base[target_read[i]] ] ] );
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

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::AddNodeCoverage(BITLENGTH &target_kmer){
    BITLENGTH target_kmer_bw=GetComplementKmer(target_kmer);
    if (junctions.find(target_kmer)!=junctions.end()){
        junctions[target_kmer].coverage++;
    }
    if (joints.find(target_kmer)!=joints.end()) {
        joints[target_kmer].coverage++;
    }
}

template<typename BITLENGTH>
void DeBruijnGraph<BITLENGTH>::PrintGraph(){
    std::ofstream writing_gfa;
    std::string output_file ="./de_bruijn_graph.gfa";
    writing_gfa.open(output_file, std::ios::out);
    //header
    writing_gfa<<"H\tVN:Z:1.0"<<std::endl;
    //nodes(straight)
    for (auto itr=straights.begin();itr!=straights.end();++itr){
        writing_gfa<<"S"<<"\t"<<"Straight_"<<(itr->first)<<"\t"<<(itr->second).sequence<<"\t"<<"KC:i:"<<(itr->second).sequence.size()<<std::endl;
    }
    //nodes(junction)
    for (auto itr=junctions.begin();itr!=junctions.end();++itr){
        writing_gfa<<"S"<<"\t"<<"Junction_"<<(itr->second).id<<"\t"<<GetStringKmer((itr->first))<<"\t"<<"KC:i:"<<(itr->second).id*kmer_length<<std::endl;
    }
    //nodes(joint)
    for (auto itr=joints.begin();itr!=joints.end();++itr){
        writing_gfa<<"S"<<"\t"<<"Joint_"<<(itr->second).id<<"\t"<<GetStringKmer((itr->first))<<"\t"<<"KC:i:"<<(itr->second).id*kmer_length<<std::endl;
    }

    //edges(straight and joint) 
    for (auto itr=straights.begin();itr!=straights.end();++itr){
        writing_gfa<<"L"<<"\t";
        //left 
        writing_gfa<<"Joint_"<<(*((itr->second).left_joint)).id<<"\t"<<"+"<<"\t";
        writing_gfa<<"Straight_"<<(itr->first)<<"\t"<<"+"<<"\t";
        writing_gfa<<kmer_length<<"M"<<std::endl;
        //right
        writing_gfa<<"Straight_"<<(itr->first)<<"\t"<<"+"<<"\t";
        writing_gfa<<"Joint_"<<(*((itr->second).right_joint)).id<<"\t"<<"+"<<"\t";
        writing_gfa<<kmer_length<<"M"<<std::endl;
    }

    //edges(junction and joint)
    for (auto itr = junctions.begin(); itr!=junctions.end(); ++itr){
        writing_gfa<<"L"<<"\t";
        //left
        for (int i = 0; i < 4 ; i++){
            if ((itr->second).left_kmers_cov[i]==0) continue;
            BITLENGTH output_left_kmer=( ((itr->first)>>2) | end_bases[i] );
            if (junctions.find(output_left_kmer)!=junctions.end()){
                writing_gfa<<"Junction_"<<junctions[output_left_kmer].id<<"\t"<<"+"<<"\t";
                writing_gfa<<"Junction_"<<(itr->second).id<<"\t"<<"+"<<"\t";
                writing_gfa<<kmer_length-1<<"M"<<std::endl;
            }
            else {
                //joint
                if (joints.find(output_left_kmer)!=joints.end()){
                    writing_gfa<<"Joint_"<<joints[output_left_kmer].id<<"\t"<<"+"<<"\t";
                    writing_gfa<<"Junction_"<<(itr->second).id<<"\t"<<"+"<<"\t";
                    writing_gfa<<kmer_length-1<<"M"<<std::endl;
                    
                }
                //complement?
                else{
                    writing_gfa<<"Junction_"<<(itr->second).id<<"\t"<<"-"<<"\t";
                    writing_gfa<<"Junction_"<<(itr->second).id<<"\t"<<"+"<<"\t";
                    writing_gfa<<kmer_length-1<<"M"<<std::endl;
                }
            }
        }
        //right
        for (int i = 0; i < 4 ; i++){
            if ((itr->second).right_kmers_cov[i]==0) continue;
            BITLENGTH output_right_kmer=( ((itr->first)<<2) | end_bases[i+4] );
            /*right edge*/
            if (junctions.find(output_right_kmer)!=junctions.end()){
                writing_gfa<<"Junction_"<<(itr->second).id<<"\t"<<"+"<<"\t";
                writing_gfa<<"Junction_"<<junctions[output_right_kmer].id<<"\t"<<"+"<<"\t";
                writing_gfa<<kmer_length-1<<"M"<<std::endl;
            }
            else {
                //joint
                if (joints.find(output_right_kmer)!=joints.end()){
                    writing_gfa<<"Junction_"<<(itr->second).id<<"\t"<<"+"<<"\t";
                    writing_gfa<<"Joint_"<<joints[output_right_kmer].id<<"\t"<<"+"<<"\t";
                    writing_gfa<<kmer_length-1<<"M"<<std::endl;
                }
                //complement?
                else{
                    writing_gfa<<"Junction_"<<(itr->second).id<<"\t"<<"+"<<"\t";
                    writing_gfa<<"Junction_"<<(itr->second).id<<"\t"<<"-"<<"\t";
                    writing_gfa<<kmer_length-1<<"M"<<std::endl;
                }
            }
        }
    }
}
#endif