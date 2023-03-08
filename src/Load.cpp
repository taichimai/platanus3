#ifndef LOAD_CPP
#define LOAD_CPP
#include"common.h"
#include"BitCalc.cpp"

class ReadFile{
    public:
        ReadSet reads;
        std::string file_name;
        uint32_t large_kmer_length;
        std::string file_type;
        Error error_code;
        uint64_t all_bases=0;
        KmerCount shortk_database;

        ReadFile(Options &parameters);
        void LoadFile();
        void LoadFasta(ReadSet *loaded_reads,std::string file_name);
        void LoadFastq(ReadSet *loaded_reads,std::string file_name);
        void CountShortKmer(int shortk_length);
};

ReadFile::ReadFile(Options &parameters){
    file_name=parameters.readfile_name;
    large_kmer_length=parameters.kmer_length;
    file_type=file_name.substr(file_name.size()-5,5);
    if (not (file_type=="fasta" or file_type=="fastq")){
        error_code=1;
    }
}

void ReadFile::LoadFile(){
    std::ifstream lines(file_name);
    std::string first_line;
    int is_name=0;
    char name_symbol;
    int  part_num;
    int line_cnt=0;

    getline(lines,first_line);
    if (first_line[0]=='>'){
        lines.close();
        LoadFasta(&reads,file_name);
    }
    else if (first_line[0]=='@') {
        lines.close();
        LoadFastq(&reads,file_name);
    }
}

void ReadFile::LoadFasta(ReadSet *loaded_reads,std::string file_name){
    std::ifstream fasta_lines(file_name);
    std::string line="";
    std::string seq="";
    std::string read_name="";
    while (getline(fasta_lines,line)){
        if (line[0]=='>') {
            if (read_name!=""){
                if (seq.size()>=large_kmer_length){
                    (*loaded_reads)[read_name]=seq;
                    all_bases+=seq.size();
                }
                seq="";
            }
            read_name=line;
        }
        else{
            seq+=line;
        }
    }
    if (seq.size()>=large_kmer_length){
        (*loaded_reads)[read_name]=seq;
        all_bases+=seq.size();
    }
}

void ReadFile::LoadFastq(ReadSet *loaded_reads,std::string file_name){
    std::ifstream fastq_lines(file_name);
    std::string line="";
    std::string seq="";
    std::string read_name="";
    int line_cnt=0;
    while (getline(fastq_lines,line)) {
        if (line_cnt%4==0){
            if (read_name!=""){
                if (seq.size()>=large_kmer_length){
                    (*loaded_reads)[read_name]=seq;
                    all_bases+=seq.size();
                }
                seq="";
            }
            read_name=line;
        }
        else if (line_cnt%4==1){
            seq+=line;
        }
        line_cnt+=1;
    }
    if (seq.size()>=large_kmer_length){
        (*loaded_reads)[read_name]=seq;
        all_bases+=seq.size();
    }
}

void ReadFile::CountShortKmer(int shortk_length){
    std::bitset<42>  A_right_short(0); std::bitset<42> A_left_short=(A_right_short<<(shortk_length*2-2)); 
    std::bitset<42>  C_right_short(1); std::bitset<42> C_left_short=(C_right_short<<(shortk_length*2-2));
    std::bitset<42>  G_right_short(2); std::bitset<42> G_left_short=(G_right_short<<(shortk_length*2-2));
    std::bitset<42>  T_right_short(3); std::bitset<42> T_left_short=(T_right_short<<(shortk_length*2-2));
    std::vector<std::bitset<42>> end_bases_short={A_left_short,C_left_short,G_left_short,T_left_short,A_right_short,C_right_short,G_right_short,T_right_short};
    
    KmerCount all_short_kmers;
    for (auto itr = reads.begin(); itr!=reads.end();++itr){
      std::string target_read = (itr->second); 
      std::bitset<42> shortk_for=GetFirstKmerForward<std::bitset<42>>(target_read.substr(0,shortk_length));
      std::bitset<42> shortk_rev=GetFirstKmerBackward<std::bitset<42>>(target_read.substr(0,shortk_length));

      for (int i=shortk_length-1;i<target_read.size();i++){
        if (i!=shortk_length-1){
          shortk_for=( (shortk_for<<2) | end_bases_short[ base_to_bit[ target_read[i] ] + 4 ]   );
          shortk_rev=( (shortk_rev>>2) | end_bases_short[ base_to_bit[ trans_base[target_read[i]] ] ] );
        }
        all_short_kmers[CompareBit(shortk_for,shortk_rev,42)]++;
      }
    }
    this->shortk_database=all_short_kmers;
}
#endif
