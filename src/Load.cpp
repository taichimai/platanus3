#ifndef LOAD_CPP
#define LOAD_CPP
#include"common.h"

class ReadFile{
    public:
        ReadSet reads;
        std::string file_name;
        std::string file_type;
        Error error_code;
        uint64_t all_bases=0;

        ReadFile(std::string input_file_name);
        void LoadFile();
        void LoadFasta(ReadSet *loaded_reads,std::string file_name);
        void LoadFastq(ReadSet *loaded_reads,std::string file_name);
        KmerSet GetSeedKmer(int kmer_length);
};

ReadFile::ReadFile(std::string input_file_name){
    file_name=input_file_name;
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
                (*loaded_reads)[read_name]=seq;
                all_bases+=seq.size();
                seq="";
            }
            read_name=line;
        }
        else{
            seq+=line;
        }
    }
    (*loaded_reads)[read_name]=seq;
    all_bases+=seq.size();
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
                (*loaded_reads)[read_name]=seq;
                all_bases+=seq.size();
                seq="";
            }
            read_name=line;
        }
        else if (line_cnt%4==1){
            seq+=line;
        }
        line_cnt+=1;
    }
    (*loaded_reads)[read_name]=seq;
    all_bases+=seq.size();
}


KmerSet ReadFile::GetSeedKmer(int kmer_length){
    KmerSet edge_kmers;
    for (auto itr=reads.begin();itr!=reads.end();++itr){
        edge_kmers.insert((itr->second).substr(0,kmer_length));
    }
    return edge_kmers;
}
#endif 