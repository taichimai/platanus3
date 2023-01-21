#ifndef LOAD_CPP
#define LOAD_CPP
#include"common.h"

class ReadFile{
    public:
        ReadSet reads;
        std::string file_name;
        std::string file_type;
        Error error_code;

        ReadFile(std::string input_file_name);
        void LoadFile();
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
    ReadSet  input_reads;
    std::ifstream lines(file_name);
    std::string line;
    std::string seq;
    std::string read_name;
    int line_cnt=0;
    char name_symbol;
    int  part_num;
    //fasta
    if (file_type=="fasta"){ 
      name_symbol='>';
      part_num=2;
    }
    //fastq
    if (file_type=="fastq"){
      name_symbol='@';
      part_num=4;
    }

    while (getline(lines,line)) {
        if (line_cnt%part_num==0){
            if (line[0]!=name_symbol) {
                error_code=2;
                return;
            }
            if (read_name!=""){
                input_reads[read_name]=seq;
                seq="";
            }
            read_name=line;
        }
        else if (line_cnt%part_num==1){
            seq+=line;
        }
        line_cnt+=1;
    }
    input_reads[read_name]=seq;
    reads=input_reads;
}



KmerSet ReadFile::GetSeedKmer(int kmer_length){
    KmerSet edge_kmers;
    for (auto itr=reads.begin();itr!=reads.end();++itr){
        edge_kmers.insert((itr->second).substr(0,kmer_length));
    }
    return edge_kmers;
}
#endif 