#ifndef LOGGING_CPP
#define LOGGING_CPP
#include"common.h"



class Logging{
    public:
        std::mutex mtx_log;
        std::ofstream progress;
        std::string log_file ="./platanus3.log";
        Logging(){};
        void WriteLog(std::string text);
};

void Logging::WriteLog(std::string text){
    std::lock_guard<std::mutex> lock(mtx_log);
    progress.open(log_file, std::ios::out);
    progress<<text<<"\n";
}


#endif