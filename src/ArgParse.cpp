#include<common.h>

bool Parse(int argc,char **argv,std::string*readfile_name){
    int opt;
    const char* options="i:";
    while((opt = getopt(argc, argv, options)) != -1){
        switch(opt)
        {
            case 'i':
                *readfile_name = optarg;
                break;
            default: 
                std::cerr << "Invalid option" << std::endl;
                return false;
                break;
        }
    }

    return true;


}