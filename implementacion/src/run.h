#ifndef ITER_RUN_H
#define ITER_RUN_H

#include "IO.h"
#include "iter.h"
#include "pagerank.h"


struct params {
    std::map<std::string, std::string> string_params {};
    std::map<std::string, double> double_params {};
    std::map<std::string, size_t> size_t_params {};
    std::map<std::string, bool>   bool_params {};
};

params get(int argc,  char** argv);

void run(const params &program);


#endif //ITER_RUN_H
