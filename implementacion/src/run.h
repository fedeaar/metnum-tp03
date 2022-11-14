#ifndef ITER_RUN_H
#define ITER_RUN_H

#include <Eigen/Sparse>
#include "IO.h"
#include "pagerank.h"


struct params {
    map<string, string> string_params {};
    map<string, double> double_params {};
    map<string, size_t> size_t_params {};
    map<string, bool>   bool_params {};
};

params get(int argc,  char** argv);

void run(const params &program);


#endif //ITER_RUN_H
