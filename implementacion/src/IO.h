#ifndef ITER_IO_H
#define ITER_IO_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>

#include "grafo.h"
#include "iter.h"


namespace IO {

    /** UTILS */

    const int PRECISION = 15;

    size_t stolcast(const std::string &val, const std::string &msg);

    double stodcast(const std::string &val, const std::string &msg);

    std::map<std::string, std::string> oparams(int argc,  char** argv);

    std::string filename(const std::string& path);

    void skip_lines(std::ifstream &file, size_t n);


    /** FILE HANDLING */

    graph read_grafo(const std::string &in, size_t start=0);

    metnum::RowMatrix read_matriz(const std::string &in, size_t start=0);
    metnum::RowMatrix read_matriz(std::ifstream &file);

    metnum::DenseVector read_vector(const std::string &in, size_t start=0);
    metnum::DenseVector read_vector(std::ifstream &file);

    void write_matriz(const std::string &out, const metnum::RowMatrix &mat, int precision=PRECISION);
    void write_matriz(std::ofstream &file, const metnum::RowMatrix &mat, int precision=PRECISION);

    void write_vector(const std::string &out, const metnum::DenseVector &vec, int precision=PRECISION);
    void write_vector(std::ofstream &file, const metnum::DenseVector &vec, int precision=PRECISION);

    void write_time(const std::string &out,
                    const std::vector<std::pair<std::string, std::chrono::microseconds>> &time_measures);

    std::pair<size_t, size_t> _shape(std::ifstream &in);
}


#endif //ITER_IO_H
