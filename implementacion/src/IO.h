#ifndef ITER_IO_H
#define ITER_IO_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>

#include <Eigen/Sparse>
#include "grafo.h"


namespace IO {

    /** UTILS */

    const int PRECISION = 15;

    size_t stolcast(const string &val, const string &msg);

    double stodcast(const string &val, const string &msg);

    map<string, string> oparams(int argc,  char** argv);

    string filename(const string& path);

    void skip_lines(ifstream &file, size_t n);


    /** FILE HANDLING */

    graph read_grafo(const string &in, size_t start=0);

    Eigen::SparseMatrix<double, Eigen::RowMajor> read_matriz(const string &in, size_t start=0);
    Eigen::SparseMatrix<double, Eigen::RowMajor> read_matriz(ifstream &file);

    Eigen::VectorXd read_vector(const string &in, size_t start=0);
    Eigen::VectorXd read_vector(ifstream &file);

    void write_matriz(const string &out, const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat, int precision=PRECISION);
    void write_matriz(ofstream &file, const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat, int precision=PRECISION);

    void write_vector(const string &out, const Eigen::VectorXd &vec, int precision=PRECISION);
    void write_vector(ofstream &file, const Eigen::VectorXd &vec, int precision=PRECISION);

    void write_time(const string &out, const vector<pair<string, chrono::microseconds>> &time_measures);

    pair<size_t, size_t> _shape(ifstream &in);
}


#endif //ITER_IO_H
