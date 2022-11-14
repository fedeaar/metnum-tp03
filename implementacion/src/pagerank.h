#ifndef ITER_PAGERANK_H
#define ITER_PAGERANK_H

#include <Eigen/Sparse>
#include <utility>
#include <vector>
#include "grafo.h"
#include "IO.h"
#include "iter.h"


namespace pagerank {

    enum metodo {EG, GS, J};

    namespace IO {

        struct in_file {
            in_file(grafo g, double _val) : p_val(_val), grafo(std::move(g)) {}

            double p_val;
            grafo grafo;
        };

        struct out_file {
            out_file(Eigen::VectorXd v, double _val): p_val(_val), solucion(std::move(v)) {}

            double p_val;
            Eigen::VectorXd solucion;
        };

        in_file read_in(const string &in, double p_val);

        out_file read_out(const string &in);

        void write_out(const string &out, const out_file &data, int precision=::IO::PRECISION);
    };

    Eigen::SparseMatrix<double> make(const IO::in_file &params);

    Eigen::VectorXd solve(Eigen::SparseMatrix<double>& mat, metodo met=EG, double tol=1e-4, size_t iter=1e5);
};


#endif //ITER_PAGERANK_H
