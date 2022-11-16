#ifndef ITER_PAGERANK_H
#define ITER_PAGERANK_H

#include <vector>
#include "grafo.h"
#include "IO.h"
#include "iter.h"


namespace pagerank {

    enum metodo {EG, GS, J};

    namespace IO {

        struct in_file {
            in_file(graph g, double _val) : p_val(_val), grafo(std::move(g)) {}

            double p_val;
            graph grafo;
        };
        struct out_file {
            out_file(metnum::DenseVector v, double _val): p_val(_val), solucion(std::move(v)) {}

            double p_val;
            metnum::DenseVector solucion;
        };

        in_file read_in(const std::string &in, double p_val);

        out_file read_out(const std::string &in);

        void write_out(const std::string &out, const out_file &data, int precision=::IO::PRECISION);
    };

    metnum::RowMatrix make(const IO::in_file &params);

    metnum::DenseVector solve(metnum::RowMatrix &mat, metodo met=EG, double tol=1e-4, size_t iter=1e5);
}


#endif //ITER_PAGERANK_H
