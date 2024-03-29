#ifndef ITER_GRAFO_H
#define ITER_GRAFO_H

#include <cstdlib>
#include <vector>
#include "iter.h"


struct coords {
    coords(size_t _i, size_t _j) : i(_i), j(_j) {}

    size_t i;
    size_t j;
};
struct graph {
    graph(size_t _n, size_t _l): nodos(_n), links(_l), relaciones() {}

    size_t nodos;
    size_t links;
    std::vector<coords> relaciones;
};

metnum::RowMatrix grafo_a_matriz(const graph &g);


#endif //ITER_GRAFO_H
