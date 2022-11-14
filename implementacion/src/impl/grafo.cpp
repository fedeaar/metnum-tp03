#include "../grafo.h"


//
// GRAFO
//

Eigen::SparseMatrix<double> grafo_a_matriz(const grafo &g) {
    std::vector<Eigen::Triplet<double>> t_list;
    t_list.reserve(g.links);
    for (auto &r: g.relaciones) {
        if (r.i != r.j) {
            t_list.emplace_back(Eigen::Triplet<double>(r.i - 1, r.j - 1, 1));
        }
    }
    Eigen::SparseMatrix<double> res(g.nodos, g.nodos);
    res.setFromTriplets(t_list.begin(), t_list.end());
    return res;
}
