#include "../grafo.h"


//
// GRAFO
//

metnum::RowMatrix grafo_a_matriz(const graph &g) {
    std::vector<Eigen::Triplet<double>> t_list;
    t_list.reserve(g.links);
    for (auto &r: g.relaciones) {
        if (r.i != r.j) {
            t_list.emplace_back(Eigen::Triplet<double>(r.i - 1, r.j - 1, 1));
        }
    }
    Eigen::SparseMatrix<double, Eigen::RowMajor> m(g.nodos, g.nodos);
    m.setFromTriplets(t_list.begin(), t_list.end());
    std::vector<Eigen::SparseVector<double>> res(g.nodos);
    for (int i = 0; i < g.nodos; ++i) {
        res[i] = m.row(i);
    }
    return res;
}
