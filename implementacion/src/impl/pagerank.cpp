#include "../pagerank.h"
using namespace std;


//
// IO
//

pagerank::IO::in_file pagerank::IO::read_in(const string &in, double p_val) {
    return {::IO::read_grafo(in), p_val};
}


pagerank::IO::out_file pagerank::IO::read_out(const string &in) {
    ifstream file{in, ios::binary};
    if (!file.is_open()) {
        throw std::invalid_argument("no se pudo leer el archivo: " + in + ".");
    }
    string _p_val;
    std::getline(file, _p_val);
    string msg = "error en lectura del valor p en el archivo: " + in + " .\n";
    double p_val = ::IO::stodcast(_p_val, msg);
    return {::IO::read_vector(file), p_val};
}


void pagerank::IO::write_out(const string &out, const out_file &data, int precision) {
    ofstream file{out, ios::binary};
    file << data.p_val << "\n";
    ::IO::write_vector(file, data.solucion, precision);
}




//
// PAGERANK
//

metnum::RowMatrix pagerank::make(const pagerank::IO::in_file &params) {
    size_t n = params.grafo.nodos;
    vector<double> grado(n);
    for (auto &r: params.grafo.relaciones) {
        ++grado[r.j - 1];
    }
    std::vector<Eigen::Triplet<double>> t_list;
    t_list.reserve(params.grafo.links + n);
    for (auto &r: params.grafo.relaciones) {
        if (r.i != r.j) {
            Eigen::Triplet<double> pos {r.i - 1, r.j - 1, -params.p_val / grado[r.j - 1]};
            t_list.emplace_back(pos);
        }
    }
    for (long long i = 0; i < n; ++i) {
        t_list.emplace_back(Eigen::Triplet<double>{i, i, 1});
    }
    // m = I - pWD
    Eigen::SparseMatrix<double, Eigen::RowMajor> m{(Eigen::Index) n, (Eigen::Index) n};
    m.setFromTriplets(t_list.begin(), t_list.end());
    // convierto
    std::vector<Eigen::SparseVector<double>> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = m.row(i);
    }
    return res;
}


metnum::DenseVector pagerank::solve(metnum::RowMatrix &mat, pagerank::metodo met, double tol, size_t iter) {
    metnum::DenseVector b = metnum::DenseVector::Ones(mat.size());
    metnum::DenseVector res;
    switch (met) {
        case pagerank::EG:
            metnum::eliminacion_gaussiana(mat, b, tol);
            res = metnum::backwards_substitution(mat, b);
            break;
        case pagerank::GS:
            res = metnum::gauss_seidel(mat, b, tol, iter);
            break;
        case pagerank::J:
            res = metnum::jacobi(mat, b, tol, iter);
            break;
    }
    res = res / res.sum();
    return res;
}
