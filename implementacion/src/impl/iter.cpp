#include <iostream>
#include "../iter.h"


//
// METODOS ITERATIVOS
//

void metnum::gauss_seidel(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b, double tol, size_t iter) {

    // TODO
}


void metnum::jacobi(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b, double tol, size_t iter) {

    // TODO
}


/*
void metnum::eliminacion_gaussiana(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b) {
    // pre: A_ii != 0 para i: 0 ... N. hasta el final de la eliminación
    //      b.size() == A.cols() == A.rows()
    size_t n = A.cols();
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(A.rows(), A.cols());

    std::vector<Eigen::Triplet<double>> t_list;
    t_list.reserve(A.rows() * 2);
    for (int k = 0; k < n; ++k) {
        t_list.emplace_back(Eigen::Triplet<double>{k, k, 1});
    }
    size_t pos = n;

    for (int i = 0; i < A.rows() - 1; ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
            // M
            int k = it.row();
            if (k <= i) {
                continue;
            }
            double mij = - A.coeff(k, i) / A.coeff(i, i);
            t_list[pos] = Eigen::Triplet<double>{k, i, mij};
            pos++;

            b[k] = b[k] + b[i] * mij;
        }
        M.setFromTriplets(t_list.begin(), t_list.begin() + pos);
        pos = n;

        A = (M * A).pruned(1e-4);
    }
}*/



void metnum::eliminacion_gaussiana(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b) {
    // pre: A_ii != 0 para i: 0 ... N. hasta el final de la eliminación
    //      b.size() == A.cols() == A.rows()
    Eigen::SparseMatrix<double, Eigen::RowMajor> B = A;

    std::vector<int> p_list;
    p_list.reserve(A.rows());
    int p_pos = 0;

    for (int i = 0; i < A.rows() - 1; ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
            // M
            int k = it.row();
            if (k <= i) {
                continue;
            }
            p_list[p_pos] = k;
            p_pos++;
        }
        while(p_pos > 0) {
            int j = p_list[p_pos - 1];
            double mij = A.coeff(j, i) / A.coeff(i, i);
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(B, i); it; ++it) {
                int k = it.col();
                if (k < i) {
                    continue;
                }
                double val = A.coeff(j, k) - A.coeff(i, k) * mij;
                A.coeffRef(j, k) = val;
            }
            b[j] = b[j] - b[i] * mij;
            --p_pos;
        }
        // B = A;
    }
}


void metnum::backwards_substitution(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b) {
    // pre: A es triangular superior cuadrada
    //      b.size() == a.rows()
    Eigen::SparseMatrix<double, Eigen::RowMajor> B = A;
    for (int i = B.rows() - 1; i >= 0; --i) {
        double parcial = 0;
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(B, i); it; ++it) {
            if (it.col() < i + 1) continue;
            parcial += it.value() * b[it.col()];
        }
        b[i] = (b[i] - parcial) / B.coeff(i, i);
    }
}
