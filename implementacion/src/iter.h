#ifndef ITER_METNUM_H
#define ITER_METNUM_H

#include <Eigen/Sparse>


namespace metnum {

    void gauss_seidel(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b, double tol, size_t iter);

    void jacobi(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b, double tol, size_t iter);

    void eliminacion_gaussiana(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b);

    void backwards_substitution(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b);
}


#endif //ITER_METNUM_H
