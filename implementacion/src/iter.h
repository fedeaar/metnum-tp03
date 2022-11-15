#ifndef ITER_METNUM_H
#define ITER_METNUM_H

#include <Eigen/Sparse>



namespace metnum {

    std::pair<bool, Eigen::VectorXd>
    gauss_seidel(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b, double tol, size_t iter);

    std::pair<bool, Eigen::VectorXd>
    jacobi(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b, double tol, size_t iter);

    void eliminacion_gaussiana(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b);

    std::pair<bool, Eigen::VectorXd>
    backwards_substitution(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b);
}


#endif //ITER_METNUM_H
