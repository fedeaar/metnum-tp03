#ifndef ITER_METNUM_H
#define ITER_METNUM_H

#include <Eigen/Sparse>
#include <random>


namespace metnum {

    typedef long long llong;
    typedef std::vector<Eigen::SparseVector<double>> RowMatrix;
    typedef Eigen::SparseVector<double>::InnerIterator RowIterator;
    typedef Eigen::VectorXd DenseVector;
    typedef Eigen::SparseVector<double> SparseVector;

    const double EPSILON = 1e-6;


    DenseVector gauss_seidel(RowMatrix &A, DenseVector &b, double tol, size_t iter);

    DenseVector jacobi(RowMatrix &A, DenseVector &b, double tol, size_t iter);

    void eliminacion_gaussiana(RowMatrix &A, DenseVector &b, double tol=EPSILON);

    DenseVector backwards_substitution(RowMatrix &A, DenseVector &b);
}


#endif //ITER_METNUM_H
