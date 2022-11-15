#include <iostream>
#include "../iter.h"

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> RowMatrix;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator RowIterator;
typedef Eigen::VectorXd EigenVector;


//
// UTILS
//

void printM(Eigen::SparseMatrix<double> &A){
    std::cout << A << std::endl;
}
void printM(RowMatrix &A){
    std::cout << A << std::endl;
}




//
// METODOS ITERATIVOS
//

std::pair<bool, EigenVector> metnum::gauss_seidel(RowMatrix &A, EigenVector &b, double tol, size_t iter) {

    EigenVector x(b.size()), y(b.size()), z(b.size());
    bool converge = false;

    for (int i = 0; i < iter && !converge; ++i) {
        // Actualizo x
        y = x;
        for (int j = 0; j < x.size(); ++j) {
            double suma = 0;
            for (RowIterator it(A, j); it; ++it) {
                // Si it.col() < j, va a tomar los valores nuevos de x
                if (it.col() != j) {
                    suma += it.value() * x[it.col()];
                }
            }
            x[j] = (b[j] - suma) / A.coeff(j, j);
        }
        // Chequeo convergencia
        z = x - y;
        converge = sqrt(z.dot(z)) < tol;
    }

    return std::make_pair(converge, x);
}


std::pair<bool, EigenVector> metnum::jacobi(RowMatrix &A, EigenVector &b, double tol, size_t iter) {

    EigenVector x(b.size()), y(b.size()), z(b.size());
    bool converge = false;

    for (int i = 0; i < iter && !converge; ++i) {
        // Actualizo x
        y = x;
        for (int j = 0; j < x.size(); ++j) {
            double suma = 0;
            for (RowIterator it(A, j); it; ++it) {
                // Recorro solo los indices no nulos de la fila j
                if (it.col() != j) {
                    suma += it.value() * y[it.col()];
                }
            }
            x[j] = (b[j] - suma) / A.coeff(j, j);
        }
        // Chequeo convergencia
        z = x - y;
        converge = sqrt(z.dot(z)) < tol;
    }

    return std::make_pair(converge, x);
}




//
// ELIMINACION GAUSSIANA
//

void metnum::eliminacion_gaussiana(RowMatrix &A, EigenVector &b) {
    // pre: A_ii != 0 para i: 0 ... N. hasta el final de la eliminación
    //      b.size() == A.cols() == A.rows()
    int n = A.rows();

    for (int i = 0; i < n-1; ++i) {
        double mii = A.coeff(i, i);
        for(int j = i+1; j < n; ++j){
            double mij = A.coeff(j, i) / mii;
            if(mij == 0) continue;
            Eigen::SparseVector<double> new_Bj(n);
            new_Bj.reserve(A.row(j).nonZeros() + A.row(i).nonZeros() - 2);

            RowIterator it_fila_i(A, i);
            RowIterator it_fila_j(A, j);
            ++it_fila_j; // ya que el primer elemento sabemos que se convierte en 0
            ++it_fila_i; // por lo que salteamos el primer elemento
            while(it_fila_i || it_fila_j){
                if(it_fila_j && (!it_fila_i || it_fila_j.col() < it_fila_i.col())){
                    // si it_fila_i > it_fila_j agrega el elem que no va a editar de la fila original
                    new_Bj.insertBack(it_fila_j.col()) = it_fila_j.value();
                    ++it_fila_j;
                } else {
                    double newVal = - it_fila_i.value() * mij;
                    if(it_fila_j && it_fila_j.col() == it_fila_i.col()) {
                        // si existe el elem it_fila_j entonces lo suma ya que sino es 0 y no afecta
                        newVal += it_fila_j.value();
                        ++it_fila_j;
                    }
                    if(abs(newVal) > 1e-4) new_Bj.insertBack(it_fila_i.col()) = newVal;
                    ++it_fila_i;
                }
            }

            A.row(j) = new_Bj;
            b[j] = b[j] - b[i] * mij;
        }
    }
}


std::pair<bool, EigenVector> metnum::backwards_substitution(RowMatrix &A, EigenVector &b) {
    // pre: A es triangular superior cuadrada
    //      b.size() == a.rows()
    EigenVector res = b;
    for (int i = A.rows() - 1; i >= 0; --i) {
        double parcial = 0;
        for (RowIterator it(A, i); it; ++it) {
            if (it.col() < i + 1) continue;
            parcial += it.value() * res[it.col()];
        }
        res[i] = (res[i] - parcial) / A.coeff(i, i);
    }
    return {true, res};
}
