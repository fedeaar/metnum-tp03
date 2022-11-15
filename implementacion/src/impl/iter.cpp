#include <iostream>
#include "../iter.h"

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



void metnum::eliminacion_gaussiana(RowMatrix &A, EigenVector &b, double tol) {
   // pre: A_ii != 0 para i: 0 ... N. hasta el final de la eliminación
   //      b.size() == A.cols() == A.rows()

    int n = A.rows();
    std::vector<Eigen::SparseVector<double>> B(n);
    for(int i = 0; i < n; ++i) B[i] = A.row(i);

    for (int i = 0; i < n-1; ++i) {
        double mii = B[i].coeff(i);
        for(int j = i+1; j < n; ++j){
            double mij = B[j].coeff(i) / mii;
            if(abs(mij) < tol) continue;

            Eigen::SparseVector<double> new_Bj(n);
            new_Bj.reserve(fmin(B[j].nonZeros() + B[i].nonZeros() - 1, n));

            Eigen::SparseVector<double>::InnerIterator it_fila_i(B[i]);
            Eigen::SparseVector<double>::InnerIterator it_fila_j(B[j]);
            ++it_fila_j; // ya que el primer elemento sabemos que se convierte en 0
            ++it_fila_i; // por lo que salteamos el primer elemento
            while(it_fila_i || it_fila_j) {
                if (!it_fila_i || (it_fila_j && it_fila_j.index() < it_fila_i.index())) {
                    // si it fila_i > it_fila_j agrega el elem que no va a editar de la fila original
                    new_Bj.insertBack(it_fila_j.index()) = it_fila_j.value();
                    ++it_fila_j;
                } else {
                    // si entramos aca si o si se cumple it_fila_i valido
                    double newVal = -it_fila_i.value() * mij;
                    if (it_fila_j && it_fila_j.index() == it_fila_i.index()) {
                        // si existe el elem it_fila_j entonces lo suma ya que sino es 0 y no afecta
                        newVal += it_fila_j.value();
                        ++it_fila_j;
                    }
                    if (abs(newVal) > tol) new_Bj.insertBack(it_fila_i.index()) = newVal;
                    ++it_fila_i;
                }
            }

            B[j] = new_Bj;
            b[j] = b[j] - b[i] * mij;
        }
    }
    for(int i = 0; i < n; ++i) A.row(i) = B[i];
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
