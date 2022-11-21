#include <iostream>
#include "../iter.h"
using namespace metnum;

//
// UTILS
//

DenseVector aleatorio(size_t n, std::pair<int, int> range={INT32_MIN, INT32_MAX}) {
    DenseVector res(n);
    std::random_device rd;
    std::mt19937 rng {rd()}; // Mersenne Twister
    std::uniform_int_distribution<int> dist(range.first, range.second);
    bool cero = true;
    for (double & re : res) {
        re = dist(rng);
        cero = cero && re == 0;
    }
    if (cero) { // si res == 0, usamos e1.
        res[0] = 1;
    }
    return res;
}


DenseVector normalizar(const DenseVector &v) {
    double n = sqrt(v.dot(v));
    return abs(n) < EPSILON ? v : v / n;
}


llong min(llong a, llong b) {
    return a <= b ? a : b;
}




//
// METODOS
//

DenseVector metnum::gauss_seidel(RowMatrix &A, DenseVector &b, double tol, size_t iter) {
    DenseVector x(b.size()), y(b.size()), z(b.size());
    x = normalizar(aleatorio(A.size()));
    bool converge = false;
    for (int i = 0; i < iter && !converge; ++i) {
        // Actualizo x
        y = x;
        for (int j = 0; j < x.size(); ++j) {
            double suma = 0;
            for (RowIterator it(A[j]); it; ++it) {
                // Si it.col() < j, va a tomar los valores nuevos de x
                if (it.row() != j) {
                    suma += it.value() * x[it.row()];
                }
            }
            x[j] = (b[j] - suma) / A[j].coeff(j);
        }
        // Chequeo convergencia
        z = x - y;
        converge = sqrt(z.dot(z)) < tol;
    }
    return x;
}


DenseVector metnum::jacobi(RowMatrix &A, DenseVector &b, double tol, size_t iter) {
    DenseVector x(b.size()), y(b.size()), z(b.size());
    x = normalizar(aleatorio(A.size()));
    bool converge = false;
    for (int i = 0; i < iter && !converge; ++i) {
        // Actualizo x
        y = x;
        for (int j = 0; j < x.size(); ++j) {
            double suma = 0;
            for (RowIterator it(A[j]); it; ++it) {
                // Recorro solo los indices no nulos de la fila j
                if (it.row() != j) {
                    suma += it.value() * y[it.row()];
                }
            }
            x[j] = (b[j] - suma) / A[j].coeff(j);
        }
        // Chequeo convergencia
        z = x - y;
        converge = sqrt(z.dot(z)) < tol;
    }
    return x;
}


void metnum::eliminacion_gaussiana(RowMatrix &A, DenseVector &b, double tol) {
    // pre: A_ii != 0 para i: 0 ... N. hasta el final de la eliminación
    //      b.size() == A.cols() == A.rows()
    llong n = A.size();
    for (int i = 0; i < n-1; ++i) {
        double mii = A[i].coeff(i);
        for (int j = i+1; j < n; ++j) {
            double mij = A[j].coeff(i) / mii;
            if (abs(mij) < tol) continue;
            A[j] = (A[j] - A[i] * mij).pruned(1, tol);
            b[j] = b[j] - b[i] * mij;
        }
    }
}


DenseVector metnum::backwards_substitution(RowMatrix &A, DenseVector &b) {
    // pre: A es triangular superior cuadrada
    //      b.size() == a.rows()
    DenseVector res = b;
    for (int i = A.size() - 1; i >= 0; --i) {
        double parcial = 0;
        for (RowIterator it(A[i]); it; ++it) {
            if (it.row() < i + 1) continue;
            parcial += it.value() * res[it.row()];
        }
        res[i] = (res[i] - parcial) / A[i].coeff(i);
    }
    return res;
}
