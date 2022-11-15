#include <iostream>
#include "../iter.h"


//
// METODOS ITERATIVOS
//

std::pair<bool, Eigen::VectorXd> metnum::gauss_seidel(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b, double tol, size_t iter) {
    Eigen::VectorXd x(b.size()), y(b.size()), z(b.size());
    bool converge = false;

    for (int i = 0; i < iter && !converge; ++i) {
        // Actualizo x
        y = x;
        for (int j = 0; j < x.size(); ++j) {
            int suma = 0;
            for (SparseMatrix<double>::InnerIterator it(mat,j); it; ++it) { 
                // Si it.col() < j, va a tomar los valores nuevos de x
                if (it.col() != j) {
                    suma += it.value() * x[it.col()];
                }
            }
            x[j] = (b[j] - sum) / A[j,j];
        }
        // Chequeo convergencia
        z = x - y;
        converge = sqrt(z.dot(z)) < tol;
    }

    return std::make_pair(converge, x);
}


std::pair<bool, Eigen::VectorXd> metnum::jacobi(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b, double tol, size_t iter) {

    Eigen::VectorXd x(b.size()), y(b.size()), z(b.size());
    bool converge = false;

    for (int i = 0; i < iter && !converge; ++i) {
        // Actualizo x
        y = x;
        for (int j = 0; j < x.size(); ++j) {
            int suma = 0;
            for (SparseMatrix<double>::InnerIterator it(mat,j); it; ++it) { 
                // Recorro solo los indices no nulos de la fila j
                if (it.col() != j) {
                    suma += it.value() * y[it.col()];
                }
            }
            x[j] = (b[j] - sum) / A[j,j];
        }
        // Chequeo convergencia
        z = x - y;
        converge = sqrt(z.dot(z)) < tol;
    }

    return std::make_pair(converge, x);
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



// void metnum::eliminacion_gaussiana(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &b) {
//     // pre: A_ii != 0 para i: 0 ... N. hasta el final de la eliminación
//     //      b.size() == A.cols() == A.rows()
//     Eigen::SparseMatrix<double, Eigen::RowMajor> B = A;

//     std::vector<int> p_list;
//     p_list.reserve(A.rows());
//     int p_pos = 0;

//     for (int i = 0; i < A.rows() - 1; ++i) {
//         for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
//             // M
//             int k = it.row();
//             if (k <= i) {
//                 continue;
//             }
//             p_list[p_pos] = k;
//             p_pos++;
//         }
//         while(p_pos > 0) {
//             int j = p_list[p_pos - 1];
//             double mij = A.coeff(j, i) / A.coeff(i, i);
//             for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(B, i); it; ++it) {
//                 int k = it.col();
//                 if (k < i) {
//                     continue;
//                 }
//                 double val = A.coeff(j, k) - A.coeff(i, k) * mij;
//                 A.coeffRef(j, k) = val;
//             }
//             b[j] = b[j] - b[i] * mij;
//             --p_pos;
//         }
//         // B = A;
//     }
// }

void printM(Eigen::SparseMatrix<double> &A){
    std::cout << A << std::endl;
}
void printM(Eigen::SparseMatrix<double, Eigen::RowMajor> &A){
    std::cout << A << std::endl;
}
// void metnum::eliminacion_gaussiana(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b) {
//     // pre: A_ii != 0 para i: 0 ... N. hasta el final de la eliminación
//     //      b.size() == A.cols() == A.rows()
//
//     Eigen::SparseMatrix<double> B = A;
//
//     int n = B.rows();
//
//     for(int i = 0; i < n; ++i){
//         std::vector<Eigen::Triplet<double>> M_list;
//         for(int j = 0; j < n; ++j)
//             M_list.emplace_back(Eigen::Triplet<double> {j, j, 1});
//
//
//         Eigen::SparseMatrix<double>::InnerIterator it(B, i);
//         while(it.row() < i) ++it;
//         double mii = it.value();
//         ++it;
//         for ( ; it ; ++it)
//             M_list.emplace_back(Eigen::Triplet<double> {it.row(), i, -it.value() / mii});
//
//         Eigen::SparseMatrix<double> M{(Eigen::Index) n, (Eigen::Index) n};
//         M.setFromTriplets(M_list.begin(), M_list.end());
//         B = (M * B).pruned(1e-3);
//         b = M * b;
//     }
//
//     A = B;
// }
//void metnum::eliminacion_gaussiana(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b) {
//    // pre: A_ii != 0 para i: 0 ... N. hasta el final de la eliminación
//    //      b.size() == A.cols() == A.rows()
//
//    for (int i = 0; i < A.rows() - 1; ++i) {
//        printM(A);
//        double mii = A.coeff(i, i);
//        for(int j = i+1; j < A.rows(); ++j){
//            double mij = A.coeff(j, i) / mii;
//            if(mij == 0) continue;
//            for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_fila_i(A, i); it_fila_i; ++it_fila_i){
//                A.coeffRef(j, it_fila_i.col()) -= it_fila_i.value() * mij;
//            }
//            b[j] = b[j] - b[i] * mij;
//        }
//        // A.makeCompressed();
//    }
//
//}
void metnum::eliminacion_gaussiana(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b) {
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

            Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_fila_i(A, i);
            Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_fila_j(A, j);
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





void metnum::backwards_substitution(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &b) {
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
