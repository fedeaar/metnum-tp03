import base.IO as IO
import base.utils as utils

import os
import numpy as np
import pandas as pd


"""
descripcion: 
    evaluación de la convergencia de las soluciones de PageRank para una serie de 
    casos de test. Se compararon dos implementaciones del algoritmo -sobre los 
    métodos iterativos de Gauss-Seidel y Jacobi- y se midió el error absoluto de 
    los resultados en función de la cantidad de iteraciones realizadas.

    Nota: correr el experimento regenera todos los archivos. Esto puede resultar 
    en pequeñas discrepancias con los resultados provistos.
"""


#
# IO
#

TESTS_DIR = "../catedra/tests-pagerank/"

EXPERIMENTO = "convergencia-iterativos"
DIR_IN, DIR_OUT, DIR = IO.createInOut(EXPERIMENTO, delete=True)

RESULTADOS  = DIR + "error-convergencia.csv"
COLS        = 'test,metodo,n_iter,error_abs,error_rel,error_max'
FMT_COLS    = "{0},{1},{2},{3},{4},{5}\n"

GRAFICOS    = DIR + "convergencia_{test}.png"


#
# VAR
#

TESTS   = [x[:-4] for x in os.listdir(TESTS_DIR) if \
           os.path.isfile(TESTS_DIR + x) and x[-4:] == ".txt"]
METODOS = ['GS', 'J']
ITER    = int(1e2)
STEP    = 1


#
# UTILS
#

def correr_pagerank():

    for test in TESTS: 

        for metodo in METODOS:

            print(f"corriendo test: {test}, metodo: {metodo}.")

            for k in range(0, ITER, STEP):

                in_file = TESTS_DIR + test + ".txt"
                p, x = IO.read_out(in_file + ".out")
                IO.run(filename=in_file, p_val=p, m=metodo, tol=0, niter=k+1,
                       o=DIR_OUT, save_as=f"{test}_m{metodo}_i{k}")


def medir_errores():

    with open(RESULTADOS, 'a', encoding="utf-8") as file:

        for test in TESTS:

            in_file = TESTS_DIR + test + ".txt"
            p, xe   = IO.read_out(in_file + ".out")

            _, _, W = IO.read_in(in_file)
            A = utils.W_to_A(W, p)

            for metodo in METODOS:

                for k in range(0, ITER, STEP):

                    # error absoluto
                    res_file = DIR_OUT + f"{test}_m{metodo}_i{k}" + ".out"
                    _p, x = IO.read_out(res_file)

                    assert(abs(p - _p) < 1e-4)

                    error_abs = np.linalg.norm(x - xe, 1)
                    error_max = np.linalg.norm(x - xe, np.inf)

                    # error relativo
                    Ax = A @ x.T
                    error_rel = np.linalg.norm(Ax - x, 1)
                    
                    cols = FMT_COLS.format(test, metodo, k, 
                                           error_abs, error_rel, error_max)
                    file.write(cols)


#
# MAIN
#

if __name__ == "__main__":

    correr_pagerank()

    IO.createCSV(RESULTADOS, COLS)
    medir_errores()

    df = pd.read_csv(RESULTADOS)

    for test in df.test.unique():        
        
        save_as = GRAFICOS.format(test=test)
        test_data = df.query(f"test == '{test}'")
        x = test_data.n_iter
        y = test_data.error_abs
        hue = test_data.metodo
        log = np.any(y > 0)
        utils.graficar(x, y, hue, "iteraciones", "error absoluto", save_as, log=log)
