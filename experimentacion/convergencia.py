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
COLS        = 'test,metodo,n_iter,trial,error_abs,error_rel,error_max'
FMT_COLS    = "{0},{1},{2},{3},{4},{5},{6}\n"

SUMARIO     = DIR + "sumario.csv"
SUMARIO_100 = DIR + "promedio-100.csv"

GRAFICOS    = DIR + "convergencia_{test}.png"


#
# VAR
#

TESTS   = [x[:-4] for x in os.listdir(TESTS_DIR) if \
           os.path.isfile(TESTS_DIR + x) and x[-4:] == ".txt"]
METODOS = ['GS', 'J']
ITER    = int(1e2)
STEP    = 1
TRIALS  = 10

#
# UTILS
#

def correr_pagerank():

    for test in TESTS: 

        in_file = TESTS_DIR + test + ".txt"
        p, x = IO.read_out(in_file + ".out")
            
        for metodo in METODOS:

            print(f"corriendo test: {test}, metodo: {metodo}.")    
            
            for k in range(0, ITER, STEP):

                for t in range(TRIALS):
                    
                    IO.run(filename=in_file, p_val=p, m=metodo, tol=0, niter=k+1,
                        o=DIR_OUT, save_as=f"{test}_m{metodo}_i{k+1}_t{t}")


def medir_errores():

    with open(RESULTADOS, 'a', encoding="utf-8") as file:

        for test in TESTS:

            in_file = TESTS_DIR + test + ".txt"
            p, xe   = IO.read_out(in_file + ".out")

            _, _, W = IO.read_in(in_file)
            A = utils.W_to_A(W, p)

            for metodo in METODOS:

                for k in range(0, ITER, STEP):

                    for t in range(TRIALS):

                        # error absoluto
                        res_file = DIR_OUT + f"{test}_m{metodo}_i{k+1}_t{t}" + ".out"
                        _p, x = IO.read_out(res_file)

                        assert(abs(p - _p) < 1e-4)

                        error_abs = np.linalg.norm(x - xe, 1)
                        error_max = np.linalg.norm(x - xe, np.inf)

                        # error relativo
                        Ax = A @ x.T
                        error_rel = np.linalg.norm(Ax - x, 1)
                        
                        cols = FMT_COLS.format(test, metodo, k, t,
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

    summary = df.groupby(["n_iter", "metodo", "test"])["error_abs"].describe()
    summary.to_csv(SUMARIO)

    summary100 = df.query(f"n_iter == 99")\
        .pivot_table(
            index="test", 
            columns="metodo", 
            values=["error_abs", "error_max"], 
            aggfunc="mean")
    summary100.to_csv(SUMARIO_100)

    for test in df.test.unique():        
        
        save_as = GRAFICOS.format(test=test)
        test_data = df.query(f"test == '{test}'")
        x = test_data.n_iter + 1
        y = test_data.error_abs
        hue = test_data.metodo.replace({"GS":"Gauss-Seidel", "J": "Jacobi"})
        log = np.any(y > 0)
        utils.graficar(x, y, hue, "iteraciones", "error absoluto", save_as, log=log)
