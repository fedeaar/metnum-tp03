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

EXPERIMENTO = "tiempo-de-ejecucion"
DIR_IN, DIR_OUT, DIR = IO.createInOut(EXPERIMENTO, delete=True)

RESULTADOS  = DIR + "tiempo-de-ejecucion.csv"
COLS        = 'test,metodo,tiempo'
FMT_COLS    = "{0},{1},{2}\n"

GRAFICOS    = DIR + "t-ejecucion_{test}.png"

TOL_EG  = 1e-6
EPSILON = 1e-7
#
# VAR
#

REPS    = 100  
TESTS   = [x[:-4] for x in os.listdir(TESTS_DIR) if \
           os.path.isfile(TESTS_DIR + x) and x[-4:] == ".txt"]
METODOS = ['GS', 'J','EG']
STEP    = 1
ITER    = 1

#
# UTILS
#

def correr_pagerank():

    for test in TESTS: 
        in_file = TESTS_DIR + test + ".txt"
        p, res_catedra = IO.read_out(in_file + ".out")

        IO.run(filename=in_file, p_val=p, m='EG', tol=EPSILON, o=DIR_OUT, save_as=f"{test}_EG")

        _, res_eg = IO.read_out(DIR_OUT + f"{test}_EG" + ".out")
        error_eg = np.linalg.norm(res_catedra - res_eg, 1)

        tol_gs = findCorrectTol(test, error_eg, 'GS')
        tol_j  = findCorrectTol(test, error_eg, 'J')
        for rep in REPS:
            IO.run(filename=in_file, p_val=p, m='EG', tol=EPSILON, o=DIR_OUT, save_as=f"{test}_EG_{rep}", time=True)
            IO.run(filename=in_file, p_val=p, m='GS', tol=tol_gs,  o=DIR_OUT, save_as=f"{test}_GS_{rep}", time=True)
            IO.run(filename=in_file, p_val=p, m='J',  tol=tol_j,   o=DIR_OUT, save_as=f"{test}_J_{rep}",  time=True)


  

def findCorrectTol(test, error_eg, metodo):
    in_file = TESTS_DIR + test + ".txt"
    p, res_catedra = IO.read_out(in_file + ".out")

    tol_i, tol_j = 1e-2, 1e-8
    while tol_i + EPSILON < tol_j :
        tol_k = (tol_i + tol_j) / 2
        IO.run(filename=in_file, p_val=p, m=metodo, tol=tol_k, o=DIR_OUT, save_as=f"{test}_{metodo}")
        _, res_metodo = IO.read_out(DIR_OUT + f"{test}_{metodo}" + ".out")
        error_gs = np.linalg.norm(res_catedra - res_metodo, 1)   
        if(abs(error_gs - error_eg) < EPSILON):
            return tol_k
        elif(error_gs > error_eg):
            tol_i = tol_k
        else:
            tol_j = tol_k
            
    return tol_i


def medir_tiempos():

    with open(RESULTADOS, 'a', encoding="utf-8") as file:

        for test in TESTS:

            in_file = TESTS_DIR + test + ".txt"
            p, res_catedra = IO.read_out(in_file + ".out")


            for metodo in METODOS:

                for k in range(0, ITER, STEP):

                    # error absoluto
                    res_file = DIR_OUT + f"{test}_m{metodo}_i{k}" + ".out"
                    _p, x = IO.read_out(res_file)

                    assert(abs(p - _p) < 1e-4)

                    error_abs = np.linalg.norm(x - xe, 1)
                    error_max = np.linalg.norm(x - xe, np.inf)

                    
                    cols = FMT_COLS.format(test, metodo, k, 
                                           error_abs, error_rel, error_max)
                    file.write(cols)











if __name__ == "__main__":

    correr_pagerank()

    # IO.createCSV(RESULTADOS, COLS)
    # medir_errores()

    # df = pd.read_csv(RESULTADOS)

    # for test in df.test.unique():        
        
    #     save_as = GRAFICOS.format(test=test)
    #     test_data = df.query(f"test == '{test}'")
    #     x = test_data.n_iter
    #     y = test_data.error_abs
    #     hue = test_data.metodo.replace({"GS":"Gauss-Seidel", "J": "Jacobi"})
    #     log = np.any(y > 0)
    #     utils.graficar(x, y, hue, "iteraciones", "error absoluto", save_as, log=log)
