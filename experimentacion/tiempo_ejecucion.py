import base.IO as IO
import base.utils as utils

import os
import numpy as np
import pandas as pd


"""
descripcion: 
    evaluación de los tiempos de ejecución de PageRank para una serie de casos 
    de test. Se compararon tres implementaciones del algoritmo -sobre los 
    métodos iterativos de Gauss-Seidel y Jacobi, y sobre el método exacto de
    Eliminiación Gaussiana-.

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
COLS        = 'test,rep,metodo,tiempo'
FMT_COLS    = "{0},{1},{2},{3}\n"

SUMMARY     = DIR + "tiempo-ejecucion-sumario.csv"
GRAFICOS    = DIR + "tiempo-ejecucion.png"


#
# VAR
#

TOL_EG  = 1e-6
EPSILON = 0.01

REPS    = 10
TESTS   = [x[:-4] for x in os.listdir(TESTS_DIR) if \
           os.path.isfile(TESTS_DIR + x) and x[-4:] == ".txt"]
METODOS = ['GS', 'J','EG']
STEP    = 1
ITER    = 1e8


#
# UTILS
#

def correr_pagerank():

    for test in TESTS: 
        print("--------------------------------")
        print("==>", test)

        in_file = TESTS_DIR + test + ".txt"
        p, res_catedra = IO.read_out(in_file + ".out")

        print("- Calculando error EG:")
        IO.run(filename=in_file, p_val=p, m='EG', tol=TOL_EG, 
            o=DIR_OUT, save_as=f"{test}_EG")

        _, res_eg = IO.read_out(DIR_OUT + f"{test}_EG" + ".out")
        error_eg = np.linalg.norm(res_catedra - res_eg, 1)

        print("- Buscando tolerancia GS:")
        tol_gs = findCorrectTol(test, error_eg, 'GS')

        print("- Buscando tolerancia J:")
        tol_j  = findCorrectTol(test, error_eg, 'J')

        for rep in range(REPS):
            print("- Corriendo iteración:", str(rep+1) + '/' + str(REPS))
            IO.run(filename=in_file, p_val=p, m='EG', tol=TOL_EG, 
                o=DIR_OUT, save_as=f"{test}_EG_{rep}", time=True)
            IO.run(filename=in_file, p_val=p, m='GS', tol=tol_gs, niter=ITER,
                o=DIR_OUT, save_as=f"{test}_GS_{rep}", time=True)
            IO.run(filename=in_file, p_val=p, m='J',  tol=tol_j, niter=ITER,
                o=DIR_OUT, save_as=f"{test}_J_{rep}",  time=True)



def findCorrectTol(test, error_eg, metodo):

    in_file = TESTS_DIR + test + ".txt"
    p, res_catedra = IO.read_out(in_file + ".out")

    tol_i, tol_j = 0, 1
    while tol_i < tol_j :
        tol_k = (tol_i + tol_j) / 2

        IO.run(filename=in_file, p_val=p, m=metodo, tol=tol_k, 
            o=DIR_OUT, save_as=f"{test}_{metodo}")

        _, res_metodo = IO.read_out(DIR_OUT + f"{test}_{metodo}" + ".out")
        error = np.linalg.norm(res_catedra - res_metodo, 1)   
        
        if(abs(error - error_eg) < EPSILON * error_eg): 
            return tol_k
        elif(error > error_eg):
            tol_j = tol_k
        else:
            tol_i = tol_k
            
    return tol_j


def medir_tiempos():

    with open(RESULTADOS, 'a', encoding="utf-8") as file:

        for test in TESTS:
            print("--------------------------------")
            print("==>", test)            

            for metodo in METODOS:
                
                total = 0
                for rep in range(REPS):
                    time_file = DIR_OUT + f"{test}_{metodo}_{rep}" + ".time"
                    t = IO.read_time(time_file)
                    total += t.microseconds[1]
                    
                    cols = FMT_COLS.format(test, rep, metodo, t.microseconds[1])
                    file.write(cols)

                print("-", metodo + ":", str(float((total/REPS)/1e6)), "segundos") 


#
# MAIN
#

if __name__ == "__main__":

    correr_pagerank()

    IO.createCSV(RESULTADOS, COLS)
    medir_tiempos()

    df = pd.read_csv(RESULTADOS)

    df.describe().to_csv(SUMMARY)

    x = df.test.str.replace("test_", "").replace("aleatorio_desordenado", "desordenado")
    y = df.tiempo / 1e3
    hue = df.metodo.replace({
        "GS":"Gauss-Seidel", 
        "J": "Jacobi", 
        "EG": "Eliminación Gaussiana"})
    utils.graficar_barra(x, y, hue, "tests", "tiempo (ms)", GRAFICOS, log=True)
