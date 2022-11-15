import base.IO as IO
import base.utils as utils

import numpy as np
import pandas as pd


"""
descripcion: 
    evaluación de la convergencia de las soluciones de PageRank (ver tp1) para una 
    serie de casos de test. Se compararon dos implementaciones del algoritmo -sobre 
    los métodos iterativos de Gauss-Seidel y Jacobi- y se midió el error absoluto de 
    los resultados en función de la cantidad de iteraciones realizadas.
"""

#
# IN
#

TESTS = "../catedra/tests-pagerank"


# 
# OUT
#

EXPERIMENTO = "convergencia-iterativos"
DIR_IN, DIR_OUT, DIR = IO.createInOut(EXPERIMENTO, delete=True)

ERRORES = DIR + "error-convergencia.csv"
