import os
import shutil
import subprocess as sub
import numpy as np


# GLOBAL
EXE_PATH = '../build/tp2'      # si se compilo de otra manera, o con otro nombre, cambiar por la direccion correcta
WSL = True                     # dejar true solo si se utilizó wsl para compilar el programa, false sino


#
# IO files
#

def read_adylist(filename):
    
    with open(filename) as file:
        data = file.read().splitlines()
        links = [int(y) for x in data for y in x.split(' ')]
        n = np.max(links) if links else 1
        matriz = np.zeros((n, n))
        for i in range(0, len(links), 2):
            matriz[links[i+1] - 1][links[i] - 1] = 1

    return matriz


def write_adylist(filename, matrix):
    
    links = np.where(matrix != 0)
    links = tuple(zip(*links))
    text  = []
    for coord in links:
        text.append(str(coord[1] + 1) + " " + str(coord[0] + 1))

    np.savetxt(filename, text, delimiter="\n", fmt="%s")


def read_time(filename):

    with open(filename) as file:
        data = file.read().splitlines()
        scale = data[0]
        time = int(data[1])
        
    return scale, time


#
# IO experimentos
#

def createInOut(filename, delete=False):
    
    path = "./resultados/" + filename + "/"
    pathIn =  path + "in/"
    pathOut = path + "out/"
    if delete and os.path.exists(path): 
        shutil.rmtree(path)
    if not os.path.exists(pathIn):
        os.makedirs(pathIn)
    if not os.path.exists(pathOut):
        os.makedirs(pathOut)
        
    return pathIn, pathOut, path


def createCSV(filename, columnas): 

    with open(filename, "w", encoding="utf-8") as file:
        file.write(columnas + "\n")


#
# IO ./tp2
#

def run(filename, niter, tol, 
        f="matriz", m="base", o="./", p=15, save_as=None, time=False, verbose=False,
        exe_path=EXE_PATH):

    call_params = [
        "wsl" if WSL else "",   
        exe_path, 
        filename, str(niter), str(tol),
        "-f", f,
        "-m", m,
        "-o", o,
        "-p", str(p)
    ]
    if save_as: call_params.extend(["-as", save_as])
    if verbose: call_params.extend(["-v"])
    if time: call_params.extend(["-time"])

    sub.check_call(call_params)
