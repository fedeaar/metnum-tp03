

import base.IO as IO
import base.utils as utils

import random
import sys
from itertools import product

import numpy as np
import pandas as pd


"""
descripcion: 
	Similar a 'tiempo_ejecucion.py'.
	Genera matrices de cierta familia y evalua tiempo de ejecución en proporción a su densidad.
"""



#
# IO
#

DIR = "./resultados/densidad/"
DIR_IN = DIR + "in/"
DIR_OUT = DIR + "out/"

EXPERIMENTO = "densidad"

COLS        = 'familia,rep,tam,dens,metodo,tiempo'
FMT_COLS    = "{0},{1},{2},{3},{4},{5}\n"


#
# VAR
#

TOL_EG  = 1e-6
EPSILON = 0.01

REPS    = 10
METODOS = ['GS', 'J','EG']
STEP    = 1
ITER    = 1e8

TAM = 1000
DENSIDADES = [0.01, 0.05, 0.1, 0.3, 0.5, 0.9, 1]

TESTS = [
	"red_sumidero",
	"uno_a_todos",
	"viborita",
	"todo_con_todo",
	"aleatorio"
]




#
# UTILS
#

def correr_pagerank(test, rep):

	print("- - - - - - - - - - - - - - - -")
	print("==>", test)

	in_file = DIR_IN + test + f"_{rep}.txt"
	p = 0.5

	print("- Calculando error EG...")
	IO.run(filename=in_file, p_val=p, m='EG', tol=TOL_EG, 
		o=DIR_OUT, save_as=f"{test}_EG")

	_, res_eg = IO.read_out(DIR_OUT + f"{test}_EG" + ".out")
	error_eg = error_relativo(in_file, p, res_eg)

	print("- Buscando tolerancia GS...")
	tol_gs = find_correct_tol(test, error_eg, 'GS', p)

	print("- Buscando tolerancia J...")
	tol_j  = find_correct_tol(test, error_eg, 'J', p)

	print("- Corriendo iteración:", str(rep+1) + '/' + str(REPS))
	IO.run(filename=in_file, p_val=p, m='EG', tol=TOL_EG, 
	o=DIR_OUT, save_as=f"{test}_EG_{rep}", time=True)
	IO.run(filename=in_file, p_val=p, m='GS', tol=tol_gs, niter=ITER,
	o=DIR_OUT, save_as=f"{test}_GS_{rep}", time=True)
	IO.run(filename=in_file, p_val=p, m='J',  tol=tol_j, niter=ITER,
	o=DIR_OUT, save_as=f"{test}_J_{rep}",  time=True)


def find_correct_tol(test, error_eg, metodo, p):

	in_file = DIR_IN + test + ".txt"

	tol_i, tol_j = 0, 1
	while tol_i < tol_j :
		tol_k = (tol_i + tol_j) / 2

		IO.run(filename=in_file, p_val=p, m=metodo, tol=tol_k, 
			o=DIR_OUT, save_as=f"{test}_{metodo}")

		_, res_metodo = IO.read_out(DIR_OUT + f"{test}_{metodo}" + ".out")
		error = error_relativo(in_file, p, res_metodo)

		if(abs(error - error_eg) < EPSILON * error_eg): 
			return tol_k
		elif(error > error_eg):
			tol_j = tol_k
		else:
			tol_i = tol_k
            
	return tol_j



def medir_tiempos(test, tam, dens):

	with open(DIR + f"densidad_{test}.csv", 'a', encoding="utf-8") as file:
		print("--------------------------------")

		for metodo in METODOS:

			total = 0
			for rep in range(REPS):
				time_file = DIR_OUT + f"{test}_{tam}_{dens}_{metodo}_{rep}" + ".time"
				t = IO.read_time(time_file)
				total += t.microseconds[1]
				
				cols = FMT_COLS.format(test, rep, tam, dens, metodo, t.microseconds[1])
				file.write(cols)

			print("-", metodo + ":", str(float((total/REPS)/1e6)), "segundos") 


def generar_test(tipo, tam, dens, rep):
	out = f"{tam}"

	if tipo == "red_sumidero":
		n = int(dens*tam)
		out += f"\n{n}"
		js = random.sample(range(1, tam+1), n)
		for i in range(0, n):
			out += f"\n{js[i]} {js[0]}"
	elif tipo == "uno_a_todos":
		n = int(dens*tam)
		out += f"\n{n}"
		js = random.sample(range(1, tam+1), n)
		for i in range(0, n):
			out += f"\n{js[0]} {js[i]}"
	elif tipo == "viborita":
		n = int(dens*tam)
		out += f"\n{n}"
		js = random.sample(range(1, tam+1), n)
		for i in range(1, n):
			out += f"\n{js[i-1]} {js[i]}"
	elif tipo == "todo_con_todo":
		#n = int(0.5 * (1 + np.sqrt(1 + 4*dens)))	# aprox. la cant. de nodos para que haya 'dens' aristas
		n = int(dens * tam * (tam - 1))					
		js = random.sample(range(1, tam+1), tam)
		ady = ""
		total = 0
		for i in range(0, tam):
			if total >= n: 
				break
			for j in range(0, i):
				ady += f"\n{js[i]} {js[j]}"
				ady += f"\n{js[j]} {js[i]}"
				total += 2
		out += f"\n{total}" + ady

	elif tipo == "aleatorio":
		n = int(dens * tam * (tam - 1))
		out += f"\n{n}"
		tups = random.sample(list(product(range(1, tam+1), range(1, tam+1))), tam * tam)
		total = 0
		i = 0
		while total < n and i < tam * tam:
			if tups[i][0] != tups[i][1]:
				out += f"\n{tups[i][0]} {tups[i][1]}"
				total += 1
			i += 1

	open(f"{DIR_IN}{tipo}_{tam}_{dens}_{rep}.txt", "w").write(out)


#
#	ERROR
#

def resolver_con_numpy(filename, p):

	_, _, W = IO.read_in(filename)
	A = utils.W_to_A(W, p)
	vals, _ = np.linalg.eig(A)

	vals = ''.join(["\n" + str(x) for x in vals])

	return vals


def error_relativo(filename, p, x):
	
	_, _, W = IO.read_in(filename)
	A = utils.W_to_A(W, p)
	Ax = A @ x
	error_rel = np.linalg.norm(Ax - x, 1)

	return error_rel
	



#
# MAIN
#

if __name__ == "__main__":

	# densidades = [ int(d * TAM) for d in DENSIDADES ]


	for test in TESTS:
		IO.createCSV(DIR + f"densidad_{test}.csv", COLS)

		for dens in DENSIDADES:
			
			for rep in range(REPS):
				generar_test(test, TAM, dens, rep)

				correr_pagerank(f"{test}_{TAM}_{dens}", rep)

			medir_tiempos(test, TAM, dens)


		df = pd.read_csv(DIR + f"densidad_{test}.csv")
		df.describe().to_csv(DIR + f"densidad_{test}_sumario.csv")

		x = df.dens / TAM
		y = df.tiempo / 1e3
		hue = df.metodo.replace({
			"GS":"Gauss-Seidel", 
			"J": "Jacobi", 
			"EG": "Eliminación Gaussiana"})
		utils.graficar(x, y, hue, "densidad", "tiempo (ms)", DIR + f"densidad_{test}.png", log=True)
