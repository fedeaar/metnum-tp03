import numpy as np
import pandas as pd

import pygraphviz
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(rc={'figure.figsize':(14, 7)}, font="Times New Roman")


# GLOBALES
EPSILON = 1e-4


def random_conectividad(n, cantidad, dirigida=True, seed=None):

    sz = (n - 1) * n if dirigida else n * (n - 1) // 2
    W = np.zeros(sz, dtype=np.double)
    W[:cantidad] = 1
    
    rng = np.random.default_rng(seed if seed else None)
    rng.shuffle(W)
    
    if dirigida:
        for i in range(0, sz, n + 1):
            W = np.concatenate((W[:i], [0], W[i:]))
        W = np.append(W, 0)
    else:
        for i in range(n):
            W = np.concatenate((W[:i*n], [0]*(i+1), W[i*n:]))
    
    W = W.reshape((n, n))
    
    if not dirigida:
        W = W + W.T

    return W 


def random_matriz(n, m=None, r=(1, 100)):
    return np.random.randint(r[0], r[1], (n, m if m else n))


def norma(x, norma):
    return np.linalg.norm(x, norma)


def metodo_potencia(A, niter=10000, epsilon=1e-6, x={}):
    n = A.shape[0]
    
    if(len(x) == 0): x = np.random.rand(n, 1)

    z = np.zeros((n, 1))
    if np.allclose(x, z, epsilon):
        x = z
        x[0] = 1
    else:
        x = x / norma(x, 2)

    for _ in range(niter):
        y = A @ x
        n = norma(y, 2)
        if norma(x - y / n, 2) < epsilon:
            break
        x = y / n
    a = np.dot(x.T, (A @ x)) / np.dot(x.T, x)

    return a[0,0], x


def metodo_deflacion(A, k, niter=10000, epsilon=1e-6):

    n = A.shape[0]
    A = A.copy()
    eigs = []
    vecs = np.zeros((n, k))

    for i in range(k):
        a, v = metodo_potencia(A, niter, epsilon)
        eigs.append(a)
        vecs[:, i] = v.T
        A = A - a * (v @ v.T)
        
    return np.array(eigs), vecs


def eig(A):

    w, V = np.linalg.eig(A)
    idx = np.argsort(np.abs(w))[::-1]
    w = w[idx]
    V = V[:, idx]
    return w, V


def corr(x, y, epsilon=1e-16):
    
    xc = x - x.mean()
    yc = y - y.mean()
    a = xc.T @ yc
    b = np.linalg.norm(xc, 2) * np.linalg.norm(yc, 2)
    return a / b if abs(b) >= epsilon else 0


def graficar(x, y, hue, xaxis, yaxis, filename, units=None):
    
    plt.figure()
    df   = pd.DataFrame({"x":x, "y":y, "hue":hue})
    kwargs = {
        "data":df, 
        "x":"x", 
        "y":"y", 
        "hue":"hue", 
        "units":units
    }
    if units:
        kwargs["estimator"] = None
    plot = sns.lineplot(**kwargs)
    

    plot.set_xlabel(xaxis, fontsize=18, labelpad=12)
    plot.set_ylabel(yaxis, fontsize= 18, labelpad=20) 
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.legend(title=None)

    fig = plot.get_figure()
    fig.savefig(filename)
    plt.close(fig)

def graficar2(x, y, hue, x2, y2, hue2, xaxis, yaxis, filename):
    
    plt.figure()
    df   = pd.DataFrame({"x":x, "y":y, "hue":hue})
    plot = sns.lineplot(data=df, x="x", y="y", hue="hue")

    df   = pd.DataFrame({"x":x, "y":y, "hue":hue})
    plot = sns.lineplot(data=df, x="x", y="y", hue="hue")
    

    plot.set_xlabel(xaxis, fontsize=18, labelpad=12)
    plot.set_ylabel(yaxis, fontsize= 18, labelpad=20) 
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.legend(title=None)

    fig = plot.get_figure()
    fig.savefig(filename)
    plt.close(fig)


def graficar_grafo(A, filename, 
    size=(10, 10),
    node_size=300,
    font_size=12, 
    node_color="lightgray", 
    edge_color='black',
    font_color="k", 
    with_labels=True):

    G = nx.from_numpy_array(A)
    f = plt.figure(figsize=size)
    layout = nx.nx_agraph.graphviz_layout(G, 'neato')
    options={
         'node_color': node_color,
         'node_size': node_size,
         'font_size': font_size,
         'font_color': font_color,
         'edge_color': edge_color,
         'with_labels': with_labels
    }
    nx.draw(G, layout, ax=f.add_subplot(), **options)

    f.savefig(filename)
    plt.close(f)


def armarMatriz(inicial, n):
    inicial.extend(np.random.randint(-n+1, n, size=n-len(inicial)))  
    D = np.diag(inicial)
    
    u = nml(np.random.rand(n, 1))
    H = np.eye(n) - 2 * (u @ u.T)
    S = H @ D @ H.T

    a, V = eig(S)
    a = a.astype(float)
    V = V.astype(float)
    if(a[0] < a[1]): V.T[[0,1]] = V.T[[1, 0]]

    return S, V, a


def armarRandom(n):
    x = np.random.randint(-100, 100, size=n)
    x = x.astype(float)
    x = nml(x)
    return x 


def n2(v):
    return np.linalg.norm(v, 2)


def nml(x):
    return x / n2(x)


def rayleigh(A, v):
    return (np.dot(v.T, (A @ v)) / np.dot(v.T, v))[0,0]


def alt_potencia(A, niter=10000, epsilon=1e-6, x={}):
    niter //= 2
    if(len(x) == 0): x = np.random.rand(A.shape[0], 1)
    
    x = nml(x)
    for _ in range(niter):
        y = nml(A @ nml(A @ x))
        if n2(x - y) < epsilon : break
        x = y

    a = rayleigh(A, x)

    x2 = nml(A @ x) - x
    x3 = nml(A @ x) + x
    if n2(x2) > EPSILON and n2(x3) > EPSILON and niter > 100:
        print("entree...")
        x2 = nml(x2)
        a2 = rayleigh(A, x2)
        return a2, x2

    return a, x


def alt_deflacion(A, k, niter=10000, epsilon=1e-6):

    n = A.shape[0]
    A = A.copy()
    eigs = []
    vecs = np.zeros((n, k))

    for i in range(k):
        if(n2(A) < EPSILON): break
        a, v = alt_potencia(A, niter, epsilon)
        eigs.append(a)
        vecs[:, i] = v.T
        A = A - a * (v @ v.T)
        
    return np.array(eigs), vecs
    