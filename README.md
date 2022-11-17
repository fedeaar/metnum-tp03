## TP3: Métodos iterativos

grupo 18 - Arienti, Vekselman, Lakowsky, Kovo


<br>

### Estructura del repo

El repositorio cuenta con los siguientes archivos y carpetas:

- 'catedra' -  Los archivos fuente ys casos de test provistos por la cátedra.

- 'implementacion' - El código fuente para la solución propuesta.

- 'experimentos' - El material correspondiente a todos los experimentos mencionados en el informe. Incluye scripts y archivos resultado. Se omitieron los archivos intermedios, los mismo se pueden regenerar a partir de los scripts.

- el informe.


<br>

### Nota de entorno

La compilación de los ejecutables se realizó por medio de `CMAKE` y `MAKE` en `WSL (windows)`. Recomendamos utilizar las mismas aplicaciones y trabajar en un entorno de linux. Más abajo explicamos cómo. 

Así también, los scripts se pensaron para ser ejecutados por medio de WSL. De compilar para windows ó de ejecutarse directo en linux se deberán modificar las variables globales al comienzo del archivo `./experimentos/utils/IO.py`. En el mismo se detalla cómo. 



<br>

### Cómo crear los archivos ejecutables

Para este procedimiento se asume que trabajaremos en bash. Desde la raiz del repo procederemos de la siguiente forma:

1. creamos la carpeta para los ejecutables
    > $ mkdir build
    
2. nos movemos adentro
    > $ cd build

3. creamos el cmake
    > $ cmake ../implementacion

4. creamos el ejecutable principal (el mismo se requerirá para los experimentos)
    > $ make ./tp3 

5. ejecutamos

    > $ ./tp3 ../catedra/tests-pagerank/test_aleatorio.txt 0.76

Notamos que la ejecución de los experimentos requiere que el ejecutable `./tp3` se encuentre en `./build`.


<br>

### IO ./tp3

El ejecutable permite trabajar con los siguientes parámetros.


Obligatorios (deben estar en orden):

- `*` (string): fuente del archivo de entrada. El mismo debe contener una lista de adjacencia (delimitada por ' ') con el siguiente encabezado: la cantidad de vértices (línea uno) y la cantidad de aristas (línea dos). Ejemplo de uso: `../catedra/tests-pagerank/test_aleatorio.txt`.

- `*` [0, 1): valor p.


Opcionales:

- `-m`: método a utilizar. `EG` (eliminación gaussiana), `GS` (gauss-seidel) ó  `J` (jacobi). Por default `EG`. Ejemplo de uso: `-m GS`. 

- `-tol`: Se define acorde a `-m`. Si `GS` ó `J`, entonces refiere a la tolerancia mínima de cambio entre pasos consecutivos a partir de los que una solución se considera válida. Si `EG`, refiere al mínimo valor a considerar no nulo dentro de la matriz a triangular. Por default 1e-4. Ejemplo de uso: `-tol 1e-6`.

- `-iter`: máxima cantidad de iteraciones a realizar. Por default 1e5. Aplica sólo para `GS` y `J`. Ejemplo de uso: `-iter 1e6`.

- `-o`: carpeta en la que se guardarán los archivos de salida. Por defecto, la misma donde se encuentra el ejecutable. Ejemplo de uso: `-o ../experimentos/`.

- `-as`: nombre con el que se guardarán los archivos de salida. Por defecto, este nombre será el mismo que el del archivo de entrada pero con una extensión distinta acorde al archivo a guardar. Para el archivo solución, esta extensión sera `.out`. Ejemplo de uso: `-as resultado`. 

- `-p`: la precisión con la que se guardaran los resultados, en el sentido de la cantidad de dígitos decimales después de la coma. Por defecto 15. Ejemplo de uso: `-p 8`.

- `-time`: Flag. Si se pasa éste parámetro, se guardará el tiempo de ejecución de distintas etapas del algoritmo (descontando las operaciones de IO) en un archivo CSV con la extensión `.time`.

- `-v`: Flag. Verbose. Si se pasa éste parámetro, se imprimirá en la consola información relevante durante la ejecución del programa.
