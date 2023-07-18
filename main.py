import numpy as np
import time
import matplotlib.pyplot as plt
from SSTree2 import SSTree

def main():
    # Parámetros
    M = 50
    k = 5
    max_dimensions = 10
    num_data = 1000000

    # Generar datos aleatorios
    data = np.random.rand(num_data, max_dimensions)
    # Rellenar las dimensiones faltantes con ceros
    data = np.pad(data, ((0, 0), (0, max_dimensions - data.shape[1])), mode='constant')

    # Crear el árbol
    tree = SSTree(M=M, m=M//2)

    # Insertar los datos en el árbol
    for i in range(num_data):
        tree.insert(data[i])
        print(i)

    # Realizar búsquedas con diferentes dimensiones
    dimensions = range(1, max_dimensions+1)
    search_times = []

    for d in dimensions:
        query = np.random.rand(d)
        # Rellenar las dimensiones faltantes del query con ceros
        query = np.pad(query, (0, max_dimensions - d), mode='constant')

        start_time = time.time()
        tree.knn(query, k=k)
        end_time = time.time()

        search_time = (end_time - start_time) * 1000 
        search_times.append(search_time)

    # Graficar los resultados
    plt.plot(dimensions, search_times, marker='o')
    plt.xlabel('Número de dimensiones')
    plt.ylabel('Tiempo de búsqueda (segundos)')
    plt.title('Tiempo de búsqueda vs. Número de dimensiones')
    plt.show()

if __name__ == "__main__":
    main()
