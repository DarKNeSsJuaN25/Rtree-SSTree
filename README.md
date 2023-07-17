# Comparacion en tiempo de busqueda del SSTree y RTree

## Compilacion
### RTree 
Para el RTree, ejecutar en el archivo rtree.cpp: `g++ rtree.cpp -o a -lSDL2` y luego `./a`.
### SSTree
Para el SSTree, ejecutar en el archivo main.py.

## Resultados

Mostramos la grafica de resultados para el total de puntos de 100000.




## Explicacion
El RTree es una estructura de datos espaciales que está diseñada para manejar eficientemente consultas en espacios multidimensionales. Utiliza técnicas de particionamiento espacial para organizar los datos y reducir el número de accesos necesarios durante las búsquedas.

Sin embargo, a medida que aumenta el número de dimensiones, el espacio de búsqueda se vuelve más disperso y el efecto del llamado "curse of dimensionality" se vuelve más pronunciado. El "curse of dimensionality" se refiere al fenómeno de que, en espacios de alta dimensionalidad, la mayoría de los puntos están muy lejos unos de otros, lo que dificulta la identificación de los vecinos más cercanos.

Por otro lado, el SSTree es una estructura de datos basada en árboles de búsqueda que se utiliza comúnmente para buscar vecinos cercanos en espacios de alta dimensionalidad. A diferencia del RTree, el SSTree utiliza técnicas de ordenación espacial y búsqueda basada en distancia para mejorar la eficiencia de las consultas en espacios de alta dimensionalidad.

En general, debido a la maldicion de la dimensionalidad y las características de las estructuras de datos, es más probable que el RTree tenga un mayor tiempo de búsqueda que el SSTree a medida que aumenta el número de dimensiones en un escenario de búsqueda de los 5 vecinos más cercanos en un millón de datos aleatorios con M=50.
