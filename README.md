# Multiplicación de Matrices: Algoritmos Secuenciales y Paralelos
## Para ejecutar los algoritmos secuenciales
# Compilar el código
# Usamos -O3 para habilitar optimizaciones de rendimiento y -std=c++17 por compatibilidad
g++ -O3 -std=c++17 -o program scriptmultis.cpp

# Ejecutar el programa
.\program

## Para ejecutar los secuenciales en el servidor Colcura (así se hizo)
# 1. Compilar el código en el entorno de alto rendimiento
g++ -O3 -std=c++17 -o bench_colcura scriptmultis.cpp

# 2. Ejecutar el experimento en segundo plano (Background)
# nohup asegura que el proceso no se detenga si se cierra la sesión
# > guarda la salida en un archivo y & lo envía a segundo plano
nohup ./bench_colcura > estadisticas_grupo3.txt &

# 3. Ver el avance del archivo de estadísticas
tail -f estadisticas_grupo3.txt

## Para ejecutar los algoritmos paralelos

```bash
# Compilar el código
gcc -O2 -fopenmp -std=c99 -o matmul2 matmul2.c -lm

# Ejecutar el programa. El formato es ./matmul [N], siendo N el tamaño de
# las matrices a multiplicar. N debe ser potencia de 2 y por defecto es 512
./matmul2
```
