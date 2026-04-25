# Multiplicación de Matrices: Algoritmos Secuenciales y Paralelos

## Para ejecutar los algoritmos paralelos

```bash
# Compilar el código
gcc -O2 -fopenmp -std=c99 -o matmul2 matmul2.c -lm

# Ejecutar el programa. El formato es ./matmul [N], siendo N el tamaño de
# las matrices a multiplicar. N debe ser potencia de 2 y por defecto es 512
./matmul2
```
