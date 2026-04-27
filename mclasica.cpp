#include <iostream>
#include <vector>
#include <chrono>

void multiplicacionClasica(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double suma = 0.0;
            for (int k = 0; k < n; k++) {
                // Acceso a matriz plana: fila * n + columna
                suma += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = suma;
        }
    }
}

int main() {
    int n = 1024; // Ejemplo de tamaño
    std::vector<double> A(n * n, 1.0);
    std::vector<double> B(n * n, 2.0);
    std::vector<double> C(n * n, 0.0);

    auto start = std::chrono::high_resolution_clock::now();
    multiplicacionClasica(A, B, C, n);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;
    std::cout << "Tiempo para n=" << n << ": " << diff.count() << " s" << std::endl;

    return 0;
}