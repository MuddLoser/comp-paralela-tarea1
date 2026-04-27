#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm> // Para std::min

void multiplicacionBloques(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C, int n, int b) {
    // Inicializar C en cero (importante porque usamos +=)
    std::fill(C.begin(), C.end(), 0.0);

    // Los tres ciclos superiores recorren las submatrices (bloques)
    for (int i0 = 0; i0 < n; i0 += b) {
        for (int k0 = 0; k0 < n; k0 += b) {
            for (int j0 = 0; j0 < n; j0 += b) {

                // Los tres ciclos internos realizan la multiplicación del bloque b x b
                // Usamos std::min para manejar casos donde n no es divisible exactamente por b
                for (int i = i0; i < std::min(i0 + b, n); i++) {
                    for (int k = k0; k < std::min(k0 + b, n); k++) {
                        double r = A[i * n + k]; // Localidad espacial: A se lee una vez por bloque
                        for (int j = j0; j < std::min(j0 + b, n); j++) {
                            // Localidad espacial: C y B se acceden por filas (contiguo)
                            C[i * n + j] += r * B[k * n + j];
                        }
                    }
                }

            }
        }
    }
}

int main() {
    int n = 1024; 
    int b = 64; // Tamaño del bloque (puedes probar con 16, 32, 64, 128)
    
    std::vector<double> A(n * n, 1.0);
    std::vector<double> B(n * n, 2.0);
    std::vector<double> C(n * n, 0.0);

    auto start = std::chrono::high_resolution_clock::now();
    multiplicacionBloques(A, B, C, n, b);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;
    std::cout << "Tiempo Bloques (n=" << n << ", b=" << b << "): " << diff.count() << " s" << std::endl;

    return 0;
}