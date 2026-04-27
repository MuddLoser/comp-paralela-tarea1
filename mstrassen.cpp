#include <iostream>
#include <vector>
#include <chrono>

// Funciones auxiliares para sumar y restar submatrices
void sumar(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C, int n) {
    for (int i = 0; i < n * n; i++) C[i] = A[i] + B[i];
}

void restar(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C, int n) {
    for (int i = 0; i < n * n; i++) C[i] = A[i] - B[i];
}

// Multiplicación clásica para el caso base (para no saturar la recursión)
void multiplicacionClasicaAux(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C, int n) {
    std::fill(C.begin(), C.end(), 0.0);
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            double r = A[i * n + k];
            for (int j = 0; j < n; j++) {
                C[i * n + j] += r * B[k * n + j];
            }
        }
    }
}

void strassen(const std::vector<double>& A, const std::vector<double>& B, std::vector<double>& C, int n, int n0) {
    // CASO BASE: Si la matriz es pequeña, usamos la clásica
    if (n <= n0) {
        multiplicacionClasicaAux(A, B, C, n);
        return;
    }

    int k = n / 2;
    int size = k * k;

    // Submatrices
    std::vector<double> a11(size), a12(size), a21(size), a22(size);
    std::vector<double> b11(size), b12(size), b21(size), b22(size);

    // División de A y B en 4 submatrices
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            a11[i * k + j] = A[i * n + j];
            a12[i * k + j] = A[i * n + j + k];
            a21[i * k + j] = A[(i + k) * n + j];
            a22[i * k + j] = A[(i + k) * n + j + k];

            b11[i * k + j] = B[i * n + j];
            b12[i * k + j] = B[i * n + j + k];
            b21[i * k + j] = B[(i + k) * n + j];
            b22[i * k + j] = B[(i + k) * n + j + k];
        }
    }

    // Los 7 productos de Strassen (M1 a M7)
    std::vector<double> m1(size), m2(size), m3(size), m4(size), m5(size), m6(size), m7(size);
    std::vector<double> tA(size), tB(size);

    // M1 = (A11 + A22) * (B11 + B22)
    sumar(a11, a22, tA, k); sumar(b11, b22, tB, k);
    strassen(tA, tB, m1, k, n0);

    // M2 = (A21 + A22) * B11
    sumar(a21, a22, tA, k);
    strassen(tA, b11, m2, k, n0);

    // M3 = A11 * (B12 - B22)
    restar(b12, b22, tB, k);
    strassen(a11, tB, m3, k, n0);

    // M4 = A22 * (B21 - B11)
    restar(b21, b11, tB, k);
    strassen(a22, tB, m4, k, n0);

    // M5 = (A11 + A12) * B22
    sumar(a11, a12, tA, k);
    strassen(tA, b22, m5, k, n0);

    // M6 = (A21 - A11) * (B11 + B12)
    restar(a21, a11, tA, k); sumar(b11, b12, tB, k);
    strassen(tA, tB, m6, k, n0);

    // M7 = (A12 - A22) * (B21 + B22)
    restar(a12, a22, tA, k); sumar(b21, b22, tB, k);
    strassen(tA, tB, m7, k, n0);

    // Combinar resultados en C
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            int idx = i * k + j;
            C[i * n + j] = m1[idx] + m4[idx] - m5[idx] + m7[idx]; // C11
            C[i * n + j + k] = m3[idx] + m5[idx];                 // C12
            C[(i + k) * n + j] = m2[idx] + m4[idx];               // C21
            C[(i + k) * n + j + k] = m1[idx] - m2[idx] + m3[idx] + m6[idx]; // C22
        }
    }
}

int main() {
    int n = 1024; // Strassen requiere potencias de 2
    int n0 = 64;  // Umbral del caso base (punto donde deja de ser recursivo)

    std::vector<double> A(n * n, 1.0);
    std::vector<double> B(n * n, 2.0);
    std::vector<double> C(n * n, 0.0);

    auto start = std::chrono::high_resolution_clock::now();
    strassen(A, B, C, n, n0);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;
    std::cout << "Tiempo Strassen (n=" << n << ", n0=" << n0 << "): " << diff.count() << " s" << std::endl;

    return 0;
}