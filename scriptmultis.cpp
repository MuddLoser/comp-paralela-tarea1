#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <iomanip>

using namespace std;

// --- FUNCIONES DE APOYO (Suma/Resta/Verificación) ---
void sumar(const vector<double>& A, const vector<double>& B, vector<double>& C, int n) {
    for (int i = 0; i < n * n; i++) C[i] = A[i] + B[i];
}
void restar(const vector<double>& A, const vector<double>& B, vector<double>& C, int n) {
    for (int i = 0; i < n * n; i++) C[i] = A[i] - B[i];
}
bool verificar(const vector<double>& Ref, const vector<double>& Test, int n) {
    for (int i = 0; i < n * n; i++) {
        if (abs(Ref[i] - Test[i]) > 1e-6) return false;
    }
    return true;
}

// --- IMPLEMENTACIONES ---

void multClasica(const vector<double>& A, const vector<double>& B, vector<double>& C, int n) {
    fill(C.begin(), C.end(), 0.0);
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            double r = A[i * n + k];
            for (int j = 0; j < n; j++) {
                C[i * n + j] += r * B[k * n + j];
            }
        }
    }
}

void multBloques(const vector<double>& A, const vector<double>& B, vector<double>& C, int n, int b) {
    fill(C.begin(), C.end(), 0.0);
    for (int i0 = 0; i0 < n; i0 += b) {
        for (int k0 = 0; k0 < n; k0 += b) {
            for (int j0 = 0; j0 < n; j0 += b) {
                for (int i = i0; i < min(i0 + b, n); i++) {
                    for (int k = k0; k < min(k0 + b, n); k++) {
                        double r = A[i * n + k];
                        for (int j = j0; j < min(j0 + b, n); j++) {
                            C[i * n + j] += r * B[k * n + j];
                        }
                    }
                }
            }
        }
    }
}

void strassen(const vector<double>& A, const vector<double>& B, vector<double>& C, int n, int n0) {
    if (n <= n0) { multClasica(A, B, C, n); return; }
    int k = n / 2; int sz = k * k;
    vector<double> a11(sz), a12(sz), a21(sz), a22(sz), b11(sz), b12(sz), b21(sz), b22(sz);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            a11[i*k+j] = A[i*n+j]; a12[i*k+j] = A[i*n+j+k]; a21[i*k+j] = A[(i+k)*n+j]; a22[i*k+j] = A[(i+k)*n+j+k];
            b11[i*k+j] = B[i*n+j]; b12[i*k+j] = B[i*n+j+k]; b21[i*k+j] = B[(i+k)*n+j]; b22[i*k+j] = B[(i+k)*n+j+k];
        }
    }
    vector<double> m1(sz), m2(sz), m3(sz), m4(sz), m5(sz), m6(sz), m7(sz), tA(sz), tB(sz);
    sumar(a11, a22, tA, k); sumar(b11, b22, tB, k); strassen(tA, tB, m1, k, n0);
    sumar(a21, a22, tA, k); strassen(tA, b11, m2, k, n0);
    restar(b12, b22, tB, k); strassen(a11, tB, m3, k, n0);
    restar(b21, b11, tB, k); strassen(a22, tB, m4, k, n0);
    sumar(a11, a12, tA, k); strassen(tA, b22, m5, k, n0);
    restar(a21, a11, tA, k); sumar(b11, b12, tB, k); strassen(tA, tB, m6, k, n0);
    restar(a12, a22, tA, k); sumar(b21, b22, tB, k); strassen(tA, tB, m7, k, n0);
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            int idx = i * k + j;
            C[i*n+j] = m1[idx] + m4[idx] - m5[idx] + m7[idx];
            C[i*n+j+k] = m3[idx] + m5[idx];
            C[(i+k)*n+j] = m2[idx] + m4[idx];
            C[(i+k)*n+j+k] = m1[idx] - m2[idx] + m3[idx] + m6[idx];
        }
    }
}

// --- ORQUESTADOR DE EXPERIMENTOS ---

int main() {
    vector<int> tamanos = {256, 512, 1024, 2048}; // Puedes agregar 4096 si el server aguanta la RAM
    int repeticiones = 3;
    int b = 64;  // Bloque para tiling
    int n0 = 64; // Caso base Strassen

    cout << fixed << setprecision(5);
    cout << "N\tAlgoritmo\tPromedio(s)\tCorrecto?" << endl;
    cout << "--------------------------------------------------------" << endl;

    for (int n : tamanos) {
        vector<double> A(n * n, 1.2), B(n * n, 0.8), C_ref(n * n), C_test(n * n);

        // --- TEST CLÁSICA ---
        double totalTime = 0;
        for(int r=0; r<repeticiones; r++) {
            auto start = chrono::high_resolution_clock::now();
            multClasica(A, B, C_ref, n);
            auto end = chrono::high_resolution_clock::now();
            totalTime += chrono::duration<double>(end - start).count();
        }
        cout << n << "\tClasica\t\t" << totalTime/repeticiones << "\tSI" << endl;

        // --- TEST BLOQUES ---
        totalTime = 0;
        for(int r=0; r<repeticiones; r++) {
            auto start = chrono::high_resolution_clock::now();
            multBloques(A, B, C_test, n, b);
            auto end = chrono::high_resolution_clock::now();
            totalTime += chrono::duration<double>(end - start).count();
        }
        cout << n << "\tBloques\t\t" << totalTime/repeticiones << "\t" << (verificar(C_ref, C_test, n) ? "SI" : "NO") << endl;

        // --- TEST STRASSEN ---
        totalTime = 0;
        for(int r=0; r<repeticiones; r++) {
            auto start = chrono::high_resolution_clock::now();
            strassen(A, B, C_test, n, n0);
            auto end = chrono::high_resolution_clock::now();
            totalTime += chrono::duration<double>(end - start).count();
        }
        cout << n << "\tStrassen\t" << totalTime/repeticiones << "\t" << (verificar(C_ref, C_test, n) ? "SI" : "NO") << endl;
        cout << "--------------------------------------------------------" << endl;
    }

    return 0;
}