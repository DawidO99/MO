#include "lu_solver.hpp"
#include <cmath>
#include <stdexcept>

void decompose_LU(std::vector<std::vector<long double>>& A, std::vector<int>& pivot) {
    int N = A.size();
    pivot.resize(N);
    for (int i = 0; i < N; ++i) pivot[i] = i;

    for (int k = 0; k < N; ++k) {
        // Szukanie elementu g³ównego (czêœciowy wybór)
        long double maxVal = std::abs(A[k][k]);
        int maxRow = k;
        for (int i = k + 1; i < N; ++i) {
            if (std::abs(A[i][k]) > maxVal) {
                maxVal = std::abs(A[i][k]);
                maxRow = i;
            }
        }
        if (maxVal == 0.0L) throw std::runtime_error("Macierz osobliwa");

        // Zamiana wierszy
        if (maxRow != k) {
            std::swap(A[k], A[maxRow]);
            std::swap(pivot[k], pivot[maxRow]);
        }

        // Eliminacja
        for (int i = k + 1; i < N; ++i) {
            A[i][k] /= A[k][k];
            for (int j = k + 1; j < N; ++j) {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }
}

void solve_LU(const std::vector<std::vector<long double>>& LU,
    const std::vector<int>& pivot,
    const std::vector<long double>& b,
    std::vector<long double>& x)
{
    int N = LU.size();
    std::vector<long double> y(N);

    // Forward substitution: L * y = Pb
    for (int i = 0; i < N; ++i) {
        y[i] = b[pivot[i]];
        for (int j = 0; j < i; ++j)
            y[i] -= LU[i][j] * y[j];
    }

    // Backward substitution: U * x = y
    x.resize(N);
    for (int i = N - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < N; ++j)
            x[i] -= LU[i][j] * x[j];
        x[i] /= LU[i][i];
    }
}
