#include "solver_crank.hpp"
#include "solver_utils.hpp"
#include "thomas_algorithm.hpp"
#include "lu_solver.hpp"
#include <vector>
#include <cmath>

void solve_crank_nicolson(int N, int M, long double t_max,
    std::vector<long double>& x_out,
    std::vector<long double>& u_out)
{
    long double h = a / N;
    long double dt = t_max / M;
    long double lambda = D * dt / (h * h);

    std::vector<long double> x(N + 1), u_curr(N + 1), u_next(N + 1);
    for (int i = 0; i <= N; ++i) {
        x[i] = r + i * h;
        u_curr[i] = 1.0L;
    }

    u_curr[0] = 0.0L;
    u_curr[N] = right_boundary_value(0.0L);

    std::vector<long double> aT(N - 1, -lambda / 2.0L);
    std::vector<long double> bT(N - 1, 1.0L + lambda);
    std::vector<long double> cT(N - 1, -lambda / 2.0L);
    std::vector<long double> d(N - 1);

    for (int t = 1; t <= M; ++t) {
        long double current_time = t * dt;

        // Zbuduj praw¹ stronê uk³adu (macierz jawna * u_curr)
        for (int i = 1; i < N; ++i) {
            d[i - 1] = (lambda / 2.0L) * u_curr[i - 1]
                + (1.0L - lambda) * u_curr[i]
                + (lambda / 2.0L) * u_curr[i + 1];
        }

        d[0] += (lambda / 2.0L) * u_curr[0]; // lewy brzeg
        d[N - 2] += (lambda / 2.0L) * right_boundary_value(current_time); // prawy brzeg

        std::vector<long double> result(N - 1);
        thomas_algorithm(aT, bT, cT, d, result);

        u_next[0] = 0.0L;
        u_next[N] = right_boundary_value(current_time);
        for (int i = 1; i < N; ++i) {
            u_next[i] = result[i - 1];
        }

        u_curr = u_next;
    }

    x_out = x;
    u_out = u_curr;
}

void solve_crank_LU(int N, int M, long double t_max,
    std::vector<long double>& x_out,
    std::vector<long double>& u_out)
{
    long double h = a / N;
    long double dt = t_max / M;
    long double lambda = D * dt / (h * h);

    std::vector<long double> x(N + 1), u_curr(N + 1), u_next(N + 1);
    for (int i = 0; i <= N; ++i) {
        x[i] = r + i * h;
        u_curr[i] = 1.0L;
    }

    u_curr[0] = 0.0L;
    u_curr[N] = right_boundary_value(0.0L);

    // Budowa pe³nej macierzy A (dla lewego cz³onu)
    std::vector<std::vector<long double>> A(N - 1, std::vector<long double>(N - 1, 0.0L));
    for (int i = 0; i < N - 1; ++i) {
        A[i][i] = 1.0L + lambda;
        if (i > 0) A[i][i - 1] = -lambda / 2.0L;
        if (i < N - 2) A[i][i + 1] = -lambda / 2.0L;
    }

    // LU dekompozycja macierzy A
    std::vector<int> pivot;
    decompose_LU(A, pivot);

    std::vector<long double> d(N - 1), result(N - 1);

    for (int t = 1; t <= M; ++t) {
        long double current_time = t * dt;

        // Budowa prawej strony
        for (int i = 1; i < N; ++i) {
            long double left = (i == 1) ? u_curr[i - 1] : u_curr[i - 1];
            long double mid = u_curr[i];
            long double right = (i == N - 1) ? right_boundary_value(current_time - dt) : u_curr[i + 1];

            d[i - 1] = (lambda / 2.0L) * left
                + (1.0L - lambda) * mid
                + (lambda / 2.0L) * right;
        }

        // Dodanie warunków brzegowych do RHS
        d[0] += (lambda / 2.0L) * u_curr[0];
        d[N - 2] += (lambda / 2.0L) * right_boundary_value(current_time);

        solve_LU(A, pivot, d, result);

        u_next[0] = 0.0L;
        u_next[N] = right_boundary_value(current_time);
        for (int i = 1; i < N; ++i)
            u_next[i] = result[i - 1];

        u_curr = u_next;
    }

    x_out = x;
    u_out = u_curr;
}
