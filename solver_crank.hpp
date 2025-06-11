#ifndef SOLVER_CRANK_HPP
#define SOLVER_CRANK_HPP

#include <vector>

// Crank-Nicolson z algorytmem Thomasa
void solve_crank_nicolson(int N, int M, long double t_max,
    std::vector<long double>& x_out,
    std::vector<long double>& u_out);

// Crank-Nicolson z pe³n¹ macierz¹ i LU
void solve_crank_LU(int N, int M, long double t_max,
    std::vector<long double>& x_out,
    std::vector<long double>& u_out);

#endif
