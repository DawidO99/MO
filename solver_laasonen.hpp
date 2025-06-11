#ifndef SOLVER_LAASONEN_HPP
#define SOLVER_LAASONEN_HPP

#include <vector>

// Metoda Laasonena z Thomasem
void solve_laasonen(int N, int M, long double t_max,
    std::vector<long double>& x_out,
    std::vector<long double>& u_out);

// Metoda Laasonena z pe³n¹ macierz¹ i LU
void solve_laasonen_LU(int N, int M, long double t_max,
    std::vector<long double>& x_out,
    std::vector<long double>& u_out);

#endif
