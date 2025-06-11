#ifndef LU_SOLVER_HPP
#define LU_SOLVER_HPP

#include <vector>

// Wykonuje dekompozycjê LU: A = L * U (macierz A zostaje nadpisana przez L i U)
void decompose_LU(std::vector<std::vector<long double>>& A, std::vector<int>& pivot);

// Rozwi¹zuje uk³ad A * x = b, zak³adaj¹c ¿e A zosta³o zLU-dekomponowane
void solve_LU(const std::vector<std::vector<long double>>& LU,
    const std::vector<int>& pivot,
    const std::vector<long double>& b,
    std::vector<long double>& x);

#endif
