#ifndef SOLVER_UTILS_HPP
#define SOLVER_UTILS_HPP

#include <vector>
#include <string>

// Parametry fizyczne
const long double r = 1.0L;
const long double a = 10.0L;
const long double D = 1.0L;

// Funkcje narzêdziowe
long double exact_solution(long double x, long double t);
long double right_boundary_value(long double t);
void save_vector_to_file(const std::string& filename, const std::vector<long double>& x, const std::vector<long double>& u);

#endif
