#include <fstream>
#include <iomanip>
#include <cmath>
#include "solver_utils.hpp"
#include "CALERF.h"

using namespace calerfpack;

long double exact_solution(long double x, long double t) {
    if (t == 0.0L) return 1.0L;
    long double numerator = erfc_LD((x - r) / (2.0L * std::sqrt(D * t)));
    return 1.0L - (r / x) * numerator;
}

long double right_boundary_value(long double t) {
    if (t == 0.0L) return 1.0L;
    long double val = erfc_LD(a / (2.0L * std::sqrt(D * t)));
    return 1.0L - (r / (r + a)) * val;
}

void save_vector_to_file(const std::string& filename, const std::vector<long double>& x, const std::vector<long double>& u) {
    std::ofstream file(filename);
    file << std::fixed << std::setprecision(10);
    for (size_t i = 0; i < x.size(); ++i) {
        file << x[i] << " " << u[i] << "\n";
    }
    file.close();
}
