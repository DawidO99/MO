#include "thomas_algorithm.hpp"

void thomas_algorithm(const std::vector<long double>& a,
    const std::vector<long double>& b,
    const std::vector<long double>& c,
    const std::vector<long double>& d,
    std::vector<long double>& result)
{
    int n = d.size();
    std::vector<long double> c_prime(n, 0.0L);
    std::vector<long double> d_prime(n, 0.0L);

    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        long double denom = b[i] - a[i] * c_prime[i - 1];
        c_prime[i] = c[i] / denom;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom;
    }

    result[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        result[i] = d_prime[i] - c_prime[i] * result[i + 1];
    }
}
