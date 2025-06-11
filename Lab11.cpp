#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "solver_utils.hpp"
#include "solver_laasonen.hpp"
#include "solver_crank.hpp"

void test_blad_vs_h();
void test_rozwiazania_w_czasie();
void test_blad_vs_t();
void test_blad_vs_h_osobno();

int main() {
    int N = 100;                  // liczba podziałów przestrzennych
    int M = 1000;                 // liczba kroków czasowych
    long double t_max = 2.0L;     // maksymalny czas


    std::vector<long double> x, u_laasonen_thomas, u_laasonen_lu;
    std::vector<long double> u_crank_thomas, u_crank_lu, u_exact;

    // Laasonen - Thomas
    solve_laasonen(N, M, t_max, x, u_laasonen_thomas);
    save_vector_to_file("output/laasonen_thomas.txt", x, u_laasonen_thomas);

    // Laasonen - LU
    solve_laasonen_LU(N, M, t_max, x, u_laasonen_lu);
    save_vector_to_file("output/laasonen_lu.txt", x, u_laasonen_lu);

    // Crank-Nicolson - Thomas
    solve_crank_nicolson(N, M, t_max, x, u_crank_thomas);
    save_vector_to_file("output/crank_thomas.txt", x, u_crank_thomas);

    // Crank-Nicolson - LU
    solve_crank_LU(N, M, t_max, x, u_crank_lu);
    save_vector_to_file("output/crank_lu.txt", x, u_crank_lu);

    // Exact solution
    for (int i = 0; i <= N; ++i) {
        u_exact.push_back(exact_solution(x[i], t_max));
    }
    save_vector_to_file("output/exact_solution.txt", x, u_exact);

    // Porównanie błędów
    std::ofstream error_file("output/errors_summary.txt");
    error_file << std::fixed << std::setprecision(12);
    error_file << "# x error_laasonen_thomas error_laasonen_lu error_crank_thomas error_crank_lu\n";

    long double max_err_LT = 0.0, max_err_LLU = 0.0, max_err_CT = 0.0, max_err_CLU = 0.0;

    for (int i = 0; i <= N; ++i) {
        long double e_LT = std::abs(u_laasonen_thomas[i] - u_exact[i]);
        long double e_LLU = std::abs(u_laasonen_lu[i] - u_exact[i]);
        long double e_CT = std::abs(u_crank_thomas[i] - u_exact[i]);
        long double e_CLU = std::abs(u_crank_lu[i] - u_exact[i]);

        error_file << x[i] << " " << e_LT << " " << e_LLU << " " << e_CT << " " << e_CLU << "\n";

        max_err_LT = std::max(max_err_LT, e_LT);
        max_err_LLU = std::max(max_err_LLU, e_LLU);
        max_err_CT = std::max(max_err_CT, e_CT);
        max_err_CLU = std::max(max_err_CLU, e_CLU);
    }

    error_file.close();

    ///test_rozwiazania_w_czasie();
    //test_blad_vs_t();
    test_blad_vs_h_osobno();


    // Podsumowanie na ekran
    std::cout << "MAX ERRORS at t = " << t_max << ":\n";
    std::cout << "  Laasonen - Thomas:        " << max_err_LT << "\n";
    std::cout << "  Laasonen - LU:            " << max_err_LLU << "\n";
    std::cout << "  Crank-Nicolson - Thomas:  " << max_err_CT << "\n";
    std::cout << "  Crank-Nicolson - LU:      " << max_err_CLU << "\n";

    return 0;
}

void test_rozwiazania_w_czasie() {
    std::cout << "Uruchomiono test: Rozwiązania w wybranych chwilach czasu...\n";
    int N = 100;
    int M = 1000;
    long double t_max = 2.0L;
    long double dt = t_max / M;

    // Wybrane momenty czasu do obserwacji
    std::vector<long double> czasy_obserwacji = { 0.2L, 0.5L, 1.0L, 1.5L, 2.0L };

    for (long double t_obs : czasy_obserwacji) {
        int m_obs = static_cast<int>(round(t_obs / dt)); // numer kroku czasowego do obserwacji

        std::vector<long double> x, u_lt, u_lu, u_ct, u_clu, u_exact;

        solve_laasonen(N, m_obs, t_obs, x, u_lt);
        solve_laasonen_LU(N, m_obs, t_obs, x, u_lu);
        solve_crank_nicolson(N, m_obs, t_obs, x, u_ct);
        solve_crank_LU(N, m_obs, t_obs, x, u_clu);

        for (int i = 0; i <= N; ++i)
            u_exact.push_back(exact_solution(x[i], t_obs));

        // Zapis do plików
        std::string suffix = "_t" + std::to_string(static_cast<int>(t_obs * 10));
        save_vector_to_file("output/laasonen_thomas" + suffix + ".txt", x, u_lt);
        save_vector_to_file("output/laasonen_lu" + suffix + ".txt", x, u_lu);
        save_vector_to_file("output/crank_thomas" + suffix + ".txt", x, u_ct);
        save_vector_to_file("output/crank_lu" + suffix + ".txt", x, u_clu);
        save_vector_to_file("output/exact_solution" + suffix + ".txt", x, u_exact);

        std::cout << "t = " << t_obs << "s -> Zapisano profile rozwiązań.\n";
    }
    std::cout << "✔️ Zakończono. Dane do zadania 2 zapisano w plikach output/*_t*.txt\n";
}

void test_blad_vs_t() {
    std::cout << "Uruchomiono test: Błąd w funkcji czasu...\n";
    std::ofstream fout("output/blad_vs_t.txt");
    fout << std::fixed << std::setprecision(10);
    fout << "# t max_err_laasonen_thomas max_err_laasonen_lu max_err_crank_thomas max_err_crank_lu\n";

    int N = 100;
    int M = 1000;
    long double t_max = 2.0L;
    long double dt = t_max / M;

    for (int step = 1; step <= M; step += 50) {
        long double t = dt * step;
        std::vector<long double> x, u_lt, u_lu, u_ct, u_clu, u_exact;

        solve_laasonen(N, step, t, x, u_lt);
        solve_laasonen_LU(N, step, t, x, u_lu);
        solve_crank_nicolson(N, step, t, x, u_ct);
        solve_crank_LU(N, step, t, x, u_clu);

        for (int i = 0; i <= N; ++i)
            u_exact.push_back(exact_solution(x[i], t));

        long double max_lt = 0, max_lu = 0, max_ct = 0, max_clu = 0;
        for (int i = 0; i <= N; ++i) {
            max_lt = std::max(max_lt, std::abs(u_lt[i] - u_exact[i]));
            max_lu = std::max(max_lu, std::abs(u_lu[i] - u_exact[i]));
            max_ct = std::max(max_ct, std::abs(u_ct[i] - u_exact[i]));
            max_clu = std::max(max_clu, std::abs(u_clu[i] - u_exact[i]));
        }

        fout << t << " " << max_lt << " " << max_lu << " " << max_ct << " " << max_clu << "\n";

        if (step % 200 == 1) std::cout << "t = " << t << "s -> Obliczono błąd.\n";
    }

    fout.close();
    std::cout << "✔️ Zakończono. Dane do zadania 3 zapisano w pliku output/blad_vs_t.txt\n";
}

void test_blad_vs_h_osobno() {
    std::cout << "Uruchomiono test: Błąd w funkcji kroku h...\n";
    std::vector<int> N_values = { 20, 40, 80, 160, 320 };
    long double t_max = 2.0L;
    long double lambda = 1.0L; // Lambda dla metod niejawnych 

    std::ofstream f_lt("output/blad_vs_h_laasonen_thomas.txt");
    std::ofstream f_lu("output/blad_vs_h_laasonen_lu.txt");
    std::ofstream f_ct("output/blad_vs_h_crank_thomas.txt");
    std::ofstream f_clu("output/blad_vs_h_crank_lu.txt");

    f_lt << "# h blad_max\n";
    f_lu << "# h blad_max\n";
    f_ct << "# h blad_max\n";
    f_clu << "# h blad_max\n";

    for (int N : N_values) {
        long double h = a / N;
        // Obliczamy M tak, aby zachować stałe lambda
        int M = static_cast<int>(round(D * t_max / (lambda * h * h)));

        std::vector<long double> x, u_lt, u_lu, u_ct, u_clu, u_exact;

        solve_laasonen(N, M, t_max, x, u_lt);
        solve_laasonen_LU(N, M, t_max, x, u_lu);
        solve_crank_nicolson(N, M, t_max, x, u_ct);
        solve_crank_LU(N, M, t_max, x, u_clu);

        for (int i = 0; i <= N; ++i)
            u_exact.push_back(exact_solution(x[i], t_max));

        long double max_lt = 0.0L, max_lu = 0.0L, max_ct = 0.0L, max_clu = 0.0L;

        for (int i = 0; i <= N; ++i) {
            max_lt = std::max(max_lt, std::abs(u_lt[i] - u_exact[i]));
            max_lu = std::max(max_lu, std::abs(u_lu[i] - u_exact[i]));
            max_ct = std::max(max_ct, std::abs(u_ct[i] - u_exact[i]));
            max_clu = std::max(max_clu, std::abs(u_clu[i] - u_exact[i]));
        }

        f_lt << std::fixed << std::setprecision(15) << h << " " << max_lt << "\n";
        f_lu << std::fixed << std::setprecision(15) << h << " " << max_lu << "\n";
        f_ct << std::fixed << std::setprecision(15) << h << " " << max_ct << "\n";
        f_clu << std::fixed << std::setprecision(15) << h << " " << max_clu << "\n";

        std::cout << "N = " << N << ", M = " << M << " -> Obliczono.\n";
    }

    f_lt.close(); f_lu.close(); f_ct.close(); f_clu.close();
    std::cout << "✔️ Zakończono. Dane do zadania 1 zapisano w plikach output/blad_vs_h_*.txt\n";
}


