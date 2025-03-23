#include <iostream>
#include <gmpxx.h>
#include <vector>
#include <cmath>

mpq_class factorial(int n) {
    mpq_class result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

std::vector<mpq_class> compute_taylor_coefficients(int order) {
    std::vector<mpq_class> coeffs(order + 1, 0);
    for (int k = 1; k <= order; k += 2) {
        coeffs[k] = ((k / 2) % 2 == 0 ? 1 : -1) / factorial(k);
    }
    return coeffs;
}

void approximation(std::vector<mpq_class> &F, std::vector<mpq_class> &Pm, std::vector<mpq_class> &Qn, int m, int n) {
    std::vector<mpq_class> A(m + n + 1, 0);
    std::vector<mpq_class> B(n + 1, 0);
    
    for (int i = 0; i <= m + n; ++i) {
        A[i] = F[i];
    }
    
    Qn.assign(n + 1, 0);
    Qn[0] = 1;
    
    for (int k = 1; k <= n; ++k) {
        Qn[k] = -A[m + k] / A[m];
    }
    
    Pm.assign(m + 1, 0);
    for (int i = 0; i <= m; ++i) {
        Pm[i] = 0;
        for (int j = 0; j <= std::min(i, n); ++j) {
            Pm[i] += F[i - j] * Qn[j];
        }
    }
    if (Qn[0] != 0) {
        mpq_class Q0 = Qn[0];
        for (auto &coef : Pm) {
            coef /= Q0;
        }
        for (auto &coef : Qn) {
            coef /= Q0;
        }
    }
}

void compute_errors(const std::vector<mpq_class> &Pm, const std::vector<mpq_class> &Qn, 
    const std::vector<mpq_class> &Taylor, double a) {
    const int N = 100;
    mpq_class max_err_Poly = 0, max_err_Taylor = 0;

    for (int i = 0; i <= N; ++i) {
        double x = -a + i * (2 * a) / N;

        mpq_class P_val = 0, Q_val = 0;
        for (size_t j = 0; j < Pm.size(); ++j) P_val += Pm[j] * pow(x, j);
        for (size_t j = 0; j < Qn.size(); ++j) Q_val += Qn[j] * pow(x, j);
        mpq_class Poly_approx = P_val / Q_val;
        mpq_class Taylor_approx = 0;
        for (size_t j = 0; j < Taylor.size(); ++j) Taylor_approx += Taylor[j] * pow(x, j);

        mpq_class true_val = sin(x);

        mpq_class err_Poly = abs(Poly_approx.get_d() - true_val);
        mpq_class err_Taylor = abs(Taylor_approx.get_d() - true_val);

        if (err_Poly > max_err_Poly) max_err_Poly = err_Poly;
        if (err_Taylor > max_err_Taylor) max_err_Taylor = err_Taylor;
    }

    std::cout << "Максимальная ошибка Тейлора: " << max_err_Taylor.get_d() << std::endl;
    std::cout << "Максимальная ошибка Полинома: " << max_err_Poly.get_d() << std::endl;

    if (max_err_Poly < max_err_Taylor) {
        std::cout << "Полином приближает лучше!" << std::endl;
    } else {
        std::cout << "Тейлор приближает лучше!" << std::endl;
    }
}


void print_polynomial(const std::vector<mpq_class> &poly, const std::string &name) {
    std::cout << name << "(x) = ";
    bool first = true;
    for (size_t i = 0; i < poly.size(); ++i) {
        if (poly[i] != 0) {
            if (!first) std::cout << " + ";
            std::cout << poly[i] << " * x^" << i;
            first = false;
        }
    }
    std::cout << std::endl;
}

int main() {
    const int m = 5, n = 2;
    const double a = M_PI / 16;

    std::vector<mpq_class> taylor_coeffs = compute_taylor_coefficients(m + n);
    
    std::vector<mpq_class> Pm, Qn;
    approximation(taylor_coeffs, Pm, Qn, m, n);

    std::cout << "Коэффициенты Тейлора для sin(x):\n";
    print_polynomial(taylor_coeffs, "F_Taylor");
    std::cout << "Полином:\n";
    print_polynomial(Pm, "Pm");
    print_polynomial(Qn, "Qn");

    compute_errors(Pm, Qn, taylor_coeffs, a);

    return 0;
}
