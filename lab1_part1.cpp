#include <iostream>
#include <gmpxx.h>

mpz_class factorial(int n) {
    mpz_class result = 1;
    for (int i = 2; i <= n; i++) {
        result *= i;
    }
    return result;
}

void taylor_series_sin(mpq_class coeffs[], int order) {
    for (int k = 0; k <= order; k++) {
        if (k % 2 == 1) {
            mpz_class fact = factorial(k);
            mpq_class term(1, fact);
            if ((k / 2) % 2 == 1) {
                term = -term;
            }
            coeffs[k] = term;
        } else {
            coeffs[k] = 0;
        }
    }
}

int main() {
    int order = 6; 
    mpq_class coeffs[order + 1];

    taylor_series_sin(coeffs, order);

    std::cout << "sin(x) = ";
    for (int i = 0; i <= order; i++) {
        if (coeffs[i].get_num() != 0){
            std::cout << " + x^" << i << " * " << coeffs[i].get_num() << "/" << coeffs[i].get_den();
        }
    }
    std::cout << '\n';
    return 0;
}
