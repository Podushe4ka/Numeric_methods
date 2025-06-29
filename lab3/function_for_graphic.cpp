#include <cmath>
#include <iostream>

double dragCoefficient(double M) {
    const double a = 0.560;
    const double b = 0.5;

    if (M <= 1.1) {
        exit(-1);
    }
    return a / std::pow(M*M - 1.0, b);
}

int main() {
    double Mach[] = {1.2, 1.5, 2.0, 2.5, 3.0};
    for (double M : Mach) {
        std::cout << "M=" << M
                  << " -> cD=" << dragCoefficient(M)
                  << std::endl;
    }
    return 0;
}
