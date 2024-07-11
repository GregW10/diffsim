#include <iostream>
#include <random>
#include <cstdlib>
#include <cstdint>

void check_numeric(const char *str) {
    while (*str) {
        if (*str < 48 || *str > 57) {
            fprintf(stderr, "Error: non-numeric character '%c' found.\n", *str);
            exit(1);
        }
        ++str;
    }
}

std::pair<long double, long double> lognormal_to_gauss(const long double &mu_ln, const long double &sigma_ln) {
    long double sig_over_mu = sigma_ln/mu_ln;
    long double squared_p1 = 1 + sig_over_mu*sig_over_mu;
    return {logl(mu_ln/sqrtl(squared_p1)), sqrtl(logl(squared_p1))};
}

int main(int argc, char **argv) {
    if (argc < 3 || argc > 4) {
        fprintf(stderr, "Error: invalid number of command-line arguments.\nUsage: ./logn <mu> <sigma> [prec]\n");
        return 1;
    }
    char *endptr;
    long double mu_ln = strtold(*(argv + 1), &endptr);
    if (endptr == *(argv + 1)) {
        fprintf(stderr, "Error: argument for mu could not be converted, \"%s\".\n", *(argv + 1));
        return 1;
    }
    if (mu_ln <= 0) {
        fprintf(stderr, "Error: lognormal mean must be positive.\n");
        return 1;
    }
    long double sig_ln = strtold(*(argv + 2), &endptr);
    if (endptr == *(argv + 2)) {
        fprintf(stderr, "Error: argument for mu could not be converted, \"%s\".\n", *(argv + 1));
        return 1;
    }
    if (sig_ln < 0) {
        fprintf(stderr, "Error: lognormal standard deviation must be non-negative.\n");
        return 1;
    }
    const char *prec = "20";
    if (argc == 4) {
        check_numeric(*(argv + 3));
        prec = *(argv + 3);
    }
    std::string format = "%.";
    format += prec;
    format += "Lf\n";
    std::mt19937_64 rng{std::random_device{}()};
    auto [mu, sig] = lognormal_to_gauss(mu_ln, sig_ln);
    std::lognormal_distribution<long double> dist{mu, sig};
    printf(format.c_str(), dist(rng));
    return 0;
}
