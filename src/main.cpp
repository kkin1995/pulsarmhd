#include <iostream>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

#define DEBUG 1

struct parameters {
    double m_e;
    double m_p;
    double m_n;
    double x_p;
    double x_e;
};

double chi(double x) {
    double first_factor = 1 / (8.0 * pow(M_PI, 2.0));
    double second_factor = x * pow((1 + pow(x, 2.0)), 1.0/2.0) * (1.0 + 2.0 * pow(x, 2.0)) - std::log( x + pow(( 1.0 + pow(x, 2.0) ), 1.0 / 2.0) );
    
    return first_factor * second_factor;
}

double phi(double x) {
    double first_factor = 1 / (8.0 * pow(M_PI, 2.0));
    double second_factor = x * (pow(1.0 + pow(x, 2.0), 1.0 / 2.0)) * ((2.0 / 3.0) * pow(x, 2.0) - 1.0) + std::log(x + pow(1.0 + pow(x, 2.0), 1.0 / 2.0));

    return first_factor * second_factor;
}

double f(double x_n, void* p) {
    struct parameters* params = (struct parameters*) p;
    double first_term = pow(1.0 + pow(params->x_e, 2.0), 1.0 / 2.0) * params->m_e;
    double second_term = pow(1.0 + pow(params->x_p, 2.0), 1.0 / 2.0) * params->m_p;
    double third_term = pow(1.0 + pow(x_n, 2.0), 1.0 / 2.0) * params->m_n;

    return first_term + second_term - third_term;
}

double f_prime(double x_n, void* p) {
    struct parameters* params = (struct parameters*) p;
    return -params->m_n * x_n / sqrt(1.0 + pow(x_n, 2.0));
}

void f_and_df(double x_n, void* p, double* f_value, double* df_value) {
    *f_value = f(x_n, p);
    *df_value = f_prime(x_n, p);
}

int main() {
    double m_p = 1.6726e-24; // grams
    double m_n = 1.6750e-24; // grams
    double m_e = 9.1094e-28; // grams
    double hbar = 1.0546e-27; // ergs * seconds
    double c = 2.9979e10; // cm / s
    double n_B_threshold = 7.37e29; // Approximate neutron appearance density
    double x_n, x_e, x_p;
    int status;
    bool converged;

    const gsl_root_fdfsolver_type* T = gsl_root_fdfsolver_newton;
    gsl_root_fdfsolver* s = gsl_root_fdfsolver_alloc(T);

    double lambda_p =  hbar / (m_p * c);
    double lambda_n =  hbar / (m_n * c);
    double lambda_e =  hbar / (m_e * c);

    std::string filename = "ideal_non_magnetic_npe_gas_eos.csv";
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
    }

    outfile << "log_n,log_P,log_rho,converged\n";

    std::vector<double> n_B_list;
    double n_B_min = 1.0e28; // Before neutron appearance
    double n_B_max = 1.0e32; // Well into neutron dominance
    int num_points = 100;    // Higher resolution

    // Logarithmic spacing for smoother transitions
    for (int i = 0; i < num_points; i++) {
        double n_B_i = n_B_min * pow(10, i * (log10(n_B_max / n_B_min) / (num_points - 1)));
        n_B_list.push_back(n_B_i);
    }

    double prev_x_e = 0.1;
    double prev_x_n = 0.1;

    double n_B_lower = 0.8 * n_B_threshold;
    double n_B_upper = 1.2 * n_B_threshold;

    for (std::size_t i = 0; i < n_B_list.size(); i++) {
        double n_B = n_B_list[i];

        std::cout << "n_B: " << n_B << std::endl;

        if (n_B < n_B_lower) {
            x_n = 0.0;
            x_p = pow((3.0 * pow(M_PI, 2.0) * pow(lambda_p, 3.0) * n_B), 1.0 / 3.0);
            x_e = x_p * lambda_e / lambda_p;
            converged = 1;
        } else if (n_B >= n_B_lower && n_B <= n_B_upper) {
            // double S = smoothstep(n_B, n_B_lower, n_B_upper);
            double excess = sqrt((n_B - n_B_threshold) / n_B_threshold); 
            x_e = pow((3.0 * pow(M_PI, 2.0) * pow(lambda_e, 3.0) * n_B), 1.0 / 3.0);
            x_p = x_e * lambda_p / lambda_e;
            x_n = 0.01 * excess;
            converged = 1;
        } else {
            x_e = prev_x_e;
            x_n = pow((3.0 * pow(M_PI, 2.0) * pow(lambda_n, 3.0) * n_B), 1.0 / 3.0);
            x_p = x_e * lambda_p / lambda_e;

            parameters params = {
                .m_e = m_e,
                .m_p = m_p,
                .m_n = m_n,
                .x_p = x_p,
                .x_e = x_e
            };

            gsl_function_fdf F;
            F.f = &f;          // Function
            F.df = &f_prime;   // Derivative
            F.fdf = &f_and_df; // Combined
            F.params = &params; // Parameters

            double x0 = x_n * 1.1;

            gsl_root_fdfsolver_set(s, &F, x0);
        
            int iter = 0, max_iter = 200;
            double x;

            do {
                iter++;
                double f_value = f(x0, &params);
                double f_prime_value = f_prime(x0, &params);

                status = gsl_root_fdfsolver_iterate(s);  // Perform an iteration
                x = gsl_root_fdfsolver_root(s);         // Current root estimate

                // Convergence test based on the change in x
                status = gsl_root_test_delta(x, x0, 0, 1e-6);

                if (DEBUG) {
                    if (iter % 10 == 0) {
                        std::cout << "Iter " << iter << ": x = " << x << ", f(x) = " << f(x, &params) << std::endl;
                    }
                }

                x0 = x;  // Update the previous value

            } while (status == GSL_CONTINUE && iter < max_iter);

            if (status == GSL_SUCCESS) {
                converged = true;
            } else {
                converged = false;
                std::cerr << "Newton-Raphson failed for n_B = " << n_B << ". Skipping..." << std::endl;
            }

            if (status != GSL_CONTINUE && iter >= max_iter) {
                std::cerr << "Exceeded max iterations. Skipping n_B = " << n_B << "." << std::endl;
                break;
            }

            x_n = x;
            prev_x_e = x_e;
            prev_x_n = x_n;

        }

        double n_e = (1.0 / (3.0 * pow(M_PI, 2.0) * pow(lambda_e, 3.0))) * pow(x_e, 3.0);
        double epsilon_e = ((m_e * pow(c, 2.0)) / (pow(lambda_e, 3.0))) * chi(x_e);
        double P_e = ((m_e * pow(c, 2.0)) / (pow(lambda_e, 3.0))) * phi(x_e);

        double n_p = (1.0 / (3.0 * pow(M_PI, 2.0) * pow(lambda_p, 3.0))) * pow(x_p, 3.0);
        double epsilon_p = ((m_p * pow(c, 2.0)) / (pow(lambda_p, 3.0))) * chi(x_p);
        double P_p = ((m_p * pow(c, 2.0)) / (pow(lambda_p, 3.0))) * phi(x_p);

        double n_n = (1.0 / (3.0 * pow(M_PI, 2.0) * pow(lambda_n, 3.0))) * pow(x_n, 3.0);
        double epsilon_n = ((m_n * pow(c, 2.0)) / (pow(lambda_n, 3.0))) * chi(x_n);
        double P_n = ((m_n * pow(c, 2.0)) / (pow(lambda_n, 3.0))) * phi(x_n);

        double P = P_e + P_p + P_n;
        double rho = (epsilon_e + epsilon_p + epsilon_n) / (pow(c, 2.0));
        double n = n_e + n_p + n_n;

        if (P < 0 || rho < 0 || n < 0) {
            std::cerr << "Warning: Negative values detected. n_B = " << n_B
                    << ", P = " << P << ", rho = " << rho << ", n = " << n << std::endl;
        }

        // std::cout << "n_e: " << n_e << std::endl;
        // std::cout << "epsilon_e: " << epsilon_e << std::endl;
        // std::cout << "P_e: " << P_e << std::endl;

        // std::cout << "n_p: " << n_p << std::endl;
        // std::cout << "epsilon_p: " << epsilon_p << std::endl;
        // std::cout << "P_p: " << P_p << std::endl;

        // std::cout << "n_n: " << n_n << std::endl;
        // std::cout << "epsilon_n: " << epsilon_n << std::endl;
        // std::cout << "P_n: " << P_n << std::endl;

        std::cout << "P: " << P << std::endl;
        std::cout << "rho: " << rho << std::endl;
        std::cout << "n: " << n << std::endl;

        outfile << std::setprecision(10) << log10(n) << "," << log10(P) << "," << log10(rho) << "," << static_cast<int>(converged) << "\n";
    }

    outfile.close();

    gsl_root_fdfsolver_free(s);

    return 0;
}