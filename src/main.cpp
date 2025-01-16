#include <iostream>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <fstream>
#include <string>
#include <vector>

#define DEBUG 1

// struct parameters {
//     double m_e;
//     double m_p;
//     double m_n;
//     double lambda_e;
//     double lambda_p;
//     double lambda_n;
//     double n_B;
// };

struct parameters {
    double m_e;
    double m_p;
    double m_n;
    double x_p;
    double x_e;
};

double Y_p(double n_B) {
    double n_B_threshold = 7.37e29; // Approximate neutron appearance density
    return 0.5 / (1 + exp((log10(n_B) - log10(n_B_threshold))));
}


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


// double f(double x, void* p) {
//     struct parameters* params = (struct parameters*) p;
//     double first_term = params->m_e * pow(1.0 + pow(x, 2.0) * ( pow(params->lambda_e, 2.0) / pow(params->lambda_p, 2.0) ), 1.0/2.0);
//     double second_term = params->m_p * pow(1.0 + pow(x, 2.0), 1.0 / 2.0);
//     double third_term;
//     double n_n, x_n;

//     if (params->n_B <= 7.37e30) {
//         third_term = params->m_n;
//         n_n = 0.0;
//     } else {
//         double second_factor_in_third_term = pow((1.0 + pow(abs(pow(params->lambda_n, 3.0) * params->n_B - ( (pow(params->lambda_n, 3.0) * pow(x, 3.0)) / (pow(params->lambda_p, 3.0)) )), 2.0 / 3.0)), 1.0 / 2.0);
//         third_term = params->m_n * second_factor_in_third_term;
//     }

//     double result = first_term + second_term - third_term;

//     return result;
// }

int main() {
    double m_p = 1.6726e-24; // grams
    double m_n = 1.6750e-24; // grams
    double m_e = 9.1094e-28; // grams
    double hbar = 1.0546e-27; // ergs * seconds
    double c = 2.9979e10; // cm / s
    int status;

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

    outfile << "log_n,log_P,log_rho\n";

    std::vector<double> n_B_list;
    double n_B_min = 1.0e10; // Before neutron appearance
    double n_B_max = 1.0e40; // Well into neutron dominance
    int num_points = 200;    // Higher resolution

    // Logarithmic spacing for smoother transitions
    for (int i = 0; i < num_points; i++) {
        double n_B_i = n_B_min * pow(10, i * (log10(n_B_max / n_B_min) / (num_points - 1)));
        n_B_list.push_back(n_B_i);
    }

    for (int i = 0; i < n_B_list.size(); i++) {
        // double n_B = 7.37e30; // baryons / cm^3 - Neutron Appearence Density
        double n_B = n_B_list[i];
        // double Y_p_current = Y_p(n_B);  // Calculate Y_p for current n_B
        // double n_e = Y_p_current * n_B;
        double n_e = n_B / 2.0; // Initial guess: n_e ≈ n_B / 2
        double x_e = pow(((3 * pow(M_PI, 2.0) * pow(hbar, 3.0)) / pow(m_e * c, 3.0)) * n_e, 1.0 / 3.0);

        double x_p = x_e * pow(lambda_e / lambda_p, 1.0 / 3.0);

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

        // double x0 = pow(3.0 * pow(M_PI, 2.0) * pow(lambda_n, 3.0) * n_B, 1.0 / 3.0);
        double x0 = 3.0;

        gsl_root_fdfsolver_set(s, &F, x0);
    
        int iter = 0, max_iter = 200;
        double x;

        do {
            iter++;
            status = gsl_root_fdfsolver_iterate(s);  // Perform an iteration
            x = gsl_root_fdfsolver_root(s);         // Current root estimate

            // Convergence test based on the change in x
            status = gsl_root_test_delta(x, x0, 0, 1e-6);

            if (DEBUG) {
                std::cout << "Iter " << iter << ": x = " << x << ", f(x) = " << f(x, &params) << std::endl;
            }

            x0 = x;  // Update the previous value

        } while (status == GSL_CONTINUE && iter < max_iter);

        if (status != GSL_SUCCESS) {
            std::cerr << "Newton-Raphson failed for n_B = " << n_B << ". Retrying with a different x0..." << std::endl;
            x0 *= 0.5;  // Adjust initial guess and retry
            gsl_root_fdfsolver_set(s, &F, x0);
        }

        if (status != GSL_CONTINUE && iter >= max_iter) {
            std::cerr << "Exceeded max iterations. Skipping n_B = " << n_B << "." << std::endl;
            break;
        }

        double x_n = x;

        // double x_e = lambda_e * (x_p / lambda_p);
        // double x_n;
        // if (n_B <= 7.37e30) {
        //     x_n = 0;
        // } else {
        //     x_n = pow(abs(pow(lambda_n, 3.0) * (n_B - (pow(x_p, 3.0))/(pow(lambda_p, 3.0)))), 1.0 / 3.0);
        // }

        // double n_e = (1.0 / (3.0 * M_PI * pow(lambda_e, 3.0))) * pow(x_e, 3.0);
        double epsilon_e = ((m_e * pow(c, 2.0)) / (pow(lambda_e, 3.0))) * chi(x_e);
        double P_e = ((m_e * pow(c, 2.0)) / (pow(lambda_e, 3.0))) * phi(x_e);

        double n_p = (1.0 / (3.0 * M_PI * pow(lambda_p, 3.0))) * pow(x_p, 3.0);
        double epsilon_p = ((m_p * pow(c, 2.0)) / (pow(lambda_p, 3.0))) * chi(x_p);
        double P_p = ((m_p * pow(c, 2.0)) / (pow(lambda_p, 3.0))) * phi(x_p);

        double n_n = (1.0 / (3.0 * M_PI * pow(lambda_n, 3.0))) * pow(x_n, 3.0);
        double epsilon_n = ((m_n * pow(c, 2.0)) / (pow(lambda_n, 3.0))) * chi(x_n);
        double P_n = ((m_n * pow(c, 2.0)) / (pow(lambda_n, 3.0))) * phi(x_n);

        double P = P_e + P_p + P_n;
        double rho = (epsilon_e + epsilon_p + epsilon_n) / (pow(c, 2.0));
        double n = n_e + n_p + n_n;

        std::cout << "n_e: " << n_e << std::endl;
        std::cout << "epsilon_e: " << epsilon_e << std::endl;
        std::cout << "P_e: " << P_e << std::endl;

        std::cout << "n_p: " << n_p << std::endl;
        std::cout << "epsilon_p: " << epsilon_p << std::endl;
        std::cout << "P_p: " << P_p << std::endl;

        std::cout << "n_n: " << n_n << std::endl;
        std::cout << "epsilon_n: " << epsilon_n << std::endl;
        std::cout << "P_n: " << P_n << std::endl;

        std::cout << "P: " << P << std::endl;
        std::cout << "rho: " << rho << std::endl;
        std::cout << "n: " << n << std::endl;

        outfile << log10(n) << "," << log10(P) << "," << log10(rho) << "\n";
    }

    outfile.close();

    gsl_root_fdfsolver_free(s);

    return 0;
}