#include <iostream>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
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

double interpolate_spline(double t, const std::vector<double>& x, const std::vector<double>& y) {
    // Ensure input sizes match
    if (x.size() != y.size()) {
        std::cerr << "Error: x and y arrays must have the same size for spline interpolation." << std::endl;
        return NAN;
    }

    // Initialize GSL spline and accelerator
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, x.size());

    // Initialize the spline
    gsl_spline_init(spline, x.data(), y.data(), x.size());

    // Evaluate the spline at t
    double result = gsl_spline_eval(spline, t, acc);

    // Free resources
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return result;
}

double chi(double x) {
    double first_factor = 1 / (8.0 * pow(M_PI, 2.0));
    double second_factor = x * pow((1 + pow(x, 2.0)), 1.0/2.0) * (1.0 + 2.0 * pow(x, 2.0)) - std::log( x + pow(( 1.0 + pow(x, 2.0) ), 1.0 / 2.0) );
    
    return first_factor * second_factor;
}

// double phi(double x) {
//     double first_factor = 1 / (8.0 * pow(M_PI, 2.0));
//     double second_factor = x * (pow(1.0 + pow(x, 2.0), 1.0 / 2.0)) * ((2.0 / 3.0) * pow(x, 2.0) - 1.0) + std::log(x + pow(1.0 + pow(x, 2.0), 1.0 / 2.0));

//     return first_factor * second_factor;
// }

double phi(double x) {
    if (x < 1e-3) {
        // For very small x, use series expansion to avoid numerical issues
        // Leading terms in Taylor series for small x
        return (x*x*x*x) / (20.0 * pow(M_PI, 2.0));  // First non-zero term
    }
    
    // For larger x, use the full expression
    double first_factor = 1 / (8.0 * pow(M_PI, 2.0));
    double second_factor = x * pow(1.0 + pow(x, 2.0), 1.0 / 2.0) * 
                          ((2.0 / 3.0) * pow(x, 2.0) - 1.0) + 
                          std::log(x + pow(1.0 + pow(x, 2.0), 1.0 / 2.0));
    
    // Ensure result is non-negative (physical requirement)
    return std::max(first_factor * second_factor, 0.0);
}

double f(double x_n, void* p) {
    if (x_n < 0) return 1e100;
    struct parameters* params = (struct parameters*) p;
    double first_term = pow(1.0 + pow(x_n, 2.0), 1.0 / 2.0) * params->m_n;
    double second_term = pow(1.0 + pow(params->x_p, 2.0), 1.0 / 2.0) * params->m_p;
    double third_term = pow(1.0 + pow(params->x_e, 2.0), 1.0 / 2.0) * params->m_e;

    return first_term - second_term - third_term;
}

double f_prime(double x_n, void* p) {
    struct parameters* params = (struct parameters*) p;
    return -params->m_n * x_n / sqrt(1.0 + pow(x_n, 2.0));
}

void f_and_df(double x_n, void* p, double* f_value, double* df_value) {
    *f_value = f(x_n, p);
    *df_value = f_prime(x_n, p);
}

double solve_beta_equilibrium(double x_n, double n_B, parameters *params) {
    int status;
    bool converged;

    const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
    gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
    gsl_function F;
    F.function = &f;          // Function
    F.params = params; // Parameters

    double x0 = x_n;

    double x_lo = 0.0;
    double x_hi = x0;
    double f_lo = f(x_lo, params);
    double f_hi = f(x_hi, params);

    std::cout << "Initial: f(x_lo) = " << f_lo << ", f(x_hi) = " << f_hi << std::endl;

    // If both values have the same sign, we need to search
    if (f_lo * f_hi > 0) {
        if (f_lo > 0) {
            double factor = 0.1;  // Start with smaller steps
            int max_tries = 100;  // Allow more iterations
            int tries = 0;
            while (f_hi > 0 && tries < max_tries) {
                x_hi *= factor;
                f_hi = f(x_hi, params);
                tries++;
            }
        }
        if (f_lo < 0) {
            double factor = 10.0; // Start with larger steps
            int max_tries = 100;
            int tries = 0;
            while (f_hi < 0 && tries < max_tries) {
                x_hi *= factor;
                f_hi = f(x_hi, params);
                tries++;
            }
        }
    }

    std::cout << "x_lo = " << x_lo << " | f(x_lo) = " << f(x_lo, params) << std::endl;
    std::cout << "x_hi = " << x_hi << " | f(x_hi) = " << f(x_hi, params) << std::endl;

    // Check if we found straddling points
    if (f_lo * f_hi >= 0) {
        std::cerr << "Could not find straddling points for root." << std::endl;
        gsl_root_fsolver_free(s);
        return 1;
    }

    status = gsl_root_fsolver_set(s, &F, x_lo, x_hi);

    int iter = 0, max_iter = 200;
    double x;

    do {
        iter++;
        double f_value = f(x0, params);
        double f_prime_value = f_prime(x0, params);

        status = gsl_root_fsolver_iterate(s);  // Perform an iteration
        x = gsl_root_fsolver_root(s);         // Current root estimate

        // Convergence test based on the change in x
        status = gsl_root_test_delta(x, x0, 0, 1e-6);

        if (DEBUG) {
            if (iter % 10 == 0) {
                std::cout << "Iter " << iter << ": x = " << x << ", f(x) = " << f(x, params) << std::endl;
            }
        }

        x0 = x;  // Update the previous value

    } while (status == GSL_CONTINUE && iter < max_iter);

    if (status == GSL_SUCCESS) {
        converged = true;
    } else {
        converged = false;
        std::cerr << "Failed for n_B = " << n_B << ". Skipping..." << std::endl;
    }

    if (status != GSL_CONTINUE && iter >= max_iter) {
        std::cerr << "Exceeded max iterations. Skipping n_B = " << n_B << "." << std::endl;
        return 1;
    }

    x_n = x;

    gsl_root_fsolver_free(s);

    return x_n;

}

int main() {
    double m_p = 1.672622e-24; // grams
    double m_n = 1.674927e-24; // grams
    double m_e = 9.1094e-28; // grams
    double hbar = 1.0546e-27; // ergs * seconds
    double c = 2.9979e10; // cm / s
    double n_B_threshold = 7.37e29; // Approximate neutron appearance density
    double x_n, x_e, x_p;
    bool converged;
    double last_low_xe, last_low_xp, last_low_xn;  // Store last point before transition
    double first_high_xe, first_high_xp, first_high_xn;  // Store first point after transition
    bool found_low_boundary = false;
    bool found_high_boundary = false;

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
    double n_B_min = 1.0e27; // Before neutron appearance
    double n_B_max = 1.0e33; // Well into neutron dominance
    int num_points = 300;    // Higher resolution

    // Logarithmic spacing for smoother transitions
    for (int i = 0; i < num_points; i++) {
        double n_B_i = n_B_min * pow(10, i * (log10(n_B_max / n_B_min) / (num_points - 1)));
        n_B_list.push_back(n_B_i);
    }

    double n_B_lower = 0.8 * n_B_threshold;
    double n_B_upper = 1.2 * n_B_threshold;

    for (std::size_t i = 0; i < n_B_list.size(); i++) {
        double n_B = n_B_list[i];

        std::cout << "n_B: " << n_B << std::endl;

        // 1. Low Density Regime
        if (n_B < n_B_lower) {
            // Pure electron-proton gas
            std::cout << "ELECTRON PROTON GAS" << std::endl;
            x_n = 0.0;
            x_e = pow(3.0 * pow(M_PI, 2.0) * pow(lambda_e, 3.0) * n_B, 1.0 / 3.0);
            x_p = (m_e / m_p) * x_e;
            converged = 1;

            // Store last points before transition
            last_low_xe = x_e;
            last_low_xp = x_p;
            last_low_xn = x_n;
            found_low_boundary = true;

            std::cout << "x_n = " << x_n << std::endl;
            std::cout << "x_e = " << x_e << std::endl;
            std::cout << "x_p = " << x_p << std::endl;
        }
        // 3. Transition Region
        else if (n_B >= n_B_lower && n_B <= n_B_upper) {
            std::cout << "TRANSITION REGION - APPEARENCE OF NEUTRONS" << std::endl;
            if (!found_high_boundary) {
                // Set high boundary
                double Q = m_n - m_p;
                double xe_min = sqrt(pow(Q / m_e, 2.0) - 1.0);
                x_e = xe_min * pow(n_B_upper / n_B_threshold, 1.0 / 3.0);
                x_p = (m_e / m_p) * x_e;
                x_n = pow((3.0 * pow(M_PI, 2.0) * pow(lambda_n, 3.0) * n_B_upper) / (1.0 + 0.0026), 1.0 / 3.0);

                parameters params = {m_e, m_p, m_n, x_p, x_e};
                x_n = solve_beta_equilibrium(x_n, n_B_upper, &params);
                if (x_n == 1) {
                    std::cerr << "Fallback: High boundary values for n_B = " << n_B_upper << std::endl;
                    x_e = xe_min * pow(n_B_upper / n_B_threshold, 1.0 / 3.0);
                    x_p = (m_e / m_p) * x_e;
                    x_n = pow((3.0 * pow(M_PI, 2.0) * pow(lambda_n, 3.0) * n_B_upper) / (1.0 + 0.0026), 1.0 / 3.0);
                }

                first_high_xe = x_e;
                first_high_xp = x_p;
                first_high_xn = x_n;
                found_high_boundary = true;
            }

            if (found_low_boundary && found_high_boundary) {
                std::cout << "Low boundary: x_e = " << last_low_xe << ", x_p = " << last_low_xp << ", x_n = " << last_low_xn << std::endl;
                std::cout << "High boundary: x_e = " << first_high_xe << ", x_p = " << first_high_xp << ", x_n = " << first_high_xn << std::endl;
                // Interpolate in the transition region
                std::vector<double> n_B_points, x_e_values, x_p_values, x_n_values;
                int num_spline_points = 300;
                for (int i = 0; i <= num_spline_points; i++) {
                    double t = i / static_cast<double>(num_spline_points);
                    double n_B_i = n_B_lower + t * (n_B_upper - n_B_lower);
                    n_B_points.push_back(n_B_i);

                    x_e_values.push_back(last_low_xe + t * (first_high_xe - last_low_xe));
                    x_p_values.push_back(last_low_xp + t * (first_high_xp - last_low_xp));
                    x_n_values.push_back(last_low_xn + t * (first_high_xn - last_low_xn));
                }

                x_e = interpolate_spline(n_B, n_B_points, x_e_values);
                x_p = interpolate_spline(n_B, n_B_points, x_p_values);
                x_n = interpolate_spline(n_B, n_B_points, x_n_values);
            }
            std::cout << "x_n = " << x_n << std::endl;
            std::cout << "x_e = " << x_e << std::endl;
            std::cout << "x_p = " << x_p << std::endl;
        }
        // 4. High Density Regime
        else {
            std::cout << "NEUTRON RICH GAS" << std::endl;
            double Q = m_n - m_p;
            double xe_min = sqrt(pow(Q/m_e, 2.0) - 1.0);
            x_e = xe_min * pow(n_B/n_B_threshold, 1.0/3.0);
            x_p = (m_e / m_p) * x_e;
            x_n = pow((3.0 * pow(M_PI, 2.0) * pow(lambda_n, 3.0) * n_B) / (1.0 + 0.0026), 1.0/3.0);

            parameters params;
            params.m_e = m_e;
            params.m_p = m_p;
            params.m_n = m_n;
            params.x_p = x_p;
            params.x_e = x_e;

            x_n = solve_beta_equilibrium(x_n, n_B, &params);
            if (x_n == 1) {
                std::cerr << "Beta equilibrium solver failed in high density region at n_B = " << n_B << std::endl;
                continue;
            }
            
            std::cout << "x_n = " << x_n << std::endl;
            std::cout << "x_e = " << x_e << std::endl;
            std::cout << "x_p = " << x_p << std::endl;
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
        double n = n_p + n_n;

        if (P < 0 || rho < 0 || n < 0) {
            std::cerr << "Warning: Negative values detected. n_B = " << n_B
                    << ", P = " << P << ", rho = " << rho << ", n = " << n << std::endl;
        }

        std::cout << "P: " << P << std::endl;
        std::cout << "rho: " << rho << std::endl;
        std::cout << "n: " << n << std::endl;

        outfile << std::setprecision(10) << log10(n) << "," << log10(P) << "," << log10(rho) << "," << static_cast<int>(converged) << "\n";
    }

    outfile.close();
    return 0;
}