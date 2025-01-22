#include "non_magnetic_ideal_npe_gas.hpp"

NonMagneticNPEGas::NonMagneticNPEGas(bool debug) :
    debug_mode(debug),
    // Computed Members Using Already Initialized Values
    lambda_e(hbar / (m_electron * c)),
    lambda_p(hbar / (m_proton * c)),
    lambda_n(hbar / (m_neutron * c))
{}

// Helper Functions
double NonMagneticNPEGas::phi(double x) const {
    if (x < 0.1) {
        return (1.0 / (15 * pow(M_PI, 2.0))) * (pow(x, 5.0) - (5.0 / 14.0) * pow(x, 7.0) + (5.0 / 24.0) * pow(x, 9.0));
    } else if (x > 10.0) {
        return (1.0 / (12 * pow(M_PI, 2.0))) * (pow(x, 4.0) - pow(x, 2.0) + (3.0 / 2.0) * log(2.0 * x));
    } else {
        double first_factor = 1.0 / (8.0 * pow(M_PI, 2.0));
        double second_factor = x * (pow(1.0 + pow(x, 2.0), 1.0 / 2.0)) * ((2.0 / 3.0) * pow(x, 2.0) - 1.0) + std::log(x + pow(1.0 + pow(x, 2.0), 1.0 / 2.0));

        return first_factor * second_factor;
    }
}

double NonMagneticNPEGas::chi(double x) const {
    if (x < 0.1) {
        return (1.0 / (3.0 * pow(M_PI, 2.0))) * (pow(x, 3.0) + (3.0 / 10.0) * pow(x, 5.0) - (3.0 / 56.0) * pow(x, 7.0));
    } else if (x > 10.0) {
        return (1.0 / (4.0 * pow(M_PI, 2.0))) * (pow(x, 4.0) + pow(x, 2.0) - (1.0 / 2.0) * log(2.0 * x));
    } else {
        double first_factor = 1 / (8.0 * pow(M_PI, 2.0));
        double second_factor = x * pow((1 + pow(x, 2.0)), 1.0/2.0) * (1.0 + 2.0 * pow(x, 2.0)) - std::log( x + pow(( 1.0 + pow(x, 2.0) ), 1.0 / 2.0) );

        return first_factor * second_factor;
    }
}

double NonMagneticNPEGas::interpolate_spline(double t, const std::vector<double>& x, const std::vector<double>& y) const {
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

double NonMagneticNPEGas::BetaEquilibriumSolver::f(double x_n, void* p) {
    if (x_n < 0) return 1e100;
    struct Parameters* params = static_cast<Parameters*>(p);
    double first_term = pow(1.0 + pow(x_n, 2.0), 1.0 / 2.0) * params->m_n;
    double second_term = pow(1.0 + pow(params->x_p, 2.0), 1.0 / 2.0) * params->m_p;
    double third_term = pow(1.0 + pow(params->x_e, 2.0), 1.0 / 2.0) * params->m_e;

    return first_term - second_term - third_term;
};

double NonMagneticNPEGas::BetaEquilibriumSolver::f_prime(double x_n, void* p) {
    struct Parameters* params = static_cast<Parameters*>(p);
    return -params->m_n * x_n / sqrt(1.0 + pow(x_n, 2.0));
}

void NonMagneticNPEGas::BetaEquilibriumSolver::f_and_df(double x_n, void* p, double* f_value, double* df_value) {
    *f_value = f(x_n, p);
    *df_value = f_prime(x_n, p);
}

double NonMagneticNPEGas::BetaEquilibriumSolver::solve(double x_n, double n_B, double m_e, double m_p, double m_n, double x_p, double x_e, bool debug) {
    Parameters params = {m_e, m_p, m_n, x_p, x_e};
    
    int status;

    const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
    gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
    gsl_function F;
    F.function = &f;
    F.params = &params;

    double x0 = x_n;

    double x_lo = 0.0;
    double x_hi = x0;
    double f_lo = f(x_lo, &params);
    double f_hi = f(x_hi, &params);

    if (debug) {
        std::cout << "Initial: f(x_lo) = " << f_lo << ", f(x_hi) = " << f_hi << std::endl;
    }

    // If both values have the same sign, we need to search
    if (f_lo * f_hi > 0) {
        if (f_lo > 0) {
            double factor = 0.1;  // Start with smaller steps
            int max_tries = 100;
            int tries = 0;
            while (f_hi > 0 && tries < max_tries) {
                x_hi *= factor;
                f_hi = f(x_hi, &params);
                tries++;
            }
        }
        if (f_lo < 0) {
            double factor = 10.0; // Start with larger steps
            int max_tries = 100;
            int tries = 0;
            while (f_hi < 0 && tries < max_tries) {
                x_hi *= factor;
                f_hi = f(x_hi, &params);
                tries++;
            }
        }
    }

    if (debug) {
        std::cout << "x_lo = " << x_lo << " | f(x_lo) = " << f(x_lo, &params) << std::endl;
        std::cout << "x_hi = " << x_hi << " | f(x_hi) = " << f(x_hi, &params) << std::endl;
    }

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
        double f_value = f(x0, &params);
        double f_prime_value = f_prime(x0, &params);

        status = gsl_root_fsolver_iterate(s);  // Perform an iteration
        x = gsl_root_fsolver_root(s);         // Current root estimate

        // Convergence test based on the change in x
        status = gsl_root_test_delta(x, x0, 0, 1e-6);

        if (debug) {
            if (iter % 10 == 0) {
                std::cout << "Iter " << iter << ": x = " << x << ", f(x) = " << f(x, &params) << std::endl;
            }
        }

        x0 = x;  // Update the previous value

    } while (status == GSL_CONTINUE && iter < max_iter);

    if (status == GSL_SUCCESS && debug) {
        std::cout << "Succeeded for n_B = " << n_B << std::endl;
    } else {
        std::cerr << "Failed for n_B = " << n_B << ". Skipping..." << std::endl;
    }

    if (status != GSL_CONTINUE && iter >= max_iter && debug) {
        std::cerr << "Exceeded max iterations. Skipping n_B = " << n_B << "." << std::endl;
        return 1;
    }

    x_n = x;

    gsl_root_fsolver_free(s);

    return x_n;

}

void NonMagneticNPEGas::calculateEOS(const std::string &filename, double nB_min, double nB_max, int num_points) {
    double x_n, x_e, x_p;
    double last_low_xe, last_low_xp, last_low_xn;  // Store last point before transition
    double first_high_xe, first_high_xp, first_high_xn;  // Store first point after transition
    bool found_low_boundary = false;
    bool found_high_boundary = false;

    //std::string filename = "ideal_non_magnetic_npe_gas_eos.csv";
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
    }

    outfile << "log_n,log_P,log_rho\n";

    std::vector<double> n_B_list;
    // double n_B_min = 1.0e27; // Before neutron appearance
    // double n_B_max = 1.0e35; // Well into neutron dominance
    // int num_points = 300;    // Higher resolution

    // Logarithmic spacing for smoother transitions
    for (int i = 0; i < num_points; i++) {
        double n_B_i = nB_min * pow(10, i * (log10(nB_max / nB_min) / (num_points - 1)));
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
            x_n = 0.0;
            x_e = lambda_e * pow(3.0 * pow(M_PI, 2.0) * n_B, 1.0 / 3.0);
            x_p = (m_electron / m_proton) * x_e;

            // Store last points before transition
            last_low_xe = x_e;
            last_low_xp = x_p;
            last_low_xn = x_n;
            found_low_boundary = true;

            if (debug_mode) {
                std::cout << "Electron Rich Gas" << std::endl;
                std::cout << "x_n = " << x_n << std::endl;
                std::cout << "x_e = " << x_e << std::endl;
                std::cout << "x_p = " << x_p << std::endl;
            }
        }
        // 3. Transition Region
        else if (n_B >= n_B_lower && n_B <= n_B_upper) {
            if (!found_high_boundary) {
                // Set high boundary
                double Q = m_neutron - m_proton;
                double xe_min = sqrt(pow(Q / m_electron, 2.0) - 1.0);
                x_e = xe_min * pow(n_B_upper / n_B_threshold, 1.0 / 3.0);
                x_p = (m_electron / m_proton) * x_e;
                x_n = pow((3.0 * pow(M_PI, 2.0) * pow(lambda_n, 3.0) * n_B_upper) / (1.0 + 0.0026), 1.0 / 3.0);

                x_n = BetaEquilibriumSolver::solve(x_n, n_B_upper, m_electron, m_proton, m_neutron, x_p, x_e, debug_mode);
                if (x_n == 1) {
                    std::cerr << "Fallback: High boundary values for n_B = " << n_B_upper << std::endl;
                    x_e = xe_min * pow(n_B_upper / n_B_threshold, 1.0 / 3.0);
                    x_p = (m_electron / m_proton) * x_e;
                    x_n = pow((3.0 * pow(M_PI, 2.0) * pow(lambda_n, 3.0) * n_B_upper) / (1.0 + 0.0026), 1.0 / 3.0);
                }

                first_high_xe = x_e;
                first_high_xp = x_p;
                first_high_xn = x_n;
                found_high_boundary = true;
            }

            if (found_low_boundary && found_high_boundary) {
                if (debug_mode) {
                    std::cout << "Low boundary: x_e = " << last_low_xe << ", x_p = " << last_low_xp << ", x_n = " << last_low_xn << std::endl;
                    std::cout << "High boundary: x_e = " << first_high_xe << ", x_p = " << first_high_xp << ", x_n = " << first_high_xn << std::endl;
                }
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
            if (debug_mode) {
                std::cout << "Transition Region - Appearence of Neutrons" << std::endl;
                std::cout << "x_n = " << x_n << std::endl;
                std::cout << "x_e = " << x_e << std::endl;
                std::cout << "x_p = " << x_p << std::endl;
            }
        }
        // 4. High Density Regime
        else {
            double Q = m_neutron - m_proton;
            double xe_min = sqrt(pow(Q/m_electron, 2.0) - 1.0);
            x_e = xe_min * pow(n_B/n_B_threshold, 1.0/3.0);
            x_p = (m_electron / m_proton) * x_e;
            x_n = pow((3.0 * pow(M_PI, 2.0) * pow(lambda_n, 3.0) * n_B) / (1.0 + 0.0026), 1.0/3.0);

            x_n = BetaEquilibriumSolver::solve(x_n, n_B, m_electron, m_proton, m_neutron, x_p, x_e, debug_mode);
            if (x_n == 1) {
                std::cerr << "Beta equilibrium solver failed in high density region at n_B = " << n_B << std::endl;
                continue;
            }
            
            if (debug_mode) {
                std::cout << "Neutron Rich Region" << std::endl;
                std::cout << "x_n = " << x_n << std::endl;
                std::cout << "x_e = " << x_e << std::endl;
                std::cout << "x_p = " << x_p << std::endl;
            }
        }

        double n_e = (1.0 / (3.0 * pow(M_PI, 2.0) * pow(lambda_e, 3.0))) * pow(x_e, 3.0);
        double epsilon_e = ((m_electron * pow(c, 2.0)) / (pow(lambda_e, 3.0))) * chi(x_e);
        double P_e = ((m_electron * pow(c, 2.0)) / (pow(lambda_e, 3.0))) * phi(x_e);

        double n_p = (1.0 / (3.0 * pow(M_PI, 2.0) * pow(lambda_p, 3.0))) * pow(x_p, 3.0);
        double epsilon_p = ((m_proton * pow(c, 2.0)) / (pow(lambda_p, 3.0))) * chi(x_p);
        double P_p = ((m_proton * pow(c, 2.0)) / (pow(lambda_p, 3.0))) * phi(x_p);

        double n_n = (1.0 / (3.0 * pow(M_PI, 2.0) * pow(lambda_n, 3.0))) * pow(x_n, 3.0);
        double epsilon_n = ((m_neutron * pow(c, 2.0)) / (pow(lambda_n, 3.0))) * chi(x_n);
        double P_n = ((m_neutron * pow(c, 2.0)) / (pow(lambda_n, 3.0))) * phi(x_n);

        double P = P_e + P_p + P_n;
        double rho = (epsilon_e + epsilon_p + epsilon_n) / (pow(c, 2.0));
        double n = n_p + n_n;

        if (P < 0 || rho < 0 || n < 0) {
            std::cerr << "Warning: Negative values detected. n_B = " << n_B
                    << ", P = " << P << ", rho = " << rho << ", n = " << n << std::endl;
        }

        if (debug_mode) {
            std::cout << "P: " << P << std::endl;
            std::cout << "rho: " << rho << std::endl;
            std::cout << "n: " << n << std::endl;
        }

        outfile << std::setprecision(10) << log10(n) << "," << log10(P) << "," << log10(rho) << "\n";
    }

    outfile.close();
}



