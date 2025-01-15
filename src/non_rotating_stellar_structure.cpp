#include "non_rotating_stellar_structure.hpp"
#include "rk4.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>

std::tuple<int, double> non_rotating_stellar_structure(std::string name, double rho_c, double r_start, double r_end, double dlogr, double k, double gamma) {
        std::string filename = get_filename(name, rho_c);
        std::ofstream outfile(filename);
        if (!outfile.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
        }
        std::cout << "Writing to file: " << filename << std::endl;

        double fraction = 4.0 / 3.0;
        double log_r_start = log10(r_start);
        double log_r_end = log10(r_end);
        double log_m0 = log10(fraction) + log10(M_PI) + 3.0*log10(r_start) + log10(rho_c);
        double log_p0 = log10(k) + (gamma * log10(rho_c)); 
        std::vector<double> state = {log_m0, log_p0};

        std::cout << "Initial conditions:" << std::endl;
        std::cout << "  log_r_start = " << log_r_start << std::endl;
        std::cout << "  log_m0 = " << log_m0 << std::endl;
        std::cout << "  log_p0 = " << log_p0 << std::endl;
        std::cout << "  k = " << k << ", gamma = " << gamma << std::endl;

        outfile << "log_r[cm],log_m[g],log_P[dyne/cm^2]\n";

        int idx = 0;
        double log_r = log_r_start;

        while (log_r < log_r_end) {
            state = rk4_step(log_r, dlogr, state, 
            [k, gamma](double r, const std::vector<double>& s) {
                return tolman_oppenheimer_volkoff_derivatives(r, s, k, gamma);
            });
            log_r += dlogr;
            outfile << log_r << "," << state[0] << "," << state[1] << "\n";

            // Stop if pressure drops below a threshold
            if (state[1] < log10(1e-8)) {
                std::cout << "Surface reached at log_r: " << log_r << std::endl;
                break;
            }

            if (!std::isfinite(state[0]) || !std::isfinite(state[1])) {
                std::cout << "Non-finite values at log_r: " << log_r << std::endl;
                break;
            }

            if (idx % 100 == 0) { 
                std::cout << "log_r: " << log_r << ", log_m: " << state[0] 
                        << ", log_P: " << state[1] << std::endl;
            }

            idx += 1;
        }
        outfile.close();
        // At end of integration
        return {idx, state[0]};
}

std::vector<double> newtonian(double log_r, const std::vector<double>& state, double k, double gamma) {
    double m = pow(10.0, state[0]);
    double P = pow(10.0, state[1]);
    double r = pow(10.0, log_r);
    double log_rho = (state[1] - log10(k)) / gamma;
    double rho = pow(10.0, log_rho);

    double dlogm_dlogr = ((4.0 * M_PI * pow(r, 3.0) * rho) / m);
    double dlogP_dlogr = (- (G * m * rho) / (P * r));

    // std::cout << "Debug: r: " << r << ", dlogm_dlogr: " << dlogm_dlogr 
    //           << ", dlogP_dlogr: " << dlogP_dlogr << std::endl;

    return {dlogm_dlogr, dlogP_dlogr};
}

std::vector<double> tolman_oppenheimer_volkoff_derivatives(double log_r, const std::vector<double>& state, double k, double gamma) {
    double m = pow(10.0, state[0]);
    double P = pow(10.0, state[1]);
    double r = pow(10.0, log_r);
    double log_rho = (state[1] - log10(k)) / gamma;
    double rho = pow(10.0, log_rho);

    // Factor calculation with detailed debugging
    // std::cout << "\nDetailed factor analysis at r = " << r << " cm:" << std::endl;
    // std::cout << "Pre-factor values:" << std::endl;
    // std::cout << "  G = " << G << std::endl;
    // std::cout << "  m = " << m << std::endl;
    // std::cout << "  rho = " << rho << std::endl;
    // std::cout << "  P = " << P << std::endl;
    // std::cout << "  r = " << r << std::endl;
    // std::cout << "  c = " << c << std::endl;

    // std::cout << "Converted values:" << std::endl;
    // std::cout << "  r = " << r << " cm" << std::endl;
    // std::cout << "  m = " << m << " g" << std::endl;
    // std::cout << "  P = " << P << " dyne/cm²" << std::endl;
    // std::cout << "  rho = " << rho << " g/cm³" << std::endl;

    // Add check for m approaching zero
    if (m < 1e-30) {
        m = 1e-30;  // Prevent division by zero
    }

    // // Add check for approaching Schwarzschild radius
    // if (2.0 * G * m >= r * pow(c, 2.0)) {
    //     // We've hit the Schwarzschild radius
    //     return {0.0, -INFINITY};
    // }

    double dlogm_dlogr = ((4.0 * M_PI * pow(r, 3.0) * rho) / m);
    double first_factor = (- (G * m * rho) / (P * r));
    // std::cout << "\nFirst factor components:" << std::endl;
    // std::cout << "  G*m = " << G*m << std::endl;
    // std::cout << "  G*m*rho = " << G*m*rho << std::endl;
    // std::cout << "  P*r = " << P*r << std::endl;
    // std::cout << "  first_factor = " << first_factor << std::endl;

    double second_factor = (1.0 + P/(rho * pow(c, 2.0)));
    // std::cout << "\nStep by step calculation:" << std::endl;
    // std::cout << "1. Basic gravitational term (-GMρ/Pr) = " << first_factor << std::endl;
    // std::cout << "2. Energy correction terms:" << std::endl;
    // std::cout << "   ρ/P = " << rho/P << std::endl;
    // std::cout << "   1/c² = " << 1.0/pow(c, 2.0) << std::endl;
    // std::cout << "   Combined = " << second_factor << std::endl;


    double third_factor = 1.0 + ( (4.0 * M_PI * P * pow(r, 3.0) ) / ( m * pow(c, 2.0) ) );
    // std::cout << "3. Pressure contribution (1 + 4πr³P/mc²) = " << third_factor << std::endl;
    
    double fourth_factor = 1.0 / ( 1.0 - ( (2.0 * G * m) / (r * pow(c, 2.0)) ) );
    // std::cout << "4. Metric factor 1/(1 - 2GM/rc²) = " << fourth_factor << std::endl;
    // std::cout << "Final dlogP/dlogr = " << first_factor * second_factor * third_factor * fourth_factor << std::endl;

    // std::cout << "\nAll factors:" << std::endl;
    // std::cout << "  first = " << first_factor << std::endl;
    // std::cout << "  second = " << second_factor << std::endl;
    // std::cout << "  third = " << third_factor << std::endl;
    // std::cout << "  fourth = " << fourth_factor << std::endl;

    // // Debug prints
    // std::cout << "Debug at r = " << r << " cm:" << std::endl;
    // std::cout << "  rho = " << rho << " g/cm³" << std::endl;
    // std::cout << "  factors: " << first_factor << ", " << second_factor 
    //           << ", " << third_factor << ", " << fourth_factor << std::endl;

    double dlogP_dlogr = (first_factor * second_factor * third_factor * fourth_factor);

    // std::cout << "Derivatives:" << std::endl;
    // std::cout << "  dlogm_dlogr = " << dlogm_dlogr << std::endl;
    // std::cout << "  dlogP_dlogr = " << dlogP_dlogr << std::endl;

    return {dlogm_dlogr, dlogP_dlogr};
}

void set_eos_parameters(DegenerateGasType type, std::string &name, double &k, double &gamma) {
    double mu_e = 2.0;
    switch(type) {
        case ELECTRON_NON_RELATIVISTIC:
            name = "electron_non_relativistic";
            k = 1.0036e13 / pow(mu_e, 5.0 / 3.0);
            gamma = 5.0 / 3.0;
            break;
        case ELECTRON_RELATIVISTIC:
            name = "electron_relativistic";
            k = 1.2435e15 / pow(mu_e, 4.0 / 3.0);
            gamma = 4.0 / 3.0;
            break;
        case NEUTRON_NON_RELATIVISTIC:
            name = "neutron_non_relativistic";
            k = 5.3802e9;
            gamma = 5.0 / 3.0;
            break;
        case NEUTRON_RELATIVISTIC:
            name = "neutron_relativistic";
            k = 1.2293e15;
            gamma = 4.0 / 3.0;
            break;
    }
}

std::string get_filename(const std::string& name, double rho_c) {
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(2) << rho_c;
    std::string rho_str = ss.str();
    
    // Only replace in the scientific notation part
    std::replace(rho_str.begin(), rho_str.end(), '+', 'p');
    std::replace(rho_str.begin(), rho_str.end(), 'e', 'p');
    
    return "data/" + name + "_rhoc_" + rho_str + ".csv";
}