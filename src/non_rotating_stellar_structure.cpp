#include "non_rotating_stellar_structure.hpp"
#include "polytropic_eos.hpp"
#include "rk4.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>

std::tuple<int, double> non_rotating_stellar_structure(PolytropicGasType eos_type, double rho_c, double r_start, double r_end, double dlogr, double mu_e) {
    /*
     * STELLAR STRUCTURE INTEGRATION USING TOV EQUATIONS
     * 
     * This function solves the Tolman-Oppenheimer-Volkoff (TOV) equations for relativistic
     * stellar structure in logarithmic coordinates. The approach integrates outward from
     * a small starting radius until the stellar surface is reached.
     * 
     * PHYSICS:
     * - Uses polytropic EOS: P = k * ρ^γ
     * - Solves TOV equations with relativistic corrections
     * - Handles both white dwarf (electron gas) and neutron star (neutron gas) regimes
     * - Supports both non-relativistic (γ=5/3) and relativistic (γ=4/3) cases
     * 
     * NUMERICAL METHOD:
     * - Logarithmic coordinates: log(r), log(m), log(P)
     * - 4th-order Runge-Kutta integration with fixed step size
     * - Surface detection via pressure threshold (P < 10^-8 dyne/cm²)
     * - Stability checks for finite values at each step
     * 
     * INITIAL CONDITIONS:
     * - Assumes uniform density sphere: m₀ = (4π/3) * r_start³ * ρ_c
     * - Initial pressure from EOS: P₀ = k * ρ_c^γ
     * - Custom μₑ support for electron gases (composition effects)
     */
    
    // NEW: Get EOS parameters from polytropic calculator
    PolytropicEOS eos_calculator;
    auto eos_data = eos_calculator.getEOSParameters(eos_type);
    double k = eos_data.k;
    double gamma = eos_data.gamma;
    std::string name = eos_data.name;
    
    // Adjust k for custom mu_e (for electron gases only)
    if (eos_type == PolytropicGasType::ELECTRON_NON_RELATIVISTIC && std::abs(mu_e - 2.0) > 1e-10) {
        k = 1.0036e13 / std::pow(mu_e, 5.0/3.0);
    } else if (eos_type == PolytropicGasType::ELECTRON_RELATIVISTIC && std::abs(mu_e - 2.0) > 1e-10) {
        k = 1.2435e15 / std::pow(mu_e, 4.0/3.0);
    }
    
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

    // Add check for m approaching zero
    if (m < 1e-30) {
        m = 1e-30;  // Prevent division by zero
    }

    double dlogm_dlogr = ((4.0 * M_PI * pow(r, 3.0) * rho) / m);
    double first_factor = (- (G * m * rho) / (P * r));
    double second_factor = (1.0 + P/(rho * pow(c, 2.0)));
    double third_factor = 1.0 + ( (4.0 * M_PI * P * pow(r, 3.0) ) / ( m * pow(c, 2.0) ) );    
    double fourth_factor = 1.0 / ( 1.0 - ( (2.0 * G * m) / (r * pow(c, 2.0)) ) );

    double dlogP_dlogr = (first_factor * second_factor * third_factor * fourth_factor);

    return {dlogm_dlogr, dlogP_dlogr};
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