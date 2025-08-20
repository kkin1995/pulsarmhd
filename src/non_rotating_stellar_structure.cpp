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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

std::tuple<int, double, double> non_rotating_stellar_structure(PolytropicGasType eos_type, double rho_c, double r_start, double r_end, double dlogr, double mu_e) {
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
        // Use 1e15 for Neutron Stars
        if (state[1] < log10(1e15)) {  // Reasonable threshold: 10^15 dyne/cm² for surface detection
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
    return {idx, state[0], log_r};
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

std::vector<double> tolman_oppenheimer_volkoff_derivatives_spline(
    double log_r, 
    const std::vector<double>& state, 
    const gsl_spline* spline_inv, 
    gsl_interp_accel* acc_inv, 
    double min_logP, 
    double max_logP
) {
    /*
     * TOV DERIVATIVES WITH SPLINE-BASED EOS
     * 
     * This function computes the same TOV derivatives as the polytropic version,
     * but uses tabulated EOS data via GSL spline interpolation to get ρ(P).
     * 
     * PHYSICS:
     * - Identical TOV equations to polytropic version
     * - Density obtained from spline: log_ρ = f^(-1)(log_P)
     * - Same relativistic corrections and factors
     * 
     * EOS HANDLING:
     * - Pressure clamping to validity range [min_logP, max_logP]
     * - GSL spline evaluation for density lookup
     * - Error handling for extrapolation cases
     */
    
    double m = pow(10.0, state[0]);
    double r = pow(10.0, log_r);
    
    // Add check for m approaching zero (same as polytropic version)
    if (m < 1e-30) {
        m = 1e-30;  // Prevent division by zero
    }
    
    // Clamp pressure to EOS validity range
    double log_P = state[1];
    if (log_P < min_logP) log_P = min_logP;
    if (log_P > max_logP) log_P = max_logP;
    
    // Get density from spline interpolation: log_rho = f^(-1)(log_P)
    double P = std::pow(10.0, log_P);
    double log_rho = gsl_spline_eval(spline_inv, log_P, acc_inv);
    double rho = pow(10.0, log_rho);
    
    // Compute TOV derivatives (identical to polytropic version)
    double dlogm_dlogr = ((4.0 * M_PI * pow(r, 3.0) * rho) / m);
    double first_factor = (- (G * m * rho) / (P * r));
    double second_factor = (1.0 + P/(rho * pow(c, 2.0)));
    double third_factor = 1.0 + ( (4.0 * M_PI * P * pow(r, 3.0) ) / ( m * pow(c, 2.0) ) );    
    double fourth_factor = 1.0 / ( 1.0 - ( (2.0 * G * m) / (r * pow(c, 2.0)) ) );

    double dlogP_dlogr = (first_factor * second_factor * third_factor * fourth_factor);

    return {dlogm_dlogr, dlogP_dlogr};
}

TovResult non_rotating_stellar_structure_spline(
    const gsl_spline* spline_inv,
    gsl_interp_accel* acc_inv,
    double min_logP,
    double max_logP,
    double rho_c,
    double r_start,
    double r_end,
    double base_dlogr,
    bool use_adaptive_stepping,
    double pressure_threshold,
    const std::string& output_filename
) {
    /*
     * STELLAR STRUCTURE INTEGRATION WITH SPLINE-BASED EOS
     * 
     * This function extends the existing TOV solver to work with tabulated EOS
     * data via GSL splines. It provides the same functionality as the polytropic
     * version but supports realistic, complex equations of state.
     * 
     * PHYSICS:
     * - Same TOV equations as polytropic version
     * - Realistic EOS via spline interpolation
     * - Support for complex phase transitions and composition changes
     * 
     * NUMERICAL ENHANCEMENTS:
     * - Optional adaptive step size control
     * - Configurable surface pressure threshold
     * - Optional file output (empty filename = no output)
     * - Enhanced error handling for EOS validity range
     */
    
    // Get forward spline for initial pressure calculation
    // We need to find the forward spline (P from ρ) - this would need to be passed
    // For now, we'll estimate initial pressure using a simple approach
    // In practice, you'd pass both forward and inverse splines
    
    // Initialize output file if filename provided
    std::ofstream outfile;
    bool write_output = !output_filename.empty();
    if (write_output) {
        outfile.open(output_filename);
        if (!outfile.is_open()) {
            std::cerr << "Error opening file: " << output_filename << std::endl;
            write_output = false;
        } else {
            std::cout << "Writing to file: " << output_filename << std::endl;
            outfile << "log_r[cm],log_m[g],log_P[dyne/cm^2]\n";
        }
    }
    
    // Initial conditions
    double fraction = 4.0 / 3.0;
    double log_r_start = log10(r_start);
    double log_r_end = log10(r_end);
    double log_m0 = log10(fraction) + log10(M_PI) + 3.0*log10(r_start) + log10(rho_c);
    
    // For initial pressure, we need to find P(ρ_c)
    // Since we only have the inverse spline ρ(P), we need to search for the correct P
    // This is a limitation - in practice you'd want both forward and inverse splines
    // For now, estimate using a reasonable range and binary search
    
    // Simple search to find pressure that gives approximately the right density
    auto rho_of_logP = [&](double lp){ return gsl_spline_eval(spline_inv, lp, acc_inv); };
    
    double a = min_logP;
    double b = max_logP;
    double target = std::log10(rho_c);
    for (int it = 0; it < 60; ++it) {
        double m = 0.5 * (a + b);
        if (rho_of_logP(m) < target) a = m; else b = m;
    }
    double log_p0 = 0.5 * (a + b);
    std::vector<double> state = {log_m0, log_p0};
    std::vector<double> state_prev = state;

    std::cout << "Initial conditions (spline-based):" << std::endl;
    std::cout << "  log_r_start = " << log_r_start << std::endl;
    std::cout << "  log_m0 = " << log_m0 << std::endl;
    std::cout << "  log_p0 = " << log_p0 << std::endl;
    std::cout << "  rho_c = " << rho_c << " g/cm³" << std::endl;
    std::cout << "  EOS range: [" << min_logP << ", " << max_logP << "]" << std::endl;

    int idx = 0;
    double log_r = log_r_start;
    double log_r_prev = log_r_start;
    double current_dlogr = base_dlogr;
    double dlogr_min = 0.2 * base_dlogr;
    double dlogr_max = 5.0 * base_dlogr;

    double logP_stop = std::log10(pressure_threshold);

    while (log_r < log_r_end) {
        // Integration step using spline-based derivatives
        state = rk4_step(log_r, current_dlogr, state, 
        [spline_inv, acc_inv, min_logP, max_logP](double r, const std::vector<double>& s) {
            return tolman_oppenheimer_volkoff_derivatives_spline(r, s, spline_inv, acc_inv, min_logP, max_logP);
        });
        
        log_r += current_dlogr;

        double g_prev = state_prev[1] - logP_stop;
        double g_curr = state[1] - logP_stop;
        if (g_prev > 0.0 && g_curr <= 0.0) {
            // Linear interpolation to the crossing
            double t = g_prev / (g_prev - g_curr);
            double log_r_surf = log_r_prev + t * (log_r - log_r_prev);
            
            double log_m_surf = state_prev[0] + t * (state[0] - state_prev[0]);

            if (write_output) {
                outfile << log_r_surf << "," << log_m_surf << "," << logP_stop << "\n";
                outfile.close();
            }
            return {idx+1, log_m_surf, log_r_surf, true};
        }
        
        // Output to file if requested
        if (write_output) {
            outfile << log_r << "," << state[0] << "," << state[1] << "\n";
        }
        
        // EOS validity check
        if (state[1] < min_logP || state[1] > max_logP) {
            std::cout << "Pressure outside EOS range at log_r: " << log_r 
                      << " (log_P = " << state[1] << ")" << std::endl;
            break;
        }

        // Numerical stability check
        if (!std::isfinite(state[0]) || !std::isfinite(state[1])) {
            std::cout << "Non-finite values at log_r: " << log_r << std::endl;
            break;
        }

        // Progress reporting
        if (idx % 1000 == 0) { 
            std::cout << "log_r: " << log_r << ", log_m: " << state[0] 
                      << ", log_P: " << state[1] << ", step: " << current_dlogr << std::endl;
        }

        // Adaptive step size control
        if (use_adaptive_stepping && idx > 0) {
            // Adjust step size based on pressure gradient
            double slope = std::abs((state[1] - state_prev[1]) / (log_r - log_r_prev));
            // larger slope -> smaller step; smaller slope -> larger step
            double scale = (slope > 1.0 ? 0.5 : (slope < 0.1 ? 2.0 : 1.0));
            // Clamp step size to reasonable bounds
            current_dlogr = std::clamp(base_dlogr * scale, dlogr_min, dlogr_max);
        }

        if (log_r + current_dlogr > log_r_end) current_dlogr = log_r_end - log_r;

        state_prev = state;
        log_r_prev = log_r;
        ++idx;
    }

    std::cerr << "[NoSurface] rho_c=" << rho_c
          << " last log_r=" << log_r
          << " last log_P=" << state[1]
          << " logP_stop=" << logP_stop << "\n";
    
    if (write_output) outfile.close();
    return {idx, state[0], log_r, false};
}