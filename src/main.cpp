#include "rk4.hpp"
#include "non_rotating_stellar_structure.hpp"
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <tuple>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h> 

bool readEoSData(const std::string& filename, std::vector<double>& log_rho, std::vector<double>& log_P) {
    // Clear any existing data
    log_rho.clear();
    log_P.clear();

    // Open the file
    std::ifstream file(filename);
    if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return false;
    }

    // Skip the header line
    std::string line;
    std::getline(file, line);

    // Read data line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string log_rho_str;
        std::string log_p_str;

        // Parse the line using comma as delimiter
        if (std::getline(ss, log_rho_str, ',') && std::getline(ss, log_p_str)) {
            try {
                double lgrho = std::stod(log_rho_str);
                double lgp = std::stod(log_p_str);
                
                log_rho.push_back(lgrho);
                log_P.push_back(lgp);

            } catch (const std::invalid_argument&) {
                std::cerr << "Error parsing line (invalid format): " << line << std::endl;
                continue;
            } catch (const std::out_of_range&) {
                std::cerr << "Error parsing line (value out of range): " << line << std::endl;
                continue;
            }
        }
    }

    file.close();

    // Check if we read any data
    if (log_rho.empty()) {
        std::cerr << "No data was read from the file." << std::endl;
        return false;
    }

    return true;
}

int main() {
    /*
     * MODULAR TOV STELLAR STRUCTURE SOLVER
     * 
     * This program has been refactored to use the modular TOV solver from
     * non_rotating_stellar_structure.cpp instead of embedded equations.
     * 
     * BENEFITS:
     * - Cleaner, more maintainable code
     * - Reuses tested TOV equation implementations
     * - Better separation of concerns
     * - Enhanced error handling and adaptive stepping
     * - Consistent physics across the codebase
     * 
     * FUNCTIONALITY:
     * - Loads realistic EOS data from CSV files
     * - Sets up GSL splines for EOS interpolation
     * - Scans over central densities to generate mass-radius relation
     * - Uses modular spline-based TOV solver for each stellar model
     * - Outputs individual stellar structure files for analysis
     */
    
    std::vector<double> log_rho;
    std::vector<double> log_P;
    std::vector<double> central_densities;
    
    // Configure central density scan
    int num_points = 25;  // Number of stellar models to compute
    double start = 1e14;  // Minimum central density (g/cm³)
    double end = 1e16;    // Maximum central density (g/cm³)
    double log_start = std::log10(start);
    double log_end = std::log10(end);
    double step = (log_end - log_start) / (num_points - 1);
    
    for (int i = 0; i < num_points; ++i) {
        double rho = std::pow(10.0, log_start + i * step);
        central_densities.push_back(rho);
    }
    
    // Integration parameters
    double r_start = 1.0;   // Starting radius (cm)
    double r_end = 1.0e8;   // Maximum radius (cm)
    double base_dlogr = 0.0001;  // Base step size in log(r)
    
    // Load EOS data
    std::string eos_filename = "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/unified_eos_magnetic_BPS-BBP-Polytrope_B_001.csv";
    
    if (readEoSData(eos_filename, log_rho, log_P)) {
        std::cout << "Successfully read " << log_rho.size() << " EOS data points." << std::endl;

        // Sort the data by log_rho for proper spline interpolation
        std::vector<std::pair<double, double>> data;
        for (size_t i = 0; i < log_rho.size(); i++) {
            data.push_back({log_rho[i], log_P[i]});
        }
        std::sort(data.begin(), data.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
            return a.first < b.first;
        });
        log_rho.clear();
        log_P.clear();
        for (const auto& [log_rho_val, log_p_val] : data) {
            log_rho.push_back(log_rho_val);
            log_P.push_back(log_p_val);

        }

        // Remove duplicate entries in log_P to ensure strict monotonicity for inverse spline
        std::vector<double> unique_log_rho;
        std::vector<double> unique_log_P;
        if (!log_P.empty()) {
            unique_log_rho.push_back(log_rho[0]);
            unique_log_P.push_back(log_P[0]);
            for (size_t i = 1; i < log_P.size(); i++) {
                if (log_P[i] > unique_log_P.back()) {
                    unique_log_P.push_back(log_P[i]);
                    unique_log_rho.push_back(log_rho[i]);
                }
            }
        }
        log_P = unique_log_P;
        log_rho = unique_log_rho;

        if (log_P.size() < 2) {
            std::cerr << "EOS has fewer than 2 unique (log_P, log_rho) points after cleaning.\n";
            return 1;
        }

        
        std::cout << "After cleaning: " << log_P.size() << " unique data points." << std::endl;
    } else {
        std::cerr << "Failed to load EOS data. Exiting." << std::endl;
        return 1;
    }

    // Determine EOS validity range
    double min_logP = log_P.front();
    double max_logP = log_P.back();
    
    std::cout << "EOS pressure range: [" << min_logP << ", " << max_logP << "] (log10 dyne/cm²)" << std::endl;

    // --- choose a physical “surface” density (typical outer envelope; tweak as you like)
    const double rho_surface = 1e6;              // g/cm^3 (try 1e6–1e8)
    const double target_log_rho = std::log10(rho_surface);
    const double min_logrho = log_rho.front();
    const double max_logrho = log_rho.back();
    double goal_log_rho = target_log_rho;
    if (goal_log_rho < min_logrho) {
        std::cerr << "[SurfaceWarn] Requested log10(rho_surface)=" << target_log_rho
                << " below EOS min=" << min_logrho << ". Clamping to min.\n";
        goal_log_rho = min_logrho + 1e-9; // tiny epsilon to stay inside
    } else if (goal_log_rho > max_logrho) {
        std::cerr << "[SurfaceWarn] Requested log10(rho_surface)=" << target_log_rho
                << " above EOS max=" << max_logrho << ". Clamping to max.\n";
        goal_log_rho = max_logrho - 1e-9;
    }

    // Set up GSL splines for EOS interpolation
    gsl_interp_accel *acc_inv = gsl_interp_accel_alloc();
    gsl_spline *spline_inv = gsl_spline_alloc(gsl_interp_steffen, log_P.size());
    
    if (gsl_spline_init(spline_inv, log_P.data(), log_rho.data(), log_P.size()) == GSL_SUCCESS) {
        std::cout << "Inverse spline ρ(P) initialized successfully." << std::endl;
    } else {
        std::cerr << "Error initializing inverse spline. Exiting." << std::endl;
        gsl_spline_free(spline_inv);
        gsl_interp_accel_free(acc_inv);
        return 1;
    }

    // Find logP_stop such that gsl_spline_eval(logP_stop) ≈ target_log_rho
    auto rho_of_logP = [&](double lp){ return gsl_spline_eval(spline_inv, lp, acc_inv); };

    double a = min_logP, b = max_logP;
    for (int it = 0; it < 80; ++it) {
        double mid = 0.5*(a+b);
        if (rho_of_logP(mid) < goal_log_rho) a = mid; else b = mid;
    }
    double logP_stop = 0.5*(a+b);
    double P_stop = std::pow(10.0, logP_stop);

    // Optional sanity print:
    std::cout << "Surface target: log10(rho)=" << goal_log_rho
            << " -> log10(P_stop)=" << logP_stop
          << " (rho_at_stop=" << rho_of_logP(logP_stop) << ")\n";

    // Generate mass-radius relation by scanning central densities
    std::cout << "\n=== GENERATING MASS-RADIUS RELATION ===" << std::endl;
    std::cout << "Computing " << central_densities.size() << " stellar models..." << std::endl;
    
    for (size_t i = 0; i < central_densities.size(); ++i) {
        double rho_c = central_densities[i];
        
        std::cout << "\n--- Model " << (i+1) << "/" << central_densities.size() 
                  << " (ρc = " << std::scientific << std::setprecision(2) << rho_c << " g/cm³) ---" << std::endl;

        // Generate output filename
        std::ostringstream oss;
        oss << "data/tov_solution_magnetic_bps_bbp_polytrope_"
            << std::scientific << std::setprecision(2) << rho_c
            << ".csv";
        std::string filename = oss.str();
        
        // Use modular TOV solver with spline-based EOS
        TovResult res = non_rotating_stellar_structure_spline(
            spline_inv,           // Inverse EOS spline ρ(P)
            acc_inv,              // GSL accelerator
            min_logP,             // EOS pressure range minimum
            max_logP,             // EOS pressure range maximum
            rho_c,                // Central density
            r_start,              // Starting radius
            r_end,                // Maximum radius
            base_dlogr,           // Base step size
            true,                 // Enable adaptive stepping
            P_stop,                 // Surface pressure threshold (dyne/cm²)
            filename              // Output file
        );
        
        if (!res.found_surface) {
            std::cout << "No surface bracketed before r_end; skipping. \n";
            continue;
        }

        // Convert results to physical units
        double mass_grams = std::pow(10.0, res.log10_m_surface);
        double radius_cm = std::pow(10.0, res.log10_r_surface);
        double mass_solar = mass_grams / 1.989e33;  // Convert to solar masses
        double radius_km = radius_cm / 1e5;         // Convert to kilometers
        
        // Report results
        std::cout << "Integration completed in " << res.steps << " steps" << std::endl;
        std::cout << "Final Mass = " << mass_grams << " g = " << mass_solar << " M☉" << std::endl;
        std::cout << "Final Radius = " << radius_cm << " cm = " << radius_km << " km" << std::endl;
        std::cout << "Compactness = " << (2.95325 * mass_solar / radius_km) << std::endl;
    }

    // Clean up GSL resources
    gsl_spline_free(spline_inv);
    gsl_interp_accel_free(acc_inv);
    
    std::cout << "\n=== COMPUTATION COMPLETE ===" << std::endl;
    std::cout << "All stellar models computed successfully!" << std::endl;
    std::cout << "Output files written to data/ directory." << std::endl;

    return 0;
}