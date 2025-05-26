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
        std::string rho_str;
        std::string p_str;

        // Parse the line using comma as delimiter
        if (std::getline(ss, rho_str, ',') && std::getline(ss, p_str)) {
            try {
                double rho = std::stod(rho_str);
                double p = std::stod(p_str);
                
                log_rho.push_back(rho);
                log_P.push_back(p);
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
    double start = 1e15;  // Minimum central density (g/cm³)
    double end = 1e18;    // Maximum central density (g/cm³)
    double log_start = std::log10(start);
    double log_end = std::log10(end);
    double step = (log_end - log_start) / (num_points - 1);
    
    for (int i = 0; i < num_points; ++i) {
        double rho = std::pow(10.0, log_start + i * step);
        central_densities.push_back(rho);
    }
    
    // Integration parameters
    double r_start = 1.0;   // Starting radius (cm)
    double r_end = 1.0e6;   // Maximum radius (cm)
    double base_dlogr = 0.0001;  // Base step size in log(r)
    
    // Load EOS data
    std::string eos_filename = "/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/unified_eos_magnetic_BPS-BBP-Polytrope_B_001.csv";
    
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
        for (const auto& [rho_val, p_val] : data) {
            log_rho.push_back(rho_val);
            log_P.push_back(p_val);
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
        
        std::cout << "After cleaning: " << log_P.size() << " unique data points." << std::endl;
    } else {
        std::cerr << "Failed to load EOS data. Exiting." << std::endl;
        return 1;
    }

    // Determine EOS validity range
    double min_logP = log_P.front();
    double max_logP = log_P.back();
    
    std::cout << "EOS pressure range: [" << min_logP << ", " << max_logP << "] (log10 dyne/cm²)" << std::endl;

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
        auto [steps, log_mass_final, log_radius_final] = non_rotating_stellar_structure_spline(
            spline_inv,           // Inverse EOS spline ρ(P)
            acc_inv,              // GSL accelerator
            min_logP,             // EOS pressure range minimum
            max_logP,             // EOS pressure range maximum
            rho_c,                // Central density
            r_start,              // Starting radius
            r_end,                // Maximum radius
            base_dlogr,           // Base step size
            true,                 // Enable adaptive stepping
            1e-8,                 // Surface pressure threshold (dyne/cm²)
            filename              // Output file
        );
        
        // Convert results to physical units
        double mass_grams = std::pow(10.0, log_mass_final);
        double radius_cm = std::pow(10.0, log_radius_final);
        double mass_solar = mass_grams / 1.989e33;  // Convert to solar masses
        double radius_km = radius_cm / 1e5;         // Convert to kilometers
        
        // Report results
        std::cout << "Integration completed in " << steps << " steps" << std::endl;
        std::cout << "Final Mass = " << mass_grams << " g = " << mass_solar << " M☉" << std::endl;
        std::cout << "Final Radius = " << radius_cm << " cm = " << radius_km << " km" << std::endl;
        std::cout << "Compactness = " << (2.95 * mass_solar / radius_km) << std::endl;
    }

    // Clean up GSL resources
    gsl_spline_free(spline_inv);
    gsl_interp_accel_free(acc_inv);
    
    std::cout << "\n=== COMPUTATION COMPLETE ===" << std::endl;
    std::cout << "All stellar models computed successfully!" << std::endl;
    std::cout << "Output files written to data/ directory." << std::endl;

    return 0;
}