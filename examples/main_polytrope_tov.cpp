/**
 * @file main.cpp
 * @brief Configurable TOV Equation Solver for Compact Stars
 * @author Karan Amit Kinariwala
 * @date 2025-01-27
 * 
 * @details
 * This program solves the Tolman-Oppenheimer-Volkoff (TOV) equations for compact stars
 * using configurable polytropic equations of state. It supports all four degenerate gas types:
 * - Electron non-relativistic (white dwarf cores)
 * - Electron relativistic (massive white dwarfs)
 * - Neutron non-relativistic (low-density neutron star regions)
 * - Neutron relativistic (high-density neutron star cores)
 * 
 * The program can be easily configured by changing the EOS_TYPE constant below.
 * 
 * Output files:
 * - Individual CSV files: [eos_name]_rhoc_[density].csv
 * - Summary file: tov_mass_radius_summary_[eos_name].csv
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <chrono>
#include <algorithm>
#include <iterator>
#include "non_rotating_stellar_structure.hpp"
#include "polytropic_eos.hpp"

// Physical constants
const double M_SUN = 1.989e33;  // Solar mass in grams
const double KM_TO_CM = 1e5;    // Conversion factor from km to cm

// ============================================================================
// CONFIGURATION SECTION - MODIFY THESE PARAMETERS AS NEEDED
// ============================================================================

// Choose the EOS type (uncomment one):
const PolytropicGasType EOS_TYPE = PolytropicGasType::NEUTRON_RELATIVISTIC;
// const PolytropicGasType EOS_TYPE = PolytropicGasType::NEUTRON_NON_RELATIVISTIC;
// const PolytropicGasType EOS_TYPE = PolytropicGasType::ELECTRON_RELATIVISTIC;
// const PolytropicGasType EOS_TYPE = PolytropicGasType::ELECTRON_NON_RELATIVISTIC;

// EOS-specific configuration parameters
struct EOSConfig {
    PolytropicGasType type;
    double mu_e;                // Mean molecular weight per electron
    double rho_min;             // Minimum central density (g/cm³)
    double rho_max;             // Maximum central density (g/cm³)
    int num_densities;          // Number of density points
    double r_end;               // Maximum radius (cm)
    std::string description;    // Description for output
};

/**
 * @brief Get configuration parameters for different EOS types
 */
EOSConfig getEOSConfig(PolytropicGasType eos_type) {
    switch(eos_type) {
        case PolytropicGasType::ELECTRON_NON_RELATIVISTIC:
            return {
                eos_type,
                2.0,                    // mu_e
                1e5,                    // rho_min (g/cm³) - white dwarf range
                1e8,                   // rho_max (g/cm³)
                30,                     // num_densities
                5e9,                    // r_end (cm) - 10,000 km for white dwarfs
                "White Dwarf Cores (Non-Relativistic Electrons)"
            };
            
        case PolytropicGasType::ELECTRON_RELATIVISTIC:
            return {
                eos_type,
                2.0,                    // mu_e
                1e10,                   // rho_min (g/cm³) - massive white dwarf range
                1e14,                   // rho_max (g/cm³)
                30,                     // num_densities
                5e8,                    // r_end (cm) - 5,000 km for massive white dwarfs
                "Massive White Dwarfs (Relativistic Electrons)"
            };
            
        case PolytropicGasType::NEUTRON_NON_RELATIVISTIC:
            return {
                eos_type,
                2.0,                    // mu_e (not used for neutrons)
                1e12,                   // rho_min (g/cm³) - low-density neutron star
                1e15,                   // rho_max (g/cm³)
                35,                     // num_densities
                3e6,                    // r_end (cm) - 30 km for neutron stars
                "Neutron Stars (Non-Relativistic Neutrons)"
            };
            
        case PolytropicGasType::NEUTRON_RELATIVISTIC:
             return {
                 eos_type,
                 2.0,                    // mu_e (not used for neutrons)
                 1e11,                   // rho_min (g/cm³) - very low for light neutron stars
                 1e13,                   // rho_max (g/cm³) - moderate to avoid overly massive stars
                 40,                     // num_densities
                 1e8,                    // r_end (cm) - 1000 km maximum (allow integration to find surface)
                 "Neutron Stars (Relativistic Neutrons)"
             };
            
        default:
            throw std::invalid_argument("Unknown EOS type");
    }
}

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

/**
 * @brief Generate logarithmically spaced central densities
 * 
 * @param rho_min Minimum central density (g/cm³)
 * @param rho_max Maximum central density (g/cm³)
 * @param num_points Number of density points
 * @return Vector of central densities
 */
std::vector<double> generateCentralDensities(double rho_min, double rho_max, int num_points) {
    std::vector<double> densities;
    
    if (num_points <= 1) {
        densities.push_back(rho_min);
        return densities;
    }
    
    double log_min = std::log10(rho_min);
    double log_max = std::log10(rho_max);
    double step = (log_max - log_min) / (num_points - 1);
    
    for (int i = 0; i < num_points; ++i) {
        double log_rho = log_min + i * step;
        densities.push_back(std::pow(10.0, log_rho));
    }
    
    return densities;
}

/**
 * @brief Write mass-radius summary to CSV file
 * 
 * @param filename Output filename
 * @param densities Vector of central densities
 * @param masses Vector of stellar masses (in solar masses)
 * @param radii Vector of stellar radii (in km)
 */
void writeMassRadiusSummary(const std::string& filename, 
                           const std::vector<double>& densities,
                           const std::vector<double>& masses,
                           const std::vector<double>& radii) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing." << std::endl;
        return;
    }
    
    // Write header
    file << "central_density[g/cm^3],mass[M_sun],radius[km],log_central_density,log_mass[g],log_radius[cm]" << std::endl;
    
    // Write data
    for (size_t i = 0; i < densities.size(); ++i) {
        if (i < masses.size() && i < radii.size()) {
            file << std::scientific << std::setprecision(6)
                 << densities[i] << ","
                 << masses[i] << ","
                 << radii[i] << ","
                 << std::log10(densities[i]) << ","
                 << std::log10(masses[i] * M_SUN) << ","
                 << std::log10(radii[i] * KM_TO_CM) << std::endl;
        }
    }
    
    file.close();
    std::cout << "Mass-radius summary written to: " << filename << std::endl;
}

/**
 * @brief Find radius corresponding to a specific mass using interpolation
 * 
 * @param target_mass Target mass in solar masses
 * @param masses Vector of stellar masses (in solar masses)
 * @param radii Vector of stellar radii (in km)
 * @param densities Vector of central densities (for reporting)
 * @return Interpolated radius in km, or -1 if not found
 */
double findRadiusForMass(double target_mass, 
                        const std::vector<double>& masses,
                        const std::vector<double>& radii,
                        const std::vector<double>& densities) {
    if (masses.size() != radii.size() || masses.empty()) {
        return -1.0;
    }
    
    // Find the closest masses above and below the target
    double closest_lower_mass = -1.0, closest_upper_mass = -1.0;
    double closest_lower_radius = -1.0, closest_upper_radius = -1.0;
    double closest_lower_density = -1.0, closest_upper_density = -1.0;
    
    for (size_t i = 0; i < masses.size(); ++i) {
        if (masses[i] <= target_mass) {
            if (closest_lower_mass < 0 || masses[i] > closest_lower_mass) {
                closest_lower_mass = masses[i];
                closest_lower_radius = radii[i];
                closest_lower_density = densities[i];
            }
        }
        if (masses[i] >= target_mass) {
            if (closest_upper_mass < 0 || masses[i] < closest_upper_mass) {
                closest_upper_mass = masses[i];
                closest_upper_radius = radii[i];
                closest_upper_density = densities[i];
            }
        }
    }
    
    // Check if we found exact match
    for (size_t i = 0; i < masses.size(); ++i) {
        if (std::abs(masses[i] - target_mass) < 1e-6) {
            std::cout << "  Exact match found!" << std::endl;
            std::cout << "  Mass: " << std::fixed << std::setprecision(3) << masses[i] << " M☉" << std::endl;
            std::cout << "  Radius: " << std::setprecision(2) << radii[i] << " km" << std::endl;
            std::cout << "  Central density: " << std::scientific << std::setprecision(2) 
                      << densities[i] << " g/cm³" << std::endl;
            return radii[i];
        }
    }
    
    // Interpolate if we have bounds
    if (closest_lower_mass >= 0 && closest_upper_mass >= 0 && 
        closest_lower_mass != closest_upper_mass) {
        
        double fraction = (target_mass - closest_lower_mass) / (closest_upper_mass - closest_lower_mass);
        double interpolated_radius = closest_lower_radius + fraction * (closest_upper_radius - closest_lower_radius);
        double interpolated_density = closest_lower_density * std::pow(closest_upper_density / closest_lower_density, fraction);
        
        std::cout << "  Interpolated result:" << std::endl;
        std::cout << "  Target mass: " << std::fixed << std::setprecision(3) << target_mass << " M☉" << std::endl;
        std::cout << "  Interpolated radius: " << std::setprecision(2) << interpolated_radius << " km" << std::endl;
        std::cout << "  Interpolated central density: " << std::scientific << std::setprecision(2) 
                  << interpolated_density << " g/cm³" << std::endl;
        std::cout << "  Bounds used:" << std::endl;
        std::cout << "    Lower: " << std::fixed << std::setprecision(3) << closest_lower_mass 
                  << " M☉ → " << std::setprecision(2) << closest_lower_radius << " km" << std::endl;
        std::cout << "    Upper: " << std::setprecision(3) << closest_upper_mass 
                  << " M☉ → " << std::setprecision(2) << closest_upper_radius << " km" << std::endl;
        
        return interpolated_radius;
    }
    
    // Check if target is outside our range
    if (closest_lower_mass >= 0 && closest_upper_mass < 0) {
        std::cout << "  Target mass (" << std::fixed << std::setprecision(3) << target_mass 
                  << " M☉) is above our maximum calculated mass (" 
                  << std::setprecision(3) << closest_lower_mass << " M☉)" << std::endl;
    } else if (closest_lower_mass < 0 && closest_upper_mass >= 0) {
        std::cout << "  Target mass (" << std::fixed << std::setprecision(3) << target_mass 
                  << " M☉) is below our minimum calculated mass (" 
                  << std::setprecision(3) << closest_upper_mass << " M☉)" << std::endl;
    }
    
    return -1.0;
}

/**
 * @brief Display program header and configuration
 */
void displayHeader(const EOSConfig& config) {
    std::cout << "========================================" << std::endl;
    std::cout << "  Configurable TOV Equation Solver" << std::endl;
    std::cout << "  " << config.description << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::endl;
    
    // Display EOS information
    PolytropicEOS eos;
    auto eos_data = eos.getEOSParameters(config.type);
    
    std::cout << "EOS Configuration:" << std::endl;
    std::cout << "  Type: " << eos_data.name << std::endl;
    std::cout << "  Polytropic index (γ): " << std::fixed << std::setprecision(5) << eos_data.gamma << std::endl;
    std::cout << "  Proportionality constant (k): " << std::scientific << std::setprecision(6) << eos_data.k << std::endl;
    if (config.type == PolytropicGasType::ELECTRON_NON_RELATIVISTIC || 
        config.type == PolytropicGasType::ELECTRON_RELATIVISTIC) {
        std::cout << "  Mean molecular weight (μₑ): " << std::fixed << std::setprecision(1) << config.mu_e << std::endl;
    }
    std::cout << std::endl;
}

// ============================================================================
// MAIN FUNCTION
// ============================================================================

/**
 * @brief Main function - Configurable TOV equation solver
 */
int main() {
    // Get configuration for selected EOS type
    EOSConfig config = getEOSConfig(EOS_TYPE);
    
    // Display program information
    displayHeader(config);
    
    // Integration parameters
    const double r_start = 10.0;      // Starting radius (cm)
    const double dlogr = 0.0001;        // Step size in log(r)
    
    std::cout << "Integration Parameters:" << std::endl;
    std::cout << "  Starting radius: " << r_start << " cm" << std::endl;
    std::cout << "  Maximum radius: " << config.r_end << " cm (" << config.r_end/KM_TO_CM << " km)" << std::endl;
    std::cout << "  Step size (dlog r): " << dlogr << std::endl;
    std::cout << std::endl;
    
    std::cout << "Central Density Range:" << std::endl;
    std::cout << "  Minimum: " << std::scientific << config.rho_min << " g/cm³" << std::endl;
    std::cout << "  Maximum: " << std::scientific << config.rho_max << " g/cm³" << std::endl;
    std::cout << "  Number of points: " << config.num_densities << std::endl;
    std::cout << std::endl;
    
    // Generate central densities
    std::vector<double> central_densities = generateCentralDensities(config.rho_min, config.rho_max, config.num_densities);
    
    // Storage for results
    std::vector<double> stellar_masses;  // in solar masses
    std::vector<double> stellar_radii;   // in km
    std::vector<double> successful_densities;
    
    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "Starting TOV integration..." << std::endl;
    std::cout << "Progress: [Density] -> [Steps] [Mass] [Radius]" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    // Solve TOV equations for each central density
    for (size_t i = 0; i < central_densities.size(); ++i) {
        double rho_c = central_densities[i];
        
        std::cout << std::fixed << std::setprecision(0) 
                  << "(" << (i+1) << "/" << central_densities.size() << ") "
                  << std::scientific << std::setprecision(2) << rho_c << " g/cm³ -> ";
        
        try {
            // Solve TOV equations
            auto result = non_rotating_stellar_structure(config.type, rho_c, r_start, config.r_end, dlogr, config.mu_e);
            
            int steps = std::get<0>(result);
            double log_mass = std::get<1>(result);
            double log_radius = std::get<2>(result);
            
            // Convert to physical units
            double mass_grams = std::pow(10.0, log_mass);
            double mass_solar = mass_grams / M_SUN;
            
            // Get actual radius from integration result
            double radius_cm = std::pow(10.0, log_radius);
            double radius_km = radius_cm / KM_TO_CM;
            
            // Store results
            successful_densities.push_back(rho_c);
            stellar_masses.push_back(mass_solar);
            stellar_radii.push_back(radius_km);
            
            std::cout << std::fixed << std::setprecision(0) << steps << " steps, "
                      << std::setprecision(3) << mass_solar << " M☉, "
                      << std::setprecision(1) << radius_km << " km" << std::endl;
                      
        } catch (const std::exception& e) {
            std::cout << "FAILED (" << e.what() << ")" << std::endl;
        } catch (...) {
            std::cout << "FAILED (unknown error)" << std::endl;
        }
    }
    
    // End timing
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Integration completed in " << duration.count() << " ms" << std::endl;
    std::cout << "Successful calculations: " << successful_densities.size() 
              << "/" << central_densities.size() << std::endl;
    std::cout << std::endl;
    
    // Write summary file
    if (!successful_densities.empty()) {
        PolytropicEOS eos;
        auto eos_data = eos.getEOSParameters(config.type);
        std::string summary_file = "data/tov_mass_radius_summary_" + eos_data.name + ".csv";
        writeMassRadiusSummary(summary_file, successful_densities, stellar_masses, stellar_radii);
        
        // Display summary statistics
        std::cout << "Results Summary:" << std::endl;
        
        auto min_mass = *std::min_element(stellar_masses.begin(), stellar_masses.end());
        auto max_mass = *std::max_element(stellar_masses.begin(), stellar_masses.end());
        auto min_radius = *std::min_element(stellar_radii.begin(), stellar_radii.end());
        auto max_radius = *std::max_element(stellar_radii.begin(), stellar_radii.end());
        
        std::cout << "  Mass range: " << std::fixed << std::setprecision(3) 
                  << min_mass << " - " << max_mass << " M☉" << std::endl;
        std::cout << "  Radius range: " << std::setprecision(1) 
                  << min_radius << " - " << max_radius << " km" << std::endl;
        
        // Find maximum mass (Chandrasekhar-like limit for this EOS)
        auto max_mass_it = std::max_element(stellar_masses.begin(), stellar_masses.end());
        if (max_mass_it != stellar_masses.end()) {
            size_t max_idx = std::distance(stellar_masses.begin(), max_mass_it);
            std::cout << "  Maximum mass: " << std::setprecision(3) << *max_mass_it << " M☉ "
                      << "at ρc = " << std::scientific << std::setprecision(2) 
                      << successful_densities[max_idx] << " g/cm³" << std::endl;
        }
        
        std::cout << std::endl;
        std::cout << "Individual stellar structure files saved in data/ directory" << std::endl;
        std::cout << "Summary file: " << summary_file << std::endl;
        
        // Find radius for 0.77 solar masses (if appropriate for this EOS)
        if (config.type == PolytropicGasType::NEUTRON_RELATIVISTIC || 
            config.type == PolytropicGasType::NEUTRON_NON_RELATIVISTIC) {
            std::cout << std::endl;
            std::cout << "========================================" << std::endl;
            std::cout << "  Finding Radius for 0.77 Solar Masses" << std::endl;
            std::cout << "========================================" << std::endl;
            
            double target_mass = 0.77;  // Solar masses
            double radius_077 = findRadiusForMass(target_mass, stellar_masses, stellar_radii, successful_densities);
            
            if (radius_077 > 0) {
                std::cout << std::endl;
                std::cout << "✓ SUCCESS: Found radius for " << std::fixed << std::setprecision(2) 
                          << target_mass << " M☉ = " << std::setprecision(2) << radius_077 << " km" << std::endl;
            } else {
                std::cout << std::endl;
                std::cout << "✗ Could not determine radius for " << std::fixed << std::setprecision(2) 
                          << target_mass << " M☉" << std::endl;
                std::cout << "  Try adjusting the central density range." << std::endl;
            }
        }
        
    } else {
        std::cout << "No successful calculations to summarize." << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << "Program completed successfully!" << std::endl;
    
    return 0;
} 