#include "polytropic_eos.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <algorithm>

// Constructor
PolytropicEOS::PolytropicEOS() {
    // Default constructor - no initialization needed
}

// Get EOS parameters for specific gas type
PolytropicEOSData PolytropicEOS::getEOSParameters(PolytropicGasType type) const {
    // Mean molecular weight per electron (used for electron gases)
    double mu_e = 2.0;
    
    switch(type) {
        case PolytropicGasType::ELECTRON_NON_RELATIVISTIC:
            return PolytropicEOSData(
                "electron_non_relativistic",
                1.0036e13 / std::pow(mu_e, 5.0 / 3.0),  // k
                5.0 / 3.0,                               // gamma
                mu_e
            );
            
        case PolytropicGasType::ELECTRON_RELATIVISTIC:
            return PolytropicEOSData(
                "electron_relativistic",
                1.2435e15 / std::pow(mu_e, 4.0 / 3.0),  // k
                4.0 / 3.0,                               // gamma
                mu_e
            );
            
        case PolytropicGasType::NEUTRON_NON_RELATIVISTIC:
            return PolytropicEOSData(
                "neutron_non_relativistic",
                5.3802e9,                                // k
                5.0 / 3.0,                               // gamma
                1.0                                      // Not applicable for neutrons
            );
            
        case PolytropicGasType::NEUTRON_RELATIVISTIC:
            return PolytropicEOSData(
                "neutron_relativistic",
                1.2293e15,                               // k
                4.0 / 3.0,                               // gamma
                1.0                                      // Not applicable for neutrons
            );
            
        default:
            throw std::invalid_argument("Unknown polytropic gas type");
    }
}

// Calculate pressure from density
double PolytropicEOS::calculatePressure(double density, double k, double gamma) const {
    if (density <= 0.0) {
        throw std::invalid_argument("Density must be positive");
    }
    
    return k * std::pow(density, gamma);
}

// Calculate density from pressure
double PolytropicEOS::calculateDensity(double pressure, double k, double gamma) const {
    if (pressure <= 0.0) {
        throw std::invalid_argument("Pressure must be positive");
    }
    if (k <= 0.0) {
        throw std::invalid_argument("EOS constant k must be positive");
    }
    
    return std::pow(pressure / k, 1.0 / gamma);
}

// Generate EOS table
bool PolytropicEOS::generateEOSTable(const PolytropicEOSParameters& params) const {
    try {
        // Get EOS parameters for the specified gas type
        PolytropicEOSData eos_data = getEOSParameters(params.gas_type);
        
        std::cout << "Generating polytropic EOS table for: " << eos_data.name << std::endl;
        std::cout << "  k = " << std::scientific << std::setprecision(4) << eos_data.k << std::endl;
        std::cout << "  γ = " << std::fixed << std::setprecision(6) << eos_data.gamma << std::endl;
        std::cout << "  ρ range: " << std::scientific << params.rho_min 
                  << " to " << params.rho_max << " g/cm³" << std::endl;
        std::cout << "  Number of points: " << params.num_points << std::endl;
        
        // Create density grid
        std::vector<double> densities = createDensityGrid(
            params.rho_min, params.rho_max, params.num_points, params.use_log_spacing);
        
        // Calculate pressures
        std::vector<double> pressures;
        pressures.reserve(densities.size());
        
        for (const double& rho : densities) {
            double pressure = calculatePressure(rho, eos_data.k, eos_data.gamma);
            pressures.push_back(pressure);
        }
        
        // Write to file
        bool success = writeEOSToFile(params.output_file, densities, pressures, 
                                     eos_data, params.output_log_values);
        
        if (success) {
            std::cout << "Successfully wrote EOS table to: " << params.output_file << std::endl;
        }
        
        return success;
        
    } catch (const std::exception& e) {
        std::cerr << "Error generating EOS table: " << e.what() << std::endl;
        return false;
    }
}

// Get all available EOS types
std::vector<std::pair<PolytropicGasType, PolytropicEOSData>> PolytropicEOS::getAllEOSTypes() const {
    std::vector<std::pair<PolytropicGasType, PolytropicEOSData>> all_types;
    
    all_types.emplace_back(PolytropicGasType::ELECTRON_NON_RELATIVISTIC, 
                          getEOSParameters(PolytropicGasType::ELECTRON_NON_RELATIVISTIC));
    all_types.emplace_back(PolytropicGasType::ELECTRON_RELATIVISTIC, 
                          getEOSParameters(PolytropicGasType::ELECTRON_RELATIVISTIC));
    all_types.emplace_back(PolytropicGasType::NEUTRON_NON_RELATIVISTIC, 
                          getEOSParameters(PolytropicGasType::NEUTRON_NON_RELATIVISTIC));
    all_types.emplace_back(PolytropicGasType::NEUTRON_RELATIVISTIC, 
                          getEOSParameters(PolytropicGasType::NEUTRON_RELATIVISTIC));
    
    return all_types;
}

// Convert gas type to string
std::string PolytropicEOS::gasTypeToString(PolytropicGasType type) {
    switch(type) {
        case PolytropicGasType::ELECTRON_NON_RELATIVISTIC:
            return "ELECTRON_NON_RELATIVISTIC";
        case PolytropicGasType::ELECTRON_RELATIVISTIC:
            return "ELECTRON_RELATIVISTIC";
        case PolytropicGasType::NEUTRON_NON_RELATIVISTIC:
            return "NEUTRON_NON_RELATIVISTIC";
        case PolytropicGasType::NEUTRON_RELATIVISTIC:
            return "NEUTRON_RELATIVISTIC";
        default:
            throw std::invalid_argument("Unknown gas type");
    }
}

// Convert string to gas type
PolytropicGasType PolytropicEOS::stringToGasType(const std::string& str) {
    if (str == "ELECTRON_NON_RELATIVISTIC") {
        return PolytropicGasType::ELECTRON_NON_RELATIVISTIC;
    } else if (str == "ELECTRON_RELATIVISTIC") {
        return PolytropicGasType::ELECTRON_RELATIVISTIC;
    } else if (str == "NEUTRON_NON_RELATIVISTIC") {
        return PolytropicGasType::NEUTRON_NON_RELATIVISTIC;
    } else if (str == "NEUTRON_RELATIVISTIC") {
        return PolytropicGasType::NEUTRON_RELATIVISTIC;
    } else {
        throw std::invalid_argument("Unknown gas type string: " + str);
    }
}

// Create density grid
std::vector<double> PolytropicEOS::createDensityGrid(double rho_min, double rho_max, 
                                                    int num_points, bool use_log_spacing) const {
    if (rho_min >= rho_max) {
        throw std::invalid_argument("rho_min must be less than rho_max");
    }
    if (num_points < 2) {
        throw std::invalid_argument("num_points must be at least 2");
    }
    
    std::vector<double> densities;
    densities.reserve(num_points);
    
    if (use_log_spacing) {
        // Logarithmic spacing
        double log_rho_min = std::log10(rho_min);
        double log_rho_max = std::log10(rho_max);
        double log_step = (log_rho_max - log_rho_min) / (num_points - 1);
        
        for (int i = 0; i < num_points; ++i) {
            double log_rho = log_rho_min + i * log_step;
            densities.push_back(std::pow(10.0, log_rho));
        }
    } else {
        // Linear spacing
        double step = (rho_max - rho_min) / (num_points - 1);
        
        for (int i = 0; i < num_points; ++i) {
            densities.push_back(rho_min + i * step);
        }
    }
    
    return densities;
}

// Write EOS to file
bool PolytropicEOS::writeEOSToFile(const std::string& filename,
                                  const std::vector<double>& densities,
                                  const std::vector<double>& pressures,
                                  const PolytropicEOSData& eos_data,
                                  bool output_log_values) const {
    if (densities.size() != pressures.size()) {
        std::cerr << "Error: densities and pressures vectors must have same size" << std::endl;
        return false;
    }
    
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing" << std::endl;
        return false;
    }
    
    // Write header with metadata
    outfile << "# Polytropic EOS: " << eos_data.name << std::endl;
    outfile << "# k = " << std::scientific << std::setprecision(6) << eos_data.k << std::endl;
    outfile << "# gamma = " << std::fixed << std::setprecision(6) << eos_data.gamma << std::endl;
    if (eos_data.name.find("electron") != std::string::npos) {
        outfile << "# mu_e = " << eos_data.mu_e << std::endl;
    }
    outfile << "# EOS relation: P = k * rho^gamma" << std::endl;
    outfile << "#" << std::endl;
    
    // Write column headers
    if (output_log_values) {
        outfile << "log10_rho[g/cm^3],log10_P[dyne/cm^2],rho[g/cm^3],P[dyne/cm^2]" << std::endl;
    } else {
        outfile << "rho[g/cm^3],P[dyne/cm^2]" << std::endl;
    }
    
    // Write data
    outfile << std::scientific << std::setprecision(6);
    
    for (size_t i = 0; i < densities.size(); ++i) {
        if (output_log_values) {
            outfile << std::log10(densities[i]) << ","
                   << std::log10(pressures[i]) << ","
                   << densities[i] << ","
                   << pressures[i] << std::endl;
        } else {
            outfile << densities[i] << "," << pressures[i] << std::endl;
        }
    }
    
    outfile.close();
    return true;
}

// Standalone function to calculate EOS
bool calculatePolytropicEOS(const std::string& gas_type_str, 
                           const PolytropicEOSParameters& params) {
    try {
        PolytropicEOS calculator;
        
        // Create a copy of params with the specified gas type
        PolytropicEOSParameters calc_params = params;
        calc_params.gas_type = PolytropicEOS::stringToGasType(gas_type_str);
        
        std::cout << "Calculating polytropic EOS for: " << gas_type_str << std::endl;
        
        return calculator.generateEOSTable(calc_params);
        
    } catch (const std::exception& e) {
        std::cerr << "Error in calculatePolytropicEOS: " << e.what() << std::endl;
        return false;
    }
}

// Optional: Main function for standalone usage
#ifdef BUILD_POLYTROPIC_EOS_STANDALONE
int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            std::cout << "Usage: " << argv[0] << " <GAS_TYPE> [options]" << std::endl;
            std::cout << "Available gas types:" << std::endl;
            std::cout << "  ELECTRON_NON_RELATIVISTIC" << std::endl;
            std::cout << "  ELECTRON_RELATIVISTIC" << std::endl;
            std::cout << "  NEUTRON_NON_RELATIVISTIC" << std::endl;
            std::cout << "  NEUTRON_RELATIVISTIC" << std::endl;
            std::cout << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  --rho-min <value>     Minimum density (g/cm³)" << std::endl;
            std::cout << "  --rho-max <value>     Maximum density (g/cm³)" << std::endl;
            std::cout << "  --num-points <value>  Number of points" << std::endl;
            std::cout << "  --output <filename>   Output file" << std::endl;
            std::cout << "  --linear-spacing      Use linear instead of log spacing" << std::endl;
            std::cout << "  --no-log-output       Don't output logarithmic values" << std::endl;
            std::cout << "  --list-all            List all available EOS types" << std::endl;
            return 1;
        }
        
        std::string gas_type = argv[1];
        
        // Handle special commands
        if (gas_type == "--list-all") {
            PolytropicEOS calculator;
            auto all_types = calculator.getAllEOSTypes();
            
            std::cout << "Available polytropic EOS types:" << std::endl;
            std::cout << std::string(60, '-') << std::endl;
            
            for (const auto& [type, data] : all_types) {
                std::cout << std::left << std::setw(30) << PolytropicEOS::gasTypeToString(type);
                std::cout << "k = " << std::scientific << std::setprecision(3) << data.k;
                std::cout << ", γ = " << std::fixed << std::setprecision(4) << data.gamma << std::endl;
            }
            return 0;
        }
        
        // Parse command line arguments
        PolytropicEOSParameters params;
        
        for (int i = 2; i < argc; i++) {
            std::string arg = argv[i];
            
            if (arg == "--rho-min" && i + 1 < argc) {
                params.rho_min = std::stod(argv[++i]);
            } else if (arg == "--rho-max" && i + 1 < argc) {
                params.rho_max = std::stod(argv[++i]);
            } else if (arg == "--num-points" && i + 1 < argc) {
                params.num_points = std::stoi(argv[++i]);
            } else if (arg == "--output" && i + 1 < argc) {
                params.output_file = argv[++i];
            } else if (arg == "--linear-spacing") {
                params.use_log_spacing = false;
            } else if (arg == "--no-log-output") {
                params.output_log_values = false;
            } else {
                std::cerr << "Unknown option: " << arg << std::endl;
                return 1;
            }
        }
        
        // Calculate EOS
        bool success = calculatePolytropicEOS(gas_type, params);
        
        return success ? 0 : 1;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
#endif 