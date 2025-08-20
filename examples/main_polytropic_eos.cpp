/**
 * @file main_polytropic_eos.cpp
 * @brief Example usage of the Polytropic EOS Calculator
 * 
 * This program demonstrates how to use the new polytropic EOS calculator
 * that extracts the previously hardcoded EOS calculations from the
 * stellar structure solver.
 * 
 * @author Karan Amit Kinariwala
 * @date 2025-01-09
 */

#include "polytropic_eos.hpp"
#include <iostream>
#include <iomanip>
#include <vector>

void demonstrateBasicUsage() {
    std::cout << "\n=== Basic Usage Demonstration ===" << std::endl;
    
    PolytropicEOS calculator;
    
    // Get EOS parameters for different gas types
    auto electron_rel = calculator.getEOSParameters(PolytropicGasType::ELECTRON_RELATIVISTIC);
    auto neutron_rel = calculator.getEOSParameters(PolytropicGasType::NEUTRON_RELATIVISTIC);
    
    std::cout << "\nElectron Relativistic EOS:" << std::endl;
    std::cout << "  Name: " << electron_rel.name << std::endl;
    std::cout << "  k = " << std::scientific << electron_rel.k << std::endl;
    std::cout << "  γ = " << std::fixed << std::setprecision(4) << electron_rel.gamma << std::endl;
    
    std::cout << "\nNeutron Relativistic EOS:" << std::endl;
    std::cout << "  Name: " << neutron_rel.name << std::endl;
    std::cout << "  k = " << std::scientific << neutron_rel.k << std::endl;
    std::cout << "  γ = " << std::fixed << std::setprecision(4) << neutron_rel.gamma << std::endl;
    
    // Calculate pressure for a given density
    double density = 1e15; // g/cm³
    double pressure_electron = calculator.calculatePressure(density, electron_rel.k, electron_rel.gamma);
    double pressure_neutron = calculator.calculatePressure(density, neutron_rel.k, neutron_rel.gamma);
    
    std::cout << "\nPressure calculations for ρ = " << std::scientific << density << " g/cm³:" << std::endl;
    std::cout << "  Electron relativistic: P = " << pressure_electron << " dyne/cm²" << std::endl;
    std::cout << "  Neutron relativistic:  P = " << pressure_neutron << " dyne/cm²" << std::endl;
}

void demonstrateEOSTableGeneration() {
    std::cout << "\n=== EOS Table Generation Demonstration ===" << std::endl;
    
    PolytropicEOS calculator;
    
    // Generate EOS table for neutron relativistic gas
    PolytropicEOSParameters params;
    params.gas_type = PolytropicGasType::NEUTRON_RELATIVISTIC;
    params.rho_min = 1e14;              // Minimum density
    params.rho_max = 1e18;              // Maximum density  
    params.num_points = 100;            // Number of points
    params.output_file = "data/neutron_relativistic_eos_example.csv";
    params.use_log_spacing = true;      // Logarithmic spacing
    params.output_log_values = true;    // Include log values in output
    
    std::cout << "\nGenerating EOS table with parameters:" << std::endl;
    std::cout << "  Gas type: " << PolytropicEOS::gasTypeToString(params.gas_type) << std::endl;
    std::cout << "  Density range: " << std::scientific << params.rho_min 
              << " to " << params.rho_max << " g/cm³" << std::endl;
    std::cout << "  Number of points: " << params.num_points << std::endl;
    std::cout << "  Output file: " << params.output_file << std::endl;
    
    bool success = calculator.generateEOSTable(params);
    
    if (success) {
        std::cout << "✓ EOS table generated successfully!" << std::endl;
    } else {
        std::cout << "✗ Failed to generate EOS table" << std::endl;
    }
}

void compareAllEOSTypes() {
    std::cout << "\n=== Comparison of All EOS Types ===" << std::endl;
    
    PolytropicEOS calculator;
    auto all_types = calculator.getAllEOSTypes();
    
    std::cout << "\nAvailable polytropic EOS types:" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::left << std::setw(30) << "Type" 
              << std::setw(15) << "k" 
              << std::setw(10) << "γ" 
              << std::setw(25) << "Description" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    for (const auto& [type, data] : all_types) {
        std::string description;
        switch(type) {
            case PolytropicGasType::ELECTRON_NON_RELATIVISTIC:
                description = "White dwarf cores";
                break;
            case PolytropicGasType::ELECTRON_RELATIVISTIC:
                description = "Massive white dwarfs";
                break;
            case PolytropicGasType::NEUTRON_NON_RELATIVISTIC:
                description = "Low-density NS regions";
                break;
            case PolytropicGasType::NEUTRON_RELATIVISTIC:
                description = "High-density NS cores";
                break;
        }
        
        std::cout << std::left << std::setw(30) << data.name;
        std::cout << std::scientific << std::setprecision(2) << std::setw(15) << data.k;
        std::cout << std::fixed << std::setprecision(4) << std::setw(10) << data.gamma;
        std::cout << std::setw(25) << description << std::endl;
    }
    
    // Compare pressures at a reference density
    double ref_density = 1e15; // g/cm³
    std::cout << "\nPressure comparison at ρ = " << std::scientific << ref_density << " g/cm³:" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    for (const auto& [type, data] : all_types) {
        double pressure = calculator.calculatePressure(ref_density, data.k, data.gamma);
        std::cout << std::left << std::setw(30) << data.name;
        std::cout << "P = " << std::scientific << std::setprecision(3) << pressure << " dyne/cm²" << std::endl;
    }
}

void generateComparisonTables() {
    std::cout << "\n=== Generating Comparison Tables ===" << std::endl;
    
    // Generate tables for all EOS types for comparison
    std::vector<PolytropicGasType> types = {
        PolytropicGasType::ELECTRON_NON_RELATIVISTIC,
        PolytropicGasType::ELECTRON_RELATIVISTIC,
        PolytropicGasType::NEUTRON_NON_RELATIVISTIC,
        PolytropicGasType::NEUTRON_RELATIVISTIC
    };
    
    for (const auto& gas_type : types) {
        PolytropicEOSParameters params;
        params.gas_type = gas_type;
        params.rho_min = 1e10;
        params.rho_max = 1e18;
        params.num_points = 200;
        
        // Create filename based on gas type
        PolytropicEOS calculator;
        auto eos_data = calculator.getEOSParameters(gas_type);
        params.output_file = "data/" + eos_data.name + "_comparison.csv";
        
        std::cout << "Generating: " << params.output_file << " ... ";
        
        bool success = calculator.generateEOSTable(params);
        if (success) {
            std::cout << "✓" << std::endl;
        } else {
            std::cout << "✗" << std::endl;
        }
    }
}

int main() {
    std::cout << "Polytropic EOS Calculator - Example Usage" << std::endl;
    std::cout << "===========================================" << std::endl;
    
    try {
        // Demonstrate basic usage
        demonstrateBasicUsage();
        
        // Demonstrate EOS table generation
        demonstrateEOSTableGeneration();
        
        // Compare all EOS types
        compareAllEOSTypes();
        
        // Generate comparison tables
        generateComparisonTables();
        
        std::cout << "\n=== Summary ===" << std::endl;
        std::cout << "The polytropic EOS calculator successfully extracted the hardcoded" << std::endl;
        std::cout << "EOS calculations from the stellar structure solver into a dedicated" << std::endl;
        std::cout << "and flexible module. This provides:" << std::endl;
        std::cout << "  • Clean separation of concerns" << std::endl;
        std::cout << "  • Reusable EOS calculations" << std::endl;
        std::cout << "  • Easy extension to new EOS types" << std::endl;
        std::cout << "  • Comprehensive documentation and examples" << std::endl;
        std::cout << "  • Both programmatic and command-line interfaces" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
} 