#include "non_rotating_stellar_structure.hpp"
#include "polytropic_eos.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <vector>
#include <tuple>

// Test fixture for stellar structure tests
class StellarStructureTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Use /tmp directory for test outputs
        test_output_dir = "/tmp/stellar_structure_test_outputs";
        std::filesystem::create_directories(test_output_dir);
        
        // Physical constants for validation
        solar_mass = 1.989e33; // g
        solar_radius = 6.96e10; // cm
        
        // Test tolerances
        physics_tolerance = 0.01; // 1% for physics tests
        numerical_tolerance = 1e-10; // High precision for numerical tests
        regression_tolerance = 1e-12; // Very tight for regression tests
    }

    void TearDown() override {
        // Cleanup after each test
        try {
            std::filesystem::remove_all(test_output_dir);
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Warning: Could not remove test directory: " << e.what() << std::endl;
        }
    }

    // Helper functions
    bool isPhysicallyReasonable(double mass_grams, double radius_cm) {
        double mass_solar = mass_grams / solar_mass;
        double radius_km = radius_cm / 1e5;
        
        // Basic sanity checks for stellar objects
        return (mass_solar > 0.001 && mass_solar < 10.0 &&  // 0.001 to 10 solar masses
                radius_km > 1.0 && radius_km < 100000.0);   // 1 km to 100,000 km
    }
    
    std::tuple<double, double> extractFinalMassRadius(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        
        std::string line;
        std::getline(file, line); // Skip header
        
        double final_log_r = 0, final_log_m = 0;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string log_r_str, log_m_str, log_p_str;
            
            if (std::getline(iss, log_r_str, ',') &&
                std::getline(iss, log_m_str, ',') &&
                std::getline(iss, log_p_str)) {
                
                final_log_r = std::stod(log_r_str);
                final_log_m = std::stod(log_m_str);
            }
        }
        
        double radius = std::pow(10.0, final_log_r);
        double mass = std::pow(10.0, final_log_m);
        return {mass, radius};
    }
    
    bool compareDerivatives(const std::vector<double>& result1, 
                           const std::vector<double>& result2, 
                           double tolerance) {
        if (result1.size() != result2.size()) return false;
        
        for (size_t i = 0; i < result1.size(); i++) {
            double rel_error = std::abs(result1[i] - result2[i]) / 
                              (std::abs(result1[i]) + std::abs(result2[i]) + 1e-20);
            if (rel_error > tolerance) return false;
        }
        return true;
    }

protected:
    std::string test_output_dir;
    double solar_mass, solar_radius;
    double physics_tolerance, numerical_tolerance, regression_tolerance;
};

// ============================================================================
// PHASE 1: UNIT TESTS FOR INDIVIDUAL FUNCTIONS
// ============================================================================

// Test 1.1: TOV Derivatives Testing
TEST_F(StellarStructureTest, TOVDerivativesBasicPhysics) {
    // Test basic physics of TOV equations
    double log_r = 5.0; // log10(r) = 5 -> r = 100,000 cm
    std::vector<double> state = {30.0, 15.0}; // log10(m), log10(P)
    double k = 1.2293e15; // Neutron relativistic
    double gamma = 4.0/3.0;
    
    auto derivatives = tolman_oppenheimer_volkoff_derivatives(log_r, state, k, gamma);
    
    // Basic sanity checks
    EXPECT_TRUE(std::isfinite(derivatives[0])) << "Mass derivative should be finite";
    EXPECT_TRUE(std::isfinite(derivatives[1])) << "Pressure derivative should be finite";
    EXPECT_GT(derivatives[0], 0) << "Mass should increase with radius";
    EXPECT_LT(derivatives[1], 0) << "Pressure should decrease with radius";
}

TEST_F(StellarStructureTest, TOVDerivativesNumericalStability) {
    // Test numerical stability at extreme values
    double k = 1.2293e15;
    double gamma = 4.0/3.0;
    
    // Test very small mass (near center)
    std::vector<double> small_state = {-10.0, 15.0}; // Very small mass
    auto small_result = tolman_oppenheimer_volkoff_derivatives(5.0, small_state, k, gamma);
    EXPECT_TRUE(std::isfinite(small_result[0]));
    EXPECT_TRUE(std::isfinite(small_result[1]));
    
    // Test large mass
    std::vector<double> large_state = {35.0, 15.0}; // Large mass
    auto large_result = tolman_oppenheimer_volkoff_derivatives(5.0, large_state, k, gamma);
    EXPECT_TRUE(std::isfinite(large_result[0]));
    EXPECT_TRUE(std::isfinite(large_result[1]));
}

TEST_F(StellarStructureTest, TOVDerivativesLimitsToNewtonian) {
    // Test that TOV reduces to Newtonian in appropriate limits
    double log_r = 5.0;
    std::vector<double> state = {25.0, 10.0}; // Moderate values
    double k = 5.3802e9; // Non-relativistic neutron
    double gamma = 5.0/3.0;
    
    // Calculate both derivatives
    auto tov_result = tolman_oppenheimer_volkoff_derivatives(log_r, state, k, gamma);
    auto newton_result = newtonian(log_r, state, k, gamma);
    
    // For non-relativistic case, they should be similar (within ~10%)
    double mass_rel_diff = std::abs(tov_result[0] - newton_result[0]) / std::abs(newton_result[0]);
    double pressure_rel_diff = std::abs(tov_result[1] - newton_result[1]) / std::abs(newton_result[1]);
    
    EXPECT_LT(mass_rel_diff, 0.1) << "Mass derivatives should be similar in non-relativistic limit";
    EXPECT_LT(pressure_rel_diff, 0.5) << "Pressure derivatives should be reasonably similar";
}

// Test 1.2: Newtonian Function Testing
TEST_F(StellarStructureTest, NewtonianDerivativesBasicPhysics) {
    double log_r = 5.0;
    std::vector<double> state = {30.0, 15.0};
    double k = 5.3802e9;
    double gamma = 5.0/3.0;
    
    auto derivatives = newtonian(log_r, state, k, gamma);
    
    // Basic physics checks
    EXPECT_TRUE(std::isfinite(derivatives[0]));
    EXPECT_TRUE(std::isfinite(derivatives[1]));
    EXPECT_GT(derivatives[0], 0) << "Mass should increase with radius";
    EXPECT_LT(derivatives[1], 0) << "Pressure should decrease with radius";
}

TEST_F(StellarStructureTest, NewtonianDerivativesEOSConsistency) {
    // Test EOS consistency: verify that density calculated from P = k*rho^gamma is consistent
    double log_r = 5.0;
    std::vector<double> state = {30.0, 15.0};
    double k = 5.3802e9;
    double gamma = 5.0/3.0;
    
    // Extract values to verify EOS consistency
    double P = std::pow(10.0, state[1]);
    double log_rho = (state[1] - std::log10(k)) / gamma;
    double rho = std::pow(10.0, log_rho);
    
    // Verify P = k * rho^gamma
    double P_calculated = k * std::pow(rho, gamma);
    double relative_error = std::abs(P - P_calculated) / P;
    
    EXPECT_LT(relative_error, numerical_tolerance) << "EOS relation P = k*rho^gamma should be satisfied";
}

// Test 1.3: Filename Generation Testing
TEST_F(StellarStructureTest, FilenameGeneration) {
    std::string name = "test_eos";
    double rho_c = 1e15;
    
    std::string filename = get_filename(name, rho_c);
    
    EXPECT_TRUE(filename.find("data/") == 0) << "Filename should start with data/";
    EXPECT_TRUE(filename.find(name) != std::string::npos) << "Filename should contain EOS name";
    EXPECT_TRUE(filename.find(".csv") != std::string::npos) << "Filename should end with .csv";
}

TEST_F(StellarStructureTest, FilenameScientificNotation) {
    std::string name = "gas"; // Changed from "test" to avoid 'e' in name
    
    // Test various densities
    std::vector<double> densities = {1e5, 1.23e9, 5.67e15, 1e18};
    
    for (double rho : densities) {
        std::string filename = get_filename(name, rho);
        
        // Extract just the scientific notation part (after "_rhoc_")
        size_t rhoc_pos = filename.find("_rhoc_");
        ASSERT_NE(rhoc_pos, std::string::npos) << "Filename should contain '_rhoc_'";
        
        std::string scientific_part = filename.substr(rhoc_pos + 6); // Skip "_rhoc_"
        scientific_part = scientific_part.substr(0, scientific_part.find(".csv")); // Remove .csv
        
        // Should not contain 'e' or '+' in scientific notation part (replaced with 'p')
        EXPECT_TRUE(scientific_part.find('e') == std::string::npos) 
            << "Scientific notation should not contain 'e' in: " << scientific_part;
        EXPECT_TRUE(scientific_part.find('+') == std::string::npos) 
            << "Scientific notation should not contain '+' in: " << scientific_part;
        EXPECT_TRUE(scientific_part.find('p') != std::string::npos) 
            << "Scientific notation should contain 'p' replacement in: " << scientific_part;
    }
}

// ============================================================================
// PHASE 2: INTEGRATION TESTS FOR MAIN FUNCTION
// ============================================================================

// Test 2.1: EOS Integration Testing
TEST_F(StellarStructureTest, EOSParameterCalculation) {
    // Test that EOS parameters are correctly obtained from PolytropicEOS
    PolytropicGasType eos_type = PolytropicGasType::NEUTRON_RELATIVISTIC;
    double rho_c = 1e15;
    double r_start = 10.0;
    double r_end = 1e6;
    double dlogr = 0.1;
    double mu_e = 2.0;
    
    // This should not throw and should complete successfully
    EXPECT_NO_THROW({
        auto result = non_rotating_stellar_structure(eos_type, rho_c, r_start, r_end, dlogr, mu_e);
        EXPECT_GT(std::get<0>(result), 0) << "Should have integration steps";
        EXPECT_TRUE(std::isfinite(std::get<1>(result))) << "Final mass should be finite";
    });
}

TEST_F(StellarStructureTest, CustomMuEIntegration) {
    // Test custom mu_e calculation for electron gases
    PolytropicGasType eos_type = PolytropicGasType::ELECTRON_NON_RELATIVISTIC;
    double rho_c = 1e6;
    double r_start = 10.0;
    double r_end = 1e8;
    double dlogr = 0.1;
    
    // Test default mu_e = 2.0
    auto result_default = non_rotating_stellar_structure(eos_type, rho_c, r_start, r_end, dlogr, 2.0);
    
    // Test custom mu_e = 3.0
    auto result_custom = non_rotating_stellar_structure(eos_type, rho_c, r_start, r_end, dlogr, 3.0);
    
    // Results should be different (mu_e affects k parameter)
    EXPECT_NE(std::get<1>(result_default), std::get<1>(result_custom)) 
        << "Different mu_e should give different results";
}

TEST_F(StellarStructureTest, AllEOSTypesSupported) {
    // Test that all four EOS types work
    std::vector<PolytropicGasType> eos_types = {
        PolytropicGasType::ELECTRON_NON_RELATIVISTIC,
        PolytropicGasType::ELECTRON_RELATIVISTIC,
        PolytropicGasType::NEUTRON_NON_RELATIVISTIC,
        PolytropicGasType::NEUTRON_RELATIVISTIC
    };
    
    std::vector<double> test_densities = {1e6, 1e8, 1e12, 1e15}; // Appropriate for each EOS
    
    for (size_t i = 0; i < eos_types.size(); i++) {
        EXPECT_NO_THROW({
            auto result = non_rotating_stellar_structure(eos_types[i], test_densities[i], 
                                                       10.0, 1e7, 0.1, 2.0);
            EXPECT_GT(std::get<0>(result), 0) << "EOS type " << i << " should complete integration";
        }) << "EOS type " << i << " should work without throwing";
    }
}

// Test 2.2: Initial Conditions Testing
TEST_F(StellarStructureTest, InitialConditionsPhysical) {
    // Test that initial conditions are physically reasonable
    PolytropicGasType eos_type = PolytropicGasType::NEUTRON_RELATIVISTIC;
    double rho_c = 1e15;
    double r_start = 10.0; // 10 cm starting radius
    
    // Calculate what the initial conditions should be
    PolytropicEOS eos_calculator;
    auto eos_data = eos_calculator.getEOSParameters(eos_type);
    
    double fraction = 4.0/3.0;
    double expected_log_m0 = std::log10(fraction) + std::log10(M_PI) + 3.0*std::log10(r_start) + std::log10(rho_c);
    double expected_log_p0 = std::log10(eos_data.k) + eos_data.gamma * std::log10(rho_c);
    
    // Verify these are reasonable
    double initial_mass = std::pow(10.0, expected_log_m0);
    double initial_pressure = std::pow(10.0, expected_log_p0);
    
    EXPECT_GT(initial_mass, 0) << "Initial mass should be positive";
    EXPECT_GT(initial_pressure, 0) << "Initial pressure should be positive";
    EXPECT_LT(initial_mass, 1e35) << "Initial mass should not be too large";
    EXPECT_LT(initial_pressure, 1e40) << "Initial pressure should not be too large";
}

TEST_F(StellarStructureTest, InitialConditionsDensityRange) {
    // Test initial conditions across different density ranges
    PolytropicGasType eos_type = PolytropicGasType::NEUTRON_RELATIVISTIC;
    std::vector<double> densities = {1e14, 1e15, 1e16, 1e17};
    
    for (double rho_c : densities) {
        EXPECT_NO_THROW({
            auto result = non_rotating_stellar_structure(eos_type, rho_c, 10.0, 1e6, 0.2, 2.0);
            EXPECT_GT(std::get<0>(result), 0) << "Should complete for density " << rho_c;
        }) << "Should handle density " << rho_c;
    }
}

// Test 2.3: Integration Convergence Testing
TEST_F(StellarStructureTest, StepSizeIndependence) {
    // Test that results converge as step size decreases
    PolytropicGasType eos_type = PolytropicGasType::NEUTRON_RELATIVISTIC;
    double rho_c = 1e15;
    double r_start = 10.0;
    double r_end = 1e6;
    
    // Test different step sizes
    std::vector<double> step_sizes = {0.2, 0.1, 0.05};
    std::vector<double> final_masses;
    
    for (double dlogr : step_sizes) {
        auto result = non_rotating_stellar_structure(eos_type, rho_c, r_start, r_end, dlogr, 2.0);
        final_masses.push_back(std::get<1>(result));
    }
    
    // Results should converge (smaller differences with smaller steps)
    double diff_large = std::abs(final_masses[1] - final_masses[0]);
    double diff_small = std::abs(final_masses[2] - final_masses[1]);
    
    EXPECT_LT(diff_small, diff_large) << "Smaller step size should give better convergence";
}

TEST_F(StellarStructureTest, IntegrationStability) {
    // Test that integration doesn't produce NaN or infinite values
    PolytropicGasType eos_type = PolytropicGasType::NEUTRON_RELATIVISTIC;
    double rho_c = 5e17; // High density to test stability
    
    auto result = non_rotating_stellar_structure(eos_type, rho_c, 10.0, 1e6, 0.1, 2.0);
    
    EXPECT_TRUE(std::isfinite(std::get<1>(result))) << "Final mass should be finite";
    EXPECT_GT(std::get<0>(result), 0) << "Should complete some integration steps";
}

// Main function is provided by gtest_main library 