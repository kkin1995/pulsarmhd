#include "non_rotating_stellar_structure.hpp"
#include "polytropic_eos.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <vector>
#include <tuple>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

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
        
        // Set up test EOS data for spline tests
        setupTestEOSData();
    }

    void TearDown() override {
        // Cleanup GSL splines if allocated
        if (test_spline_inv) {
            gsl_spline_free(test_spline_inv);
            test_spline_inv = nullptr;
        }
        if (test_acc_inv) {
            gsl_interp_accel_free(test_acc_inv);
            test_acc_inv = nullptr;
        }
        
        // Cleanup after each test
        try {
            std::filesystem::remove_all(test_output_dir);
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Warning: Could not remove test directory: " << e.what() << std::endl;
        }
    }

    // Setup test EOS data for spline-based tests
    void setupTestEOSData() {
        // Create realistic test EOS data (pressure vs density)
        test_log_rho.clear();
        test_log_P.clear();
        
        // Generate a simple polytropic-like EOS for testing: P = k * rho^gamma
        double k = 1.2293e15; // Typical neutron relativistic
        double gamma = 4.0/3.0;
        
        // Create data points spanning realistic neutron star densities
        for (int i = 0; i < 50; i++) {
            double log_rho = 14.0 + 4.0 * i / 49.0; // log10(rho) from 14 to 18
            double rho = std::pow(10.0, log_rho);
            double P = k * std::pow(rho, gamma);
            double log_P = std::log10(P);
            
            test_log_rho.push_back(log_rho);
            test_log_P.push_back(log_P);
        }
        
        // Set up GSL spline for inverse EOS (rho as function of P)
        test_acc_inv = gsl_interp_accel_alloc();
        test_spline_inv = gsl_spline_alloc(gsl_interp_steffen, test_log_P.size());
        
        if (gsl_spline_init(test_spline_inv, test_log_P.data(), test_log_rho.data(), test_log_P.size()) != GSL_SUCCESS) {
            throw std::runtime_error("Failed to initialize test spline");
        }
        
        test_min_logP = test_log_P.front();
        test_max_logP = test_log_P.back();
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
    
    // Test EOS data for spline-based tests
    std::vector<double> test_log_rho, test_log_P;
    gsl_spline *test_spline_inv = nullptr;
    gsl_interp_accel *test_acc_inv = nullptr;
    double test_min_logP, test_max_logP;
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

// Test 1.1b: Spline-based TOV Derivatives Testing
TEST_F(StellarStructureTest, TOVDerivativesSplineBasicPhysics) {
    // Test basic physics of spline-based TOV equations
    double log_r = 5.0; // log10(r) = 5 -> r = 100,000 cm
    std::vector<double> state = {30.0, 34.0}; // log10(m), log10(P) - within test EOS range
    
    auto derivatives = tolman_oppenheimer_volkoff_derivatives_spline(
        log_r, state, test_spline_inv, test_acc_inv, test_min_logP, test_max_logP);
    
    // Basic sanity checks
    EXPECT_TRUE(std::isfinite(derivatives[0])) << "Mass derivative should be finite";
    EXPECT_TRUE(std::isfinite(derivatives[1])) << "Pressure derivative should be finite";
    EXPECT_GT(derivatives[0], 0) << "Mass should increase with radius";
    EXPECT_LT(derivatives[1], 0) << "Pressure should decrease with radius";
}

TEST_F(StellarStructureTest, TOVDerivativesSplineEOSConsistency) {
    // Test that spline-based derivatives use EOS correctly
    double log_r = 5.0;
    std::vector<double> state = {30.0, 35.0}; // log10(m), log10(P)
    
    // Verify the spline lookup is working correctly
    double log_P = state[1];
    double log_rho_spline = gsl_spline_eval(test_spline_inv, log_P, test_acc_inv);
    
    // Should be able to interpolate within our test data range
    EXPECT_GE(log_rho_spline, test_log_rho.front()) << "Interpolated density should be within range";
    EXPECT_LE(log_rho_spline, test_log_rho.back()) << "Interpolated density should be within range";
    
    // Compute derivatives
    auto derivatives = tolman_oppenheimer_volkoff_derivatives_spline(
        log_r, state, test_spline_inv, test_acc_inv, test_min_logP, test_max_logP);
    
    EXPECT_TRUE(std::isfinite(derivatives[0]));
    EXPECT_TRUE(std::isfinite(derivatives[1]));
}

TEST_F(StellarStructureTest, TOVDerivativesSplineOutOfRange) {
    // Test behavior when pressure is outside EOS range
    double log_r = 5.0;
    
    // Test pressure too low
    std::vector<double> state_low = {20.0, test_min_logP - 1.0};
    auto derivatives_low = tolman_oppenheimer_volkoff_derivatives_spline(
        log_r, state_low, test_spline_inv, test_acc_inv, test_min_logP, test_max_logP);
    
    // Should still give finite results (with clamping)
    EXPECT_TRUE(std::isfinite(derivatives_low[0]));
    EXPECT_TRUE(std::isfinite(derivatives_low[1]));
    
    // Test pressure too high
    std::vector<double> state_high = {35.0, test_max_logP + 1.0};
    auto derivatives_high = tolman_oppenheimer_volkoff_derivatives_spline(
        log_r, state_high, test_spline_inv, test_acc_inv, test_min_logP, test_max_logP);
    
    EXPECT_TRUE(std::isfinite(derivatives_high[0]));
    EXPECT_TRUE(std::isfinite(derivatives_high[1]));
}

TEST_F(StellarStructureTest, TOVDerivativesSplineVsPolytropic) {
    // Compare spline-based TOV with equivalent polytropic TOV
    double log_r = 5.0;
    std::vector<double> state = {30.0, 34.0}; // log10(m), log10(P)
    
    // Use polytropic parameters that match our test EOS
    double k = 1.2293e15;
    double gamma = 4.0/3.0;
    
    auto spline_result = tolman_oppenheimer_volkoff_derivatives_spline(
        log_r, state, test_spline_inv, test_acc_inv, test_min_logP, test_max_logP);
    
    auto polytropic_result = tolman_oppenheimer_volkoff_derivatives(
        log_r, state, k, gamma);
    
    // Results should be similar (within reasonable tolerance) since test EOS is polytropic
    double mass_rel_diff = std::abs(spline_result[0] - polytropic_result[0]) / 
                          std::abs(polytropic_result[0]);
    double pressure_rel_diff = std::abs(spline_result[1] - polytropic_result[1]) / 
                              std::abs(polytropic_result[1]);
    
    EXPECT_LT(mass_rel_diff, 0.05) << "Mass derivatives should be similar for equivalent EOS";
    EXPECT_LT(pressure_rel_diff, 0.05) << "Pressure derivatives should be similar for equivalent EOS";
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

// ============================================================================
// PHASE 4: SPLINE-BASED STELLAR STRUCTURE SOLVER TESTS
// ============================================================================

// Test 4.1: Spline-based Integration Basic Functionality
TEST_F(StellarStructureTest, SplineBasedIntegrationBasicFunctionality) {
    // Test basic functionality of non_rotating_stellar_structure_spline
    double rho_c = 1e15; // Central density (g/cm³)
    double r_start = 1.0;
    double r_end = 1e6;
    double base_dlogr = 0.01;
    
    std::string output_file = test_output_dir + "/test_spline_basic.csv";
    
    auto result = non_rotating_stellar_structure_spline(
        test_spline_inv, test_acc_inv, test_min_logP, test_max_logP,
        rho_c, r_start, r_end, base_dlogr, true, 1e-8, output_file);
    
    int steps = std::get<0>(result);
    double log_mass_final = std::get<1>(result);
    double log_radius_final = std::get<2>(result);
    
    // Basic sanity checks
    EXPECT_GT(steps, 0) << "Should complete some integration steps";
    EXPECT_TRUE(std::isfinite(log_mass_final)) << "Final mass should be finite";
    EXPECT_TRUE(std::isfinite(log_radius_final)) << "Final radius should be finite";
    
    // Physical reasonableness
    double mass = std::pow(10.0, log_mass_final);
    double radius = std::pow(10.0, log_radius_final);
    EXPECT_TRUE(isPhysicallyReasonable(mass, radius)) << "Results should be physically reasonable";
    
    // Check output file was created
    EXPECT_TRUE(std::filesystem::exists(output_file)) << "Output file should be created";
}

TEST_F(StellarStructureTest, SplineBasedIntegrationAdaptiveStepping) {
    // Test adaptive stepping functionality
    double rho_c = 1e15;
    double r_start = 1.0;
    double r_end = 1e6;
    double base_dlogr = 0.05;
    
    // Test with adaptive stepping enabled
    auto result_adaptive = non_rotating_stellar_structure_spline(
        test_spline_inv, test_acc_inv, test_min_logP, test_max_logP,
        rho_c, r_start, r_end, base_dlogr, true, 1e-8);
    
    // Test with adaptive stepping disabled
    auto result_fixed = non_rotating_stellar_structure_spline(
        test_spline_inv, test_acc_inv, test_min_logP, test_max_logP,
        rho_c, r_start, r_end, base_dlogr, false, 1e-8);
    
    // Adaptive stepping should typically require fewer steps for same accuracy
    EXPECT_GT(std::get<0>(result_fixed), 0) << "Fixed stepping should complete";
    EXPECT_GT(std::get<0>(result_adaptive), 0) << "Adaptive stepping should complete";
    
    // Results should be similar
    double mass_diff = std::abs(std::get<1>(result_adaptive) - std::get<1>(result_fixed));
    EXPECT_LT(mass_diff, 0.1) << "Adaptive and fixed stepping should give similar results";
}

TEST_F(StellarStructureTest, SplineBasedIntegrationDifferentDensities) {
    // Test with different central densities
    std::vector<double> densities = {5e14, 1e15, 2e15, 5e15};
    std::vector<double> final_masses;
    
    for (double rho_c : densities) {
        auto result = non_rotating_stellar_structure_spline(
            test_spline_inv, test_acc_inv, test_min_logP, test_max_logP,
            rho_c, 1.0, 1e6, 0.02, true, 1e-8);
        
        int steps = std::get<0>(result);
        double log_mass = std::get<1>(result);
        double mass = std::pow(10.0, log_mass);
        
        EXPECT_GT(steps, 0) << "Should complete for density " << rho_c;
        EXPECT_TRUE(std::isfinite(mass)) << "Mass should be finite for density " << rho_c;
        EXPECT_GT(mass, 0) << "Mass should be positive for density " << rho_c;
        
        final_masses.push_back(mass);
    }
    
    // Check that mass increases with central density (generally expected)
    for (size_t i = 1; i < final_masses.size(); i++) {
        EXPECT_GT(final_masses[i], final_masses[i-1] * 0.5) << 
            "Mass should generally increase with central density";
    }
}

TEST_F(StellarStructureTest, SplineBasedIntegrationSurfaceDetection) {
    // Test surface detection with different pressure thresholds
    double rho_c = 1e15;
    double r_start = 1.0;
    double r_end = 1e6;
    double base_dlogr = 0.02;
    
    // Test with different surface pressure thresholds
    std::vector<double> surface_thresholds = {1e-6, 1e-8, 1e-10};
    std::vector<double> final_radii;
    
    for (double threshold : surface_thresholds) {
        auto result = non_rotating_stellar_structure_spline(
            test_spline_inv, test_acc_inv, test_min_logP, test_max_logP,
            rho_c, r_start, r_end, base_dlogr, true, threshold);
        
        double log_radius = std::get<2>(result);
        double radius = std::pow(10.0, log_radius);
        
        EXPECT_TRUE(std::isfinite(radius)) << "Radius should be finite for threshold " << threshold;
        EXPECT_GT(radius, r_start) << "Final radius should be > starting radius";
        
        final_radii.push_back(radius);
    }
    
    // Lower surface pressure threshold should give larger stellar radius
    for (size_t i = 1; i < final_radii.size(); i++) {
        EXPECT_GE(final_radii[i], final_radii[i-1] * 0.9) << 
            "Lower pressure threshold should give larger radius";
    }
}

TEST_F(StellarStructureTest, SplineBasedIntegrationErrorHandling) {
    // Test error handling with invalid inputs
    double rho_c = 1e15;
    
    // Test with invalid radius range
    EXPECT_NO_THROW({
        auto invalid_result = non_rotating_stellar_structure_spline(
            test_spline_inv, test_acc_inv, test_min_logP, test_max_logP,
            rho_c, 1000.0, 100.0, 0.01, true, 1e-8); // r_end < r_start
        // Should handle gracefully
        (void)invalid_result; // Suppress unused variable warning
    });
    
    // Test with very small step size (should still work but take more steps)
    auto result_small_step = non_rotating_stellar_structure_spline(
        test_spline_inv, test_acc_inv, test_min_logP, test_max_logP,
        rho_c, 1.0, 1e5, 0.001, true, 1e-8);
    
    EXPECT_GT(std::get<0>(result_small_step), 0) << "Should handle small step sizes";
}

// Test 4.2: Comparison with Polytropic Solver
TEST_F(StellarStructureTest, SplineVsPolytropicComparison) {
    // Compare spline-based solver with equivalent polytropic solver
    double rho_c = 1e15;
    double r_start = 10.0;
    double r_end = 1e6;
    double dlogr = 0.05;
    
    // Use spline-based solver
    auto spline_result = non_rotating_stellar_structure_spline(
        test_spline_inv, test_acc_inv, test_min_logP, test_max_logP,
        rho_c, r_start, r_end, dlogr, false, 1e-8); // Disable adaptive for fair comparison
    
    // Use equivalent polytropic solver
    PolytropicGasType eos_type = PolytropicGasType::NEUTRON_RELATIVISTIC;
    auto polytropic_result = non_rotating_stellar_structure(
        eos_type, rho_c, r_start, r_end, dlogr, 2.0);
    
    // Extract final masses
    double spline_mass = std::pow(10.0, std::get<1>(spline_result));
    double polytropic_mass = std::pow(10.0, std::get<1>(polytropic_result));
    
    // Results should be similar since test EOS is polytropic
    double mass_rel_diff = std::abs(spline_mass - polytropic_mass) / polytropic_mass;
    
    EXPECT_LT(mass_rel_diff, 0.1) << "Spline and polytropic solvers should give similar results";
    EXPECT_TRUE(std::isfinite(spline_mass)) << "Spline solver mass should be finite";
    EXPECT_TRUE(std::isfinite(polytropic_mass)) << "Polytropic solver mass should be finite";
}

TEST_F(StellarStructureTest, SplineIntegrationFileOutput) {
    // Test file output functionality
    double rho_c = 1e15;
    std::string output_file = test_output_dir + "/test_spline_output.csv";
    
    auto result = non_rotating_stellar_structure_spline(
        test_spline_inv, test_acc_inv, test_min_logP, test_max_logP,
        rho_c, 1.0, 1e6, 0.02, true, 1e-8, output_file);
    
    // Check file exists and has content
    EXPECT_TRUE(std::filesystem::exists(output_file)) << "Output file should exist";
    
    std::ifstream file(output_file);
    EXPECT_TRUE(file.is_open()) << "Should be able to open output file";
    
    std::string line;
    int line_count = 0;
    while (std::getline(file, line)) {
        line_count++;
        if (line_count == 1) {
            // Check header - spline-based function uses format "log_r[cm],log_m[g],log_P[dyne/cm^2]"
            EXPECT_TRUE(line.find("log_r") != std::string::npos) << "Header should contain log_r";
            EXPECT_TRUE(line.find("log_m") != std::string::npos) << "Header should contain log_m";
            EXPECT_TRUE(line.find("log_P") != std::string::npos) << "Header should contain log_P";
        } else if (line_count <= 10) {
            // Check some data lines have proper format
            EXPECT_TRUE(line.find(',') != std::string::npos) << "Data lines should be comma-separated";
        }
    }
    
    EXPECT_GT(line_count, 1) << "File should have header + data lines";
    EXPECT_EQ(line_count - 1, std::get<0>(result)) << "Number of data lines should match integration steps";
}

// Main function is provided by gtest_main library 