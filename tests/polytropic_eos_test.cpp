#include "polytropic_eos.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <set>

// Test fixture for Polytropic EOS calculator tests
class PolytropicEOSTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Use /tmp directory for test outputs
        test_output_dir = "/tmp/polytropic_eos_test_outputs";
        std::filesystem::create_directories(test_output_dir);
        calculator = std::make_unique<PolytropicEOS>();
    }

    void TearDown() override {
        // Cleanup after each test
        try {
            std::filesystem::remove_all(test_output_dir);
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Warning: Could not remove test directory: " << e.what() << std::endl;
        }
    }

    std::string test_output_dir;
    std::unique_ptr<PolytropicEOS> calculator;
    
    PolytropicEOSParameters getDefaultParams() {
        PolytropicEOSParameters params;
        params.rho_min = 1e10;
        params.rho_max = 1e15;
        params.num_points = 50;
        return params;
    }
};

// Test basic EOS parameter retrieval
TEST_F(PolytropicEOSTest, GetEOSParameters) {
    // Test electron non-relativistic
    auto electron_nr = calculator->getEOSParameters(PolytropicGasType::ELECTRON_NON_RELATIVISTIC);
    EXPECT_EQ(electron_nr.name, "electron_non_relativistic");
    EXPECT_NEAR(electron_nr.gamma, 5.0/3.0, 1e-10);
    EXPECT_GT(electron_nr.k, 0.0);
    EXPECT_EQ(electron_nr.mu_e, 2.0);

    // Test electron relativistic
    auto electron_r = calculator->getEOSParameters(PolytropicGasType::ELECTRON_RELATIVISTIC);
    EXPECT_EQ(electron_r.name, "electron_relativistic");
    EXPECT_NEAR(electron_r.gamma, 4.0/3.0, 1e-10);
    EXPECT_GT(electron_r.k, 0.0);

    // Test neutron non-relativistic
    auto neutron_nr = calculator->getEOSParameters(PolytropicGasType::NEUTRON_NON_RELATIVISTIC);
    EXPECT_EQ(neutron_nr.name, "neutron_non_relativistic");
    EXPECT_NEAR(neutron_nr.gamma, 5.0/3.0, 1e-10);
    EXPECT_DOUBLE_EQ(neutron_nr.k, 5.3802e9);

    // Test neutron relativistic
    auto neutron_r = calculator->getEOSParameters(PolytropicGasType::NEUTRON_RELATIVISTIC);
    EXPECT_EQ(neutron_r.name, "neutron_relativistic");
    EXPECT_NEAR(neutron_r.gamma, 4.0/3.0, 1e-10);
    EXPECT_DOUBLE_EQ(neutron_r.k, 1.2293e15);
}

// Test pressure calculation
TEST_F(PolytropicEOSTest, PressureCalculation) {
    double density = 1e15; // g/cm³
    double k = 1e10;
    double gamma = 5.0/3.0;
    
    double expected_pressure = k * std::pow(density, gamma);
    double calculated_pressure = calculator->calculatePressure(density, k, gamma);
    
    EXPECT_DOUBLE_EQ(calculated_pressure, expected_pressure);
    
    // Test invalid density
    EXPECT_THROW(calculator->calculatePressure(-1.0, k, gamma), std::invalid_argument);
    EXPECT_THROW(calculator->calculatePressure(0.0, k, gamma), std::invalid_argument);
}

// Test density calculation
TEST_F(PolytropicEOSTest, DensityCalculation) {
    double pressure = 1e20; // dyne/cm²
    double k = 1e10;
    double gamma = 5.0/3.0;
    
    double expected_density = std::pow(pressure / k, 1.0 / gamma);
    double calculated_density = calculator->calculateDensity(pressure, k, gamma);
    
    EXPECT_DOUBLE_EQ(calculated_density, expected_density);
    
    // Test round-trip calculation
    double test_density = 1e14;
    double test_pressure = calculator->calculatePressure(test_density, k, gamma);
    double back_to_density = calculator->calculateDensity(test_pressure, k, gamma);
    EXPECT_NEAR(back_to_density, test_density, 1e-10 * test_density);
    
    // Test invalid inputs
    EXPECT_THROW(calculator->calculateDensity(-1.0, k, gamma), std::invalid_argument);
    EXPECT_THROW(calculator->calculateDensity(pressure, -1.0, gamma), std::invalid_argument);
}

// Test EOS table generation
TEST_F(PolytropicEOSTest, EOSTableGeneration) {
    PolytropicEOSParameters params = getDefaultParams();
    params.gas_type = PolytropicGasType::NEUTRON_RELATIVISTIC;
    params.output_file = test_output_dir + "/neutron_rel_test.csv";
    params.use_log_spacing = true;
    params.output_log_values = true;
    
    bool result = calculator->generateEOSTable(params);
    
    EXPECT_TRUE(result);
    EXPECT_TRUE(std::filesystem::exists(params.output_file));
    
    // Verify file contents
    std::ifstream file(params.output_file);
    EXPECT_TRUE(file.is_open());
    
    std::string line;
    int line_count = 0;
    bool has_header = false;
    bool has_data = false;
    
    while (std::getline(file, line)) {
        line_count++;
        if (line.find("# Polytropic EOS:") != std::string::npos) {
            has_header = true;
        }
        if (line.find("log10_rho") != std::string::npos) {
            has_data = true;
        }
    }
    
    EXPECT_TRUE(has_header);
    EXPECT_TRUE(has_data);
    EXPECT_GT(line_count, params.num_points); // Should have header + data
}

// Test different gas types table generation
TEST_F(PolytropicEOSTest, AllGasTypesTableGeneration) {
    std::vector<PolytropicGasType> types = {
        PolytropicGasType::ELECTRON_NON_RELATIVISTIC,
        PolytropicGasType::ELECTRON_RELATIVISTIC,
        PolytropicGasType::NEUTRON_NON_RELATIVISTIC,
        PolytropicGasType::NEUTRON_RELATIVISTIC
    };
    
    for (const auto& type : types) {
        PolytropicEOSParameters params = getDefaultParams();
        params.gas_type = type;
        params.output_file = test_output_dir + "/test_" + 
                            PolytropicEOS::gasTypeToString(type) + ".csv";
        
        bool result = calculator->generateEOSTable(params);
        EXPECT_TRUE(result);
        EXPECT_TRUE(std::filesystem::exists(params.output_file));
    }
}

// Test string conversions
TEST_F(PolytropicEOSTest, StringConversions) {
    // Test gas type to string
    EXPECT_EQ(PolytropicEOS::gasTypeToString(PolytropicGasType::ELECTRON_NON_RELATIVISTIC), 
              "ELECTRON_NON_RELATIVISTIC");
    EXPECT_EQ(PolytropicEOS::gasTypeToString(PolytropicGasType::ELECTRON_RELATIVISTIC), 
              "ELECTRON_RELATIVISTIC");
    EXPECT_EQ(PolytropicEOS::gasTypeToString(PolytropicGasType::NEUTRON_NON_RELATIVISTIC), 
              "NEUTRON_NON_RELATIVISTIC");
    EXPECT_EQ(PolytropicEOS::gasTypeToString(PolytropicGasType::NEUTRON_RELATIVISTIC), 
              "NEUTRON_RELATIVISTIC");
    
    // Test string to gas type
    EXPECT_EQ(PolytropicEOS::stringToGasType("ELECTRON_NON_RELATIVISTIC"), 
              PolytropicGasType::ELECTRON_NON_RELATIVISTIC);
    EXPECT_EQ(PolytropicEOS::stringToGasType("ELECTRON_RELATIVISTIC"), 
              PolytropicGasType::ELECTRON_RELATIVISTIC);
    EXPECT_EQ(PolytropicEOS::stringToGasType("NEUTRON_NON_RELATIVISTIC"), 
              PolytropicGasType::NEUTRON_NON_RELATIVISTIC);
    EXPECT_EQ(PolytropicEOS::stringToGasType("NEUTRON_RELATIVISTIC"), 
              PolytropicGasType::NEUTRON_RELATIVISTIC);
    
    // Test invalid string
    EXPECT_THROW(PolytropicEOS::stringToGasType("INVALID_TYPE"), std::invalid_argument);
}

// Test get all EOS types
TEST_F(PolytropicEOSTest, GetAllEOSTypes) {
    auto all_types = calculator->getAllEOSTypes();
    
    EXPECT_EQ(all_types.size(), 4);
    
    // Check that all expected types are present
    std::set<PolytropicGasType> expected_types = {
        PolytropicGasType::ELECTRON_NON_RELATIVISTIC,
        PolytropicGasType::ELECTRON_RELATIVISTIC,
        PolytropicGasType::NEUTRON_NON_RELATIVISTIC,
        PolytropicGasType::NEUTRON_RELATIVISTIC
    };
    
    std::set<PolytropicGasType> actual_types;
    for (const auto& [type, data] : all_types) {
        actual_types.insert(type);
        EXPECT_FALSE(data.name.empty());
        EXPECT_GT(data.k, 0.0);
        EXPECT_GT(data.gamma, 0.0);
    }
    
    EXPECT_EQ(actual_types, expected_types);
}

// Test parameter validation
TEST_F(PolytropicEOSTest, ParameterValidation) {
    PolytropicEOSParameters params = getDefaultParams();
    
    // Test invalid density range
    params.rho_min = params.rho_max;
    params.output_file = test_output_dir + "/invalid_range.csv";
    EXPECT_FALSE(calculator->generateEOSTable(params));
    
    // Test invalid number of points
    params = getDefaultParams();
    params.num_points = 1;
    params.output_file = test_output_dir + "/invalid_points.csv";
    EXPECT_FALSE(calculator->generateEOSTable(params));
}

// Test standalone function
TEST_F(PolytropicEOSTest, StandaloneFunction) {
    PolytropicEOSParameters params = getDefaultParams();
    params.output_file = test_output_dir + "/standalone_test.csv";
    
    // Test valid gas type
    bool result = calculatePolytropicEOS("NEUTRON_RELATIVISTIC", params);
    EXPECT_TRUE(result);
    EXPECT_TRUE(std::filesystem::exists(params.output_file));
    
    // Test invalid gas type
    params.output_file = test_output_dir + "/standalone_invalid.csv";
    result = calculatePolytropicEOS("INVALID_TYPE", params);
    EXPECT_FALSE(result);
}

// Test linear vs logarithmic spacing
TEST_F(PolytropicEOSTest, SpacingModes) {
    PolytropicEOSParameters params = getDefaultParams();
    params.num_points = 10;
    
    // Test logarithmic spacing
    params.use_log_spacing = true;
    params.output_file = test_output_dir + "/log_spacing.csv";
    EXPECT_TRUE(calculator->generateEOSTable(params));
    
    // Test linear spacing
    params.use_log_spacing = false;
    params.output_file = test_output_dir + "/linear_spacing.csv";
    EXPECT_TRUE(calculator->generateEOSTable(params));
    
    // Both files should exist and have content
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/log_spacing.csv"));
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/linear_spacing.csv"));
}

// Test consistency with hardcoded values from original stellar structure code
TEST_F(PolytropicEOSTest, ConsistencyWithOriginalCode) {
    double mu_e = 2.0;
    
    // Test electron non-relativistic
    auto electron_nr = calculator->getEOSParameters(PolytropicGasType::ELECTRON_NON_RELATIVISTIC);
    double expected_k_enr = 1.0036e13 / std::pow(mu_e, 5.0/3.0);
    EXPECT_NEAR(electron_nr.k, expected_k_enr, 1e-10 * expected_k_enr);
    
    // Test electron relativistic  
    auto electron_r = calculator->getEOSParameters(PolytropicGasType::ELECTRON_RELATIVISTIC);
    double expected_k_er = 1.2435e15 / std::pow(mu_e, 4.0/3.0);
    EXPECT_NEAR(electron_r.k, expected_k_er, 1e-10 * expected_k_er);
    
    // Test neutron non-relativistic
    auto neutron_nr = calculator->getEOSParameters(PolytropicGasType::NEUTRON_NON_RELATIVISTIC);
    EXPECT_DOUBLE_EQ(neutron_nr.k, 5.3802e9);
    EXPECT_NEAR(neutron_nr.gamma, 5.0/3.0, 1e-15);
    
    // Test neutron relativistic
    auto neutron_r = calculator->getEOSParameters(PolytropicGasType::NEUTRON_RELATIVISTIC);
    EXPECT_DOUBLE_EQ(neutron_r.k, 1.2293e15);
    EXPECT_NEAR(neutron_r.gamma, 4.0/3.0, 1e-15);
}

// Main function is provided by gtest_main library 