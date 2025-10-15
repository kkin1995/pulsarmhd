#include "eos_calculator.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>

// Test fixture for EOS calculator tests
class EOSCalculatorTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Use /tmp directory for test outputs
    test_output_dir = "/tmp/eos_test_outputs";
    std::filesystem::create_directories(test_output_dir);
  }

  void TearDown() override {
    // Cleanup after each test
    try {
      std::filesystem::remove_all(test_output_dir);
    } catch (const std::filesystem::filesystem_error &e) {
      std::cerr << "Warning: Could not remove test directory: " << e.what() << std::endl;
    }
  }

  std::string test_output_dir;
  EOSParameters getDefaultParams() {
    EOSParameters params;
    params.nB_min = 1e22;
    params.nB_max = 1e37;
    params.num_points = 100;
    return params;
  }
};

// Test Magnetic BPS EOS calculation
TEST_F(EOSCalculatorTest, MagneticBPSCalculation) {
  // Arrange
  EOSParameters params = getDefaultParams();
  params.output_file = test_output_dir + "/magnetic_bps_test.csv";
  params.B_ratio_electron = 0.01;

  // Act
  bool result = calculateEOS("MAGNETIC_BPS", params);

  // Assert
  EXPECT_TRUE(result);
  EXPECT_TRUE(std::filesystem::exists(params.output_file));

  // Verify file contents
  std::ifstream file(params.output_file);
  EXPECT_TRUE(file.is_open());
  std::string first_line;
  std::getline(file, first_line);
  EXPECT_FALSE(first_line.empty());
}

// Test Non-magnetic NPE Gas EOS calculation
TEST_F(EOSCalculatorTest, NonMagneticNPECalculation) {
  // Arrange
  EOSParameters params = getDefaultParams();
  params.output_file = test_output_dir + "/non_magnetic_npe_test.csv";
  params.debug_mode = false;

  // Act
  bool result = calculateEOS("NON_MAGNETIC_NPE_GAS", params);

  // Assert
  EXPECT_TRUE(result);
  EXPECT_TRUE(std::filesystem::exists(params.output_file));

  // Verify file contents
  std::ifstream file(params.output_file);
  EXPECT_TRUE(file.is_open());
  std::string first_line;
  std::getline(file, first_line);
  EXPECT_FALSE(first_line.empty());
}

// Test invalid EOS type
TEST_F(EOSCalculatorTest, InvalidEOSType) {
  // Arrange
  EOSParameters params = getDefaultParams();

  // Act
  bool result = calculateEOS("INVALID_TYPE", params);

  // Assert
  EXPECT_FALSE(result);
}

// Test parameter validation
TEST_F(EOSCalculatorTest, ParameterValidation) {
  // Arrange
  EOSParameters params = getDefaultParams();
  params.nB_min = params.nB_max; // Invalid range

  // Act
  bool result = calculateEOS("MAGNETIC_BPS", params);

  // Assert
  EXPECT_FALSE(result);
}

// Test factory creation
TEST_F(EOSCalculatorTest, FactoryCreation) {
  // Test each valid type
  EXPECT_NO_THROW(EOSCalculatorFactory::createCalculator(EOSType::MAGNETIC_BPS));
  EXPECT_NO_THROW(EOSCalculatorFactory::createCalculator(EOSType::NON_MAGNETIC_NPE_GAS));
  EXPECT_NO_THROW(EOSCalculatorFactory::createCalculator(EOSType::POLYTROPIC_ELECTRON_NON_REL));
  EXPECT_NO_THROW(EOSCalculatorFactory::createCalculator(EOSType::POLYTROPIC_ELECTRON_REL));
  EXPECT_NO_THROW(EOSCalculatorFactory::createCalculator(EOSType::POLYTROPIC_NEUTRON_NON_REL));
  EXPECT_NO_THROW(EOSCalculatorFactory::createCalculator(EOSType::POLYTROPIC_NEUTRON_REL));

  // Test invalid type
  EXPECT_THROW(EOSCalculatorFactory::createCalculator(EOSType::CUSTOM_EOS), std::runtime_error);
}

// Test polytropic EOS integration with eos_calculator
TEST_F(EOSCalculatorTest, PolytropicElectronNonRelativistic) {
  EOSParameters params;
  params.output_file = test_output_dir + "/polytropic_electron_nonrel.csv";
  params.rho_min = 1e6;
  params.rho_max = 1e10;
  params.num_points = 50;
  params.mu_e = 2.0;

  bool success = calculateEOS("POLYTROPIC_ELECTRON_NON_REL", params);
  ASSERT_TRUE(success) << "Polytropic electron non-relativistic EOS calculation should succeed";

  // Check that output file was created
  ASSERT_TRUE(std::filesystem::exists(params.output_file)) << "Output file should be created";

  // Check file has expected content
  std::ifstream file(params.output_file);
  std::string line;
  int line_count = 0;
  while (std::getline(file, line)) {
    line_count++;
  }
  // Should have header lines + data lines
  EXPECT_GT(line_count, 50) << "Output file should have header and data lines";
}

TEST_F(EOSCalculatorTest, PolytropicElectronRelativistic) {
  EOSParameters params;
  params.output_file = test_output_dir + "/polytropic_electron_rel.csv";
  params.rho_min = 1e8;
  params.rho_max = 1e12;
  params.num_points = 30;
  params.mu_e = 2.0;

  bool success = calculateEOS("POLYTROPIC_ELECTRON_REL", params);
  ASSERT_TRUE(success) << "Polytropic electron relativistic EOS calculation should succeed";

  // Check that output file was created
  ASSERT_TRUE(std::filesystem::exists(params.output_file)) << "Output file should be created";
}

TEST_F(EOSCalculatorTest, PolytropicNeutronNonRelativistic) {
  EOSParameters params;
  params.output_file = test_output_dir + "/polytropic_neutron_nonrel.csv";
  params.rho_min = 1e12;
  params.rho_max = 1e15;
  params.num_points = 25;
  params.mu_e = 2.0; // Not used for neutron gas but should not cause issues

  bool success = calculateEOS("POLYTROPIC_NEUTRON_NON_REL", params);
  ASSERT_TRUE(success) << "Polytropic neutron non-relativistic EOS calculation should succeed";

  // Check that output file was created
  ASSERT_TRUE(std::filesystem::exists(params.output_file)) << "Output file should be created";
}

TEST_F(EOSCalculatorTest, PolytropicNeutronRelativistic) {
  EOSParameters params;
  params.output_file = test_output_dir + "/polytropic_neutron_rel.csv";
  params.rho_min = 1e13;
  params.rho_max = 1e16;
  params.num_points = 20;
  params.mu_e = 2.0; // Not used for neutron gas but should not cause issues

  bool success = calculateEOS("POLYTROPIC_NEUTRON_REL", params);
  ASSERT_TRUE(success) << "Polytropic neutron relativistic EOS calculation should succeed";

  // Check that output file was created
  ASSERT_TRUE(std::filesystem::exists(params.output_file)) << "Output file should be created";
}

TEST_F(EOSCalculatorTest, PolytropicCustomMuE) {
  EOSParameters params;
  params.output_file = test_output_dir + "/polytropic_custom_mue.csv";
  params.rho_min = 1e6;
  params.rho_max = 1e10;
  params.num_points = 30;
  params.mu_e = 3.0; // Custom value

  bool success = calculateEOS("POLYTROPIC_ELECTRON_NON_REL", params);
  ASSERT_TRUE(success) << "Polytropic EOS with custom mu_e should succeed";

  // Check that output file was created
  ASSERT_TRUE(std::filesystem::exists(params.output_file)) << "Output file should be created";
}

TEST_F(EOSCalculatorTest, PolytropicInvalidParameters) {
  EOSParameters params;
  params.output_file = test_output_dir + "/invalid_test.csv";
  params.rho_min = 1e10;
  params.rho_max = 1e6; // Invalid: min > max
  params.num_points = 10;

  bool success = calculateEOS("POLYTROPIC_ELECTRON_NON_REL", params);
  ASSERT_FALSE(success) << "Polytropic EOS with invalid density range should fail";

  // Test with invalid mu_e
  params.rho_min = 1e6;
  params.rho_max = 1e10;
  params.mu_e = -1.0; // Invalid: negative mu_e

  success = calculateEOS("POLYTROPIC_ELECTRON_NON_REL", params);
  ASSERT_FALSE(success) << "Polytropic EOS with negative mu_e should fail";
}

TEST_F(EOSCalculatorTest, PolytropicUnknownType) {
  EOSParameters params;
  params.output_file = test_output_dir + "/unknown_test.csv";

  bool success = calculateEOS("POLYTROPIC_UNKNOWN_TYPE", params);
  ASSERT_FALSE(success) << "Unknown polytropic EOS type should fail";
}

// Main function is provided by gtest_main library
