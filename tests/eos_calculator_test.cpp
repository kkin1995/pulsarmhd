#include "eos_calculator.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>

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
        } catch (const std::filesystem::filesystem_error& e) {
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
    params.debug_mode = true;

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
    
    // Test invalid type
    EXPECT_THROW(EOSCalculatorFactory::createCalculator(EOSType::CUSTOM_EOS), 
                 std::runtime_error);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 