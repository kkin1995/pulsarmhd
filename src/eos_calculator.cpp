#include "eos_calculator.hpp"
#include "magnetic_bps.hpp"
#include "non_magnetic_ideal_npe_gas.hpp"
#include <stdexcept>
#include <iostream>

// Magnetic BPS EOS Calculator
class MagneticBPSCalculator : public EOSCalculator {
public:
    bool calculateEOS(const EOSParameters& params) override {
        try {
            MagneticBPSEOS eos(params.atomic_mass_file, 
                              params.B_ratio_electron,
                              params.rel_tolerance,
                              params.abs_tolerance);
            
            auto results = eos.runSimulation(params.nB_min, 
                                           params.nB_max, 
                                           params.num_points);
            
            eos.writeEOSResults(params.output_file, results);
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error in Magnetic BPS EOS calculation: " << e.what() << std::endl;
            return false;
        }
    }
    
    std::string getType() const override { return "Magnetic BPS"; }
};

// Non-magnetic NPE Gas Calculator
class NonMagneticNPECalculator : public EOSCalculator {
public:
    bool calculateEOS(const EOSParameters& params) override {
        try {
            NonMagneticNPEGas gas(params.debug_mode);
            gas.calculateEOS(params.output_file,
                           params.nB_min,
                           params.nB_max,
                           params.num_points);
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error in Non-magnetic NPE Gas calculation: " << e.what() << std::endl;
            return false;
        }
    }
    
    std::string getType() const override { return "Non-magnetic NPE Gas"; }
};

// Factory implementation
std::unique_ptr<EOSCalculator> EOSCalculatorFactory::createCalculator(EOSType type) {
    switch (type) {
        case EOSType::MAGNETIC_BPS:
            return std::make_unique<MagneticBPSCalculator>();
        case EOSType::NON_MAGNETIC_NPE_GAS:
            return std::make_unique<NonMagneticNPECalculator>();
        case EOSType::CUSTOM_EOS:
            throw std::runtime_error("Custom EOS not implemented yet");
        default:
            throw std::runtime_error("Unknown EOS type");
    }
}

// Function to calculate EOS based on parameters
bool calculateEOS(const std::string& eos_type, const EOSParameters& params) {
    // Validate parameters
    if (params.nB_min >= params.nB_max) {
        std::cerr << "Error: nB_min must be less than nB_max." << std::endl;
        return false;
    }

    try {
        // Parse EOS type
        EOSType type;
        if (eos_type == "MAGNETIC_BPS") {
            type = EOSType::MAGNETIC_BPS;
        } else if (eos_type == "NON_MAGNETIC_NPE_GAS") {
            type = EOSType::NON_MAGNETIC_NPE_GAS;
        } else {
            std::cerr << "Unknown EOS type: " << eos_type << std::endl;
            return false;
        }

        // Create calculator
        auto calculator = EOSCalculatorFactory::createCalculator(type);

        // Calculate EOS
        std::cout << "Calculating " << calculator->getType() << " EOS..." << std::endl;
        std::cout << "Parameters:" << std::endl;
        std::cout << "  nB range: " << params.nB_min << " to " << params.nB_max << std::endl;
        std::cout << "  Number of points: " << params.num_points << std::endl;
        std::cout << "  Output file: " << params.output_file << std::endl;
        if (type == EOSType::MAGNETIC_BPS) {
            std::cout << "  B ratio: " << params.B_ratio_electron << std::endl;
        }
        if (type == EOSType::NON_MAGNETIC_NPE_GAS) {
            std::cout << "  Debug mode: " << (params.debug_mode ? "enabled" : "disabled") << std::endl;
        }

        if (calculator->calculateEOS(params)) {
            std::cout << "Successfully calculated " << calculator->getType() << " EOS" << std::endl;
            std::cout << "Results written to: " << params.output_file << std::endl;
            return true;
        } else {
            std::cerr << "Failed to calculate EOS" << std::endl;
            return false;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return false;
    }
}

// Example usage in main (can be removed or kept for testing)
#ifdef TEST_EOS_CALCULATOR
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <EOS_TYPE> [options]" << std::endl;
        return 1;
    }

    EOSParameters params;
    // Parse command line options
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--nB-min") == 0 && i + 1 < argc) {
            params.nB_min = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--nB-max") == 0 && i + 1 < argc) {
            params.nB_max = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--num-points") == 0 && i + 1 < argc) {
            params.num_points = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            params.output_file = argv[++i];
        } else if (strcmp(argv[i], "--B-ratio") == 0 && i + 1 < argc) {
            params.B_ratio_electron = std::stod(argv[++i]);
        } else if (strcmp(argv[i], "--debug") == 0) {
            params.debug_mode = true;
        }
    }

    return calculateEOS(argv[1], params) ? 0 : 1;
}
#endif 