#include "eos_calculator.hpp"

#include "magnetic_bps.hpp"
#include "non_magnetic_ideal_npe_gas.hpp"
#include "polytropic_eos.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Magnetic BPS EOS Calculator
class MagneticBPSCalculator : public EOSCalculator {
public:
  bool calculateEOS(const EOSParameters &params) override {
    try {
      MagneticBPSEOS eos(params.atomic_mass_file, params.B_ratio_electron, params.rel_tolerance,
                         params.abs_tolerance);

      auto results = eos.runSimulation(params.nB_min, params.nB_max, params.num_points);

      eos.writeEOSResults(params.output_file, results);
      return true;
    } catch (const std::exception &e) {
      std::cerr << "Error in Magnetic BPS EOS calculation: " << e.what() << '\n';
      return false;
    }
  }

  [[nodiscard]] std::string getType() const override { return "Magnetic BPS"; }
};

// Non-magnetic NPE Gas Calculator
class NonMagneticNPECalculator : public EOSCalculator {
public:
  bool calculateEOS(const EOSParameters &params) override {
    try {
      NonMagneticNPEGas gas(params.debug_mode);
      gas.calculateEOS(params.output_file, params.nB_min, params.nB_max, params.num_points);
      return true;
    } catch (const std::exception &e) {
      std::cerr << "Error in Non-magnetic NPE Gas calculation: " << e.what() << '\n';
      return false;
    }
  }

  [[nodiscard]] std::string getType() const override { return "Non-magnetic NPE Gas"; }
};

// Polytropic EOS Calculator
class PolytropicEOSCalculator : public EOSCalculator {
private:
  PolytropicGasType gas_type;
  std::string type_name;

public:
  PolytropicEOSCalculator(PolytropicGasType type, const std::string &name)
      : gas_type(type), type_name(name) {}

  bool calculateEOS(const EOSParameters &params) override {
    try {
      PolytropicEOS calculator;

      // Get EOS parameters for the specific gas type
      auto eos_data = calculator.getEOSParameters(gas_type);

      // Validate density range
      if (params.rho_min >= params.rho_max) {
        std::cerr << "Error: rho_min must be less than rho_max for polytropic EOS." << '\n';
        return false;
      }

      // Open output file
      std::ofstream outfile(params.output_file);
      if (!outfile.is_open()) {
        std::cerr << "Error opening output file: " << params.output_file << '\n';
        return false;
      }

      // Write header
      outfile << "# Polytropic EOS: " << eos_data.name << '\n';
      outfile << "# k = " << std::scientific << std::setprecision(6) << eos_data.k
              << ", gamma = " << eos_data.gamma << '\n';
      outfile << "# mu_e = " << params.mu_e << '\n';
      outfile << "rho[g/cm^3],P[dyne/cm^2],h[erg/g],u[erg/cm^3]" << '\n';

      // Calculate EOS over density range
      double log_rho_min = std::log10(params.rho_min);
      double log_rho_max = std::log10(params.rho_max);
      double dlog_rho = (log_rho_max - log_rho_min) / (params.num_points - 1);

      for (int i = 0; i < params.num_points; i++) {
        double log_rho = log_rho_min + i * dlog_rho;
        double rho = std::pow(10.0, log_rho);

        // Use the custom mu_e if different from default
        double k_adjusted = eos_data.k;
        if (std::abs(params.mu_e - 2.0) > 1e-10) {
          // Adjust k for custom mu_e
          if (gas_type == PolytropicGasType::ELECTRON_NON_RELATIVISTIC) {
            k_adjusted = 1.0036e13 / std::pow(params.mu_e, 5.0 / 3.0);
          } else if (gas_type == PolytropicGasType::ELECTRON_RELATIVISTIC) {
            k_adjusted = 1.2435e15 / std::pow(params.mu_e, 4.0 / 3.0);
          }
        }

        // Calculate pressure using P = k * rho^gamma
        double pressure = k_adjusted * std::pow(rho, eos_data.gamma);

        // Calculate specific enthalpy h = (gamma * P) / ((gamma - 1) * rho)
        double specific_enthalpy = (eos_data.gamma * pressure) / ((eos_data.gamma - 1.0) * rho);

        // Calculate energy density u = P / (gamma - 1)
        double energy_density = pressure / (eos_data.gamma - 1.0);

        // Write to file
        outfile << std::scientific << std::setprecision(10) << rho << "," << pressure << ","
                << specific_enthalpy << "," << energy_density << '\n';
      }

      outfile.close();
      return true;

    } catch (const std::exception &e) {
      std::cerr << "Error in Polytropic EOS calculation: " << e.what() << '\n';
      return false;
    }
  }

  [[nodiscard]] std::string getType() const override {
    return "Polytropic EOS (" + type_name + ")";
  }
};

// Factory implementation
std::unique_ptr<EOSCalculator> EOSCalculatorFactory::createCalculator(EOSType type) {
  switch (type) {
  case EOSType::MAGNETIC_BPS:
    return std::make_unique<MagneticBPSCalculator>();
  case EOSType::NON_MAGNETIC_NPE_GAS:
    return std::make_unique<NonMagneticNPECalculator>();
  case EOSType::POLYTROPIC_ELECTRON_NON_REL:
    return std::make_unique<PolytropicEOSCalculator>(PolytropicGasType::ELECTRON_NON_RELATIVISTIC,
                                                     "Non-relativistic Electron Gas");
  case EOSType::POLYTROPIC_ELECTRON_REL:
    return std::make_unique<PolytropicEOSCalculator>(PolytropicGasType::ELECTRON_RELATIVISTIC,
                                                     "Relativistic Electron Gas");
  case EOSType::POLYTROPIC_NEUTRON_NON_REL:
    return std::make_unique<PolytropicEOSCalculator>(PolytropicGasType::NEUTRON_NON_RELATIVISTIC,
                                                     "Non-relativistic Neutron Gas");
  case EOSType::POLYTROPIC_NEUTRON_REL:
    return std::make_unique<PolytropicEOSCalculator>(PolytropicGasType::NEUTRON_RELATIVISTIC,
                                                     "Relativistic Neutron Gas");
  case EOSType::CUSTOM_EOS:
    throw std::runtime_error("Custom EOS not implemented yet");
  default:
    throw std::runtime_error("Unknown EOS type");
  }
}

// Function to calculate EOS based on parameters
bool calculateEOS(const std::string &eos_type, const EOSParameters &params) {
  // Validate parameters
  if (params.nB_min >= params.nB_max) {
    std::cerr << "Error: nB_min must be less than nB_max." << '\n';
    return false;
  }

  // Additional validation for polytropic EOS
  if (eos_type.find("POLYTROPIC") != std::string::npos) {
    if (params.rho_min >= params.rho_max) {
      std::cerr << "Error: rho_min must be less than rho_max for polytropic EOS." << '\n';
      return false;
    }
    if (params.mu_e <= 0.0) {
      std::cerr << "Error: mu_e must be positive for polytropic EOS." << '\n';
      return false;
    }
  }

  try {
    // Parse EOS type
    EOSType type = EOSType::MAGNETIC_BPS; // or any safe default; overwritten below
    if (eos_type == "MAGNETIC_BPS") {
      type = EOSType::MAGNETIC_BPS;
    } else if (eos_type == "NON_MAGNETIC_NPE_GAS") {
      type = EOSType::NON_MAGNETIC_NPE_GAS;
    } else if (eos_type == "POLYTROPIC_ELECTRON_NON_REL") {
      type = EOSType::POLYTROPIC_ELECTRON_NON_REL;
    } else if (eos_type == "POLYTROPIC_ELECTRON_REL") {
      type = EOSType::POLYTROPIC_ELECTRON_REL;
    } else if (eos_type == "POLYTROPIC_NEUTRON_NON_REL") {
      type = EOSType::POLYTROPIC_NEUTRON_NON_REL;
    } else if (eos_type == "POLYTROPIC_NEUTRON_REL") {
      type = EOSType::POLYTROPIC_NEUTRON_REL;
    } else {
      std::cerr << "Unknown EOS type: " << eos_type << '\n';
      std::cerr << "Supported types: MAGNETIC_BPS, NON_MAGNETIC_NPE_GAS, " << '\n';
      std::cerr << "                 POLYTROPIC_ELECTRON_NON_REL, POLYTROPIC_ELECTRON_REL," << '\n';
      std::cerr << "                 POLYTROPIC_NEUTRON_NON_REL, POLYTROPIC_NEUTRON_REL" << '\n';
      return false;
    }

    // Create calculator
    auto calculator = EOSCalculatorFactory::createCalculator(type);

    // Calculate EOS
    std::cout << "Calculating " << calculator->getType() << " EOS..." << '\n';
    std::cout << "Parameters:" << '\n';

    if (eos_type.find("POLYTROPIC") != std::string::npos) {
      std::cout << "  Density range: " << params.rho_min << " to " << params.rho_max << " g/cm^3"
                << '\n';
      std::cout << "  Number of points: " << params.num_points << '\n';
      std::cout << "  Mean molecular weight per electron: " << params.mu_e << '\n';
    } else {
      std::cout << "  nB range: " << params.nB_min << " to " << params.nB_max << '\n';
      std::cout << "  Number of points: " << params.num_points << '\n';
    }

    std::cout << "  Output file: " << params.output_file << '\n';

    if (type == EOSType::MAGNETIC_BPS) {
      std::cout << "  B ratio: " << params.B_ratio_electron << '\n';
    }
    if (type == EOSType::NON_MAGNETIC_NPE_GAS) {
      std::cout << "  Debug mode: " << (params.debug_mode ? "enabled" : "disabled") << '\n';
    }

    if (calculator->calculateEOS(params)) {
      std::cout << "Successfully calculated " << calculator->getType() << " EOS" << '\n';
      std::cout << "Results written to: " << params.output_file << '\n';
      return true;
    }
    std::cerr << "Failed to calculate EOS" << '\n';
    return false;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << '\n';
    return false;
  }
}

// Example usage in main (can be removed or kept for testing)
#ifdef TEST_EOS_CALCULATOR
int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <EOS_TYPE> [options]" << '\n';
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
