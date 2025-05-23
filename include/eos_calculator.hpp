#ifndef EOS_CALCULATOR_HPP
#define EOS_CALCULATOR_HPP

#include <string>
#include <vector>
#include <memory>

// Forward declarations
class MagneticBPSEOS;
class NonMagneticNPEGas;

// Enum for different EOS types
enum class EOSType {
    MAGNETIC_BPS,
    NON_MAGNETIC_NPE_GAS,
    CUSTOM_EOS
};

// Structure to hold EOS calculation parameters
struct EOSParameters {
    // Common parameters
    double nB_min;
    double nB_max;
    int num_points;
    std::string output_file;
    
    // Magnetic BPS specific parameters
    std::string atomic_mass_file;
    double B_ratio_electron;
    double rel_tolerance;
    double abs_tolerance;
    
    // Non-magnetic NPE gas specific parameters
    bool debug_mode;
    
    // Constructor with default values
    EOSParameters() : 
        nB_min(1e22),
        nB_max(1e37),
        num_points(250),
        output_file("eos_output.csv"),
        atomic_mass_file("data/atomic_masses.csv"),
        B_ratio_electron(0.01),
        rel_tolerance(1e-6),
        abs_tolerance(1e-8),
        debug_mode(false) {}
};

// Abstract base class for EOS calculators
class EOSCalculator {
public:
    virtual ~EOSCalculator() = default;
    virtual bool calculateEOS(const EOSParameters& params) = 0;
    virtual std::string getType() const = 0;
};

// Factory class to create EOS calculators
class EOSCalculatorFactory {
public:
    static std::unique_ptr<EOSCalculator> createCalculator(EOSType type);
};

// Main function to calculate EOS
bool calculateEOS(const std::string& eos_type, const EOSParameters& params);

#endif // EOS_CALCULATOR_HPP 