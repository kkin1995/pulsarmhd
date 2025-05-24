/**
 * @file eos_calculator.hpp
 * @brief Equation of State (EOS) Calculator Abstraction Layer
 * 
 * This file provides an abstraction layer for calculating different types of
 * Equations of State (EOS) for neutron stars. It implements the Factory pattern
 * to create appropriate EOS calculators based on the type of EOS required.
 * 
 * The system supports:
 * - Magnetic BPS EOS
 * - Non-magnetic NPE Gas EOS
 * - (Future) Custom EOS
 * 
 * @author Karan Amit Kinariwala
 * @date 2025-05-23
 * 
 * @example Basic Usage
 * ```cpp
 * // Create parameters with default values
 * EOSParameters params;
 * 
 * // Calculate Magnetic BPS EOS
 * bool success = calculateEOS("MAGNETIC_BPS", params);
 * ```
 * 
 * @example Custom Parameters
 * ```cpp
 * // Create parameters with custom values
 * EOSParameters params;
 * params.nB_min = 1e22;          // Minimum baryon density
 * params.nB_max = 1e37;          // Maximum baryon density
 * params.num_points = 500;       // Higher resolution
 * params.output_file = "custom_eos.csv";
 * 
 * // For Magnetic BPS
 * params.B_ratio_electron = 0.01;
 * params.atomic_mass_file = "data/custom_atomic_masses.csv";
 * params.rel_tolerance = 1e-8;   // Stricter tolerance
 * params.abs_tolerance = 1e-10;
 * 
 * // For Non-magnetic NPE Gas
 * params.debug_mode = true;      // Enable debug output
 * 
 * // Calculate EOS
 * bool success = calculateEOS("MAGNETIC_BPS", params);
 * ```
 * 
 * @example Error Handling
 * ```cpp
 * try {
 *     EOSParameters params;
 *     params.nB_min = 1e37;  // Invalid: min > max
 *     params.nB_max = 1e22;
 *     
 *     bool success = calculateEOS("MAGNETIC_BPS", params);
 *     if (!success) {
 *         std::cerr << "EOS calculation failed" << std::endl;
 *     }
 * } catch (const std::exception& e) {
 *     std::cerr << "Error: " << e.what() << std::endl;
 * }
 * ```
 * 
 * @example Command Line Usage
 * ```bash
 * # Calculate Magnetic BPS EOS with default parameters
 * ./eos_calculator MAGNETIC_BPS
 * 
 * # Calculate with custom parameters
 * ./eos_calculator MAGNETIC_BPS --nB-min 1e22 --nB-max 1e37 --num-points 500 --output custom_eos.csv --B-ratio 0.01
 * 
 * # Calculate Non-magnetic NPE Gas with debug output
 * ./eos_calculator NON_MAGNETIC_NPE_GAS --debug
 * ```
 */

#ifndef EOS_CALCULATOR_HPP
#define EOS_CALCULATOR_HPP

#include <string>
#include <vector>
#include <memory>

// Forward declarations
class MagneticBPSEOS;
class NonMagneticNPEGas;
class PolytropicEOS;

/**
 * @brief Enumeration of supported EOS types
 */
enum class EOSType {
    MAGNETIC_BPS,           ///< Magnetic Baym-Pethick-Sutherland EOS
    NON_MAGNETIC_NPE_GAS,   ///< Non-magnetic Neutron-Proton-Electron Gas EOS
    POLYTROPIC_ELECTRON_NON_REL,  ///< Non-relativistic electron polytropic EOS
    POLYTROPIC_ELECTRON_REL,      ///< Relativistic electron polytropic EOS
    POLYTROPIC_NEUTRON_NON_REL,   ///< Non-relativistic neutron polytropic EOS
    POLYTROPIC_NEUTRON_REL,       ///< Relativistic neutron polytropic EOS
    CUSTOM_EOS             ///< (Future) Custom EOS implementation
};

/**
 * @brief Structure holding parameters for EOS calculations
 * 
 * This structure contains all parameters needed for EOS calculations,
 * including common parameters and type-specific parameters.
 */
struct EOSParameters {
    // Common parameters
    double nB_min;          ///< Minimum baryon density (cm^-3)
    double nB_max;          ///< Maximum baryon density (cm^-3)
    int num_points;         ///< Number of points in the calculation
    std::string output_file;///< Output file path for results
    
    // Magnetic BPS specific parameters
    std::string atomic_mass_file;  ///< Path to atomic mass data file
    double B_ratio_electron;       ///< Magnetic field strength ratio for electrons
    double rel_tolerance;          ///< Relative tolerance for calculations
    double abs_tolerance;          ///< Absolute tolerance for calculations
    
    // Non-magnetic NPE gas specific parameters
    bool debug_mode;               ///< Enable debug output
    
    // Polytropic EOS specific parameters
    double rho_min;                ///< Minimum density for polytropic EOS (g/cm^3)
    double rho_max;                ///< Maximum density for polytropic EOS (g/cm^3)
    double mu_e;                   ///< Mean molecular weight per electron
    
    /**
     * @brief Constructor with default values
     * 
     * Sets reasonable default values for all parameters:
     * - nB range: 1e22 to 1e37 cm^-3
     * - 250 calculation points
     * - Default output file: eos_output.csv
     * - Default atomic mass file: data/atomic_masses.csv
     * - Default B ratio: 0.01
     * - Default tolerances: 1e-6 (relative), 1e-8 (absolute)
     * - Debug mode disabled by default
     */
    EOSParameters() : 
        nB_min(1e22),
        nB_max(1e37),
        num_points(250),
        output_file("eos_output.csv"),
        atomic_mass_file("data/atomic_masses.csv"),
        B_ratio_electron(0.01),
        rel_tolerance(1e-6),
        abs_tolerance(1e-8),
        debug_mode(false),
        rho_min(1e6),
        rho_max(1e15),
        mu_e(2.0) {}
};

/**
 * @brief Abstract base class for EOS calculators
 * 
 * This class defines the interface that all EOS calculators must implement.
 * It uses the Strategy pattern to allow different EOS calculation methods
 * to be used interchangeably.
 */
class EOSCalculator {
public:
    virtual ~EOSCalculator() = default;

    /**
     * @brief Calculate the EOS using the provided parameters
     * @param params The parameters for the EOS calculation
     * @return true if calculation was successful, false otherwise
     */
    virtual bool calculateEOS(const EOSParameters& params) = 0;

    /**
     * @brief Get the type of EOS calculator
     * @return String describing the EOS type
     */
    virtual std::string getType() const = 0;
};

/**
 * @brief Factory class for creating EOS calculators
 * 
 * This class implements the Factory pattern to create appropriate
 * EOS calculator instances based on the requested type.
 */
class EOSCalculatorFactory {
public:
    /**
     * @brief Create an EOS calculator of the specified type
     * @param type The type of EOS calculator to create
     * @return Unique pointer to the created calculator
     * @throw std::runtime_error if the type is not supported
     */
    static std::unique_ptr<EOSCalculator> createCalculator(EOSType type);
};

/**
 * @brief Main function to calculate EOS
 * 
 * This is the primary interface for calculating EOS. It handles:
 * - Parameter validation
 * - EOS type parsing
 * - Calculator creation
 * - Progress reporting
 * - Error handling
 * 
 * @param eos_type String specifying the type of EOS to calculate
 * @param params Parameters for the EOS calculation
 * @return true if calculation was successful, false otherwise
 * 
 * @example
 * ```cpp
 * EOSParameters params;
 * params.nB_min = 1e22;
 * params.nB_max = 1e37;
 * params.num_points = 250;
 * params.output_file = "magnetic_bps_eos.csv";
 * params.B_ratio_electron = 0.01;
 * 
 * bool success = calculateEOS("MAGNETIC_BPS", params);
 * ```
 */
bool calculateEOS(const std::string& eos_type, const EOSParameters& params);

#endif // EOS_CALCULATOR_HPP 