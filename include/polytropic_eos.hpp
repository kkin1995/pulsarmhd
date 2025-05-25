/**
 * @file polytropic_eos.hpp
 * @brief Polytropic Equation of State Calculator
 * 
 * This file provides a dedicated calculator for polytropic equations of state
 * that were previously hardcoded in the stellar structure solver. It implements
 * the four degenerate gas types used in white dwarf and neutron star modeling.
 * 
 * The polytropic EOS has the form: P = k * ρ^γ
 * 
 * Supported EOS types:
 * - Non-relativistic electron gas (white dwarf cores)
 * - Relativistic electron gas (massive white dwarfs)  
 * - Non-relativistic neutron gas (low-density neutron star regions)
 * - Relativistic neutron gas (high-density neutron star cores)
 * 
 * @author Karan Amit Kinariwala
 * @date 2025-05-25
 * 
 * @example Basic Usage
 * ```cpp
 * PolytropicEOS eos;
 * PolytropicEOSData data = eos.getEOSParameters(PolytropicGasType::ELECTRON_RELATIVISTIC);
 * 
 * // Calculate pressure for given density
 * double pressure = eos.calculatePressure(1e15, data.k, data.gamma);
 * ```
 * 
 * @example Generate EOS Table
 * ```cpp
 * PolytropicEOSParameters params;
 * params.gas_type = PolytropicGasType::NEUTRON_RELATIVISTIC;
 * params.rho_min = 1e14;
 * params.rho_max = 1e18;
 * params.num_points = 1000;
 * params.output_file = "neutron_eos.csv";
 * 
 * PolytropicEOS eos;
 * bool success = eos.generateEOSTable(params);
 * ```
 */

#ifndef POLYTROPIC_EOS_HPP
#define POLYTROPIC_EOS_HPP

#include <string>
#include <vector>
#include <memory>
#include <tuple>

/**
 * @brief Types of polytropic degenerate gases
 * 
 * These correspond to the four hardcoded types from the original
 * stellar structure code, representing different physical regimes
 * of degenerate matter in compact objects.
 */
enum class PolytropicGasType {
    ELECTRON_NON_RELATIVISTIC,  ///< Non-relativistic electron gas (white dwarf cores)
    ELECTRON_RELATIVISTIC,      ///< Relativistic electron gas (massive white dwarfs)
    NEUTRON_NON_RELATIVISTIC,   ///< Non-relativistic neutron gas (low-density neutron star)
    NEUTRON_RELATIVISTIC        ///< Relativistic neutron gas (high-density neutron star)
};

/**
 * @brief Structure holding polytropic EOS parameters
 * 
 * Contains the k and γ values that define the polytropic relation P = k ρ^γ
 */
struct PolytropicEOSData {
    std::string name;           ///< Descriptive name of the EOS
    double k;                   ///< EOS proportionality constant
    double gamma;               ///< Polytropic index
    double mu_e;                ///< Mean molecular weight per electron (where applicable)
    
    PolytropicEOSData() : name(""), k(0.0), gamma(0.0), mu_e(2.0) {}
    
    PolytropicEOSData(const std::string& n, double k_val, double gamma_val, double mu_e_val = 2.0) 
        : name(n), k(k_val), gamma(gamma_val), mu_e(mu_e_val) {}
};

/**
 * @brief Parameters for generating EOS tables
 */
struct PolytropicEOSParameters {
    PolytropicGasType gas_type;     ///< Type of polytropic gas
    double rho_min;                 ///< Minimum density (g/cm³)
    double rho_max;                 ///< Maximum density (g/cm³)
    int num_points;                 ///< Number of points in the table
    std::string output_file;        ///< Output file path
    bool use_log_spacing;           ///< Use logarithmic spacing for density grid
    bool output_log_values;         ///< Output logarithmic values
    
    /**
     * @brief Constructor with default values
     */
    PolytropicEOSParameters() : 
        gas_type(PolytropicGasType::NEUTRON_RELATIVISTIC),
        rho_min(1e10),
        rho_max(1e18),
        num_points(1000),
        output_file("polytropic_eos.csv"),
        use_log_spacing(true),
        output_log_values(true) {}
};

/**
 * @brief Main class for polytropic EOS calculations
 * 
 * This class provides methods to:
 * - Get EOS parameters for different gas types
 * - Calculate pressure for given density
 * - Generate EOS tables
 * - Convert between different units
 */
class PolytropicEOS {
public:
    /**
     * @brief Default constructor
     */
    PolytropicEOS();

    /**
     * @brief Get EOS parameters for a specific gas type
     * 
     * Returns the k and γ values that define the polytropic EOS.
     * These are the same values that were hardcoded in the original
     * stellar structure solver.
     * 
     * @param type The type of polytropic gas
     * @return PolytropicEOSData containing k, γ, and descriptive information
     */
    PolytropicEOSData getEOSParameters(PolytropicGasType type) const;

    /**
     * @brief Calculate pressure for given density using polytropic EOS
     * 
     * Computes P = k * ρ^γ
     * 
     * @param density Density in g/cm³
     * @param k EOS proportionality constant
     * @param gamma Polytropic index
     * @return Pressure in dyne/cm²
     */
    double calculatePressure(double density, double k, double gamma) const;

    /**
     * @brief Calculate density for given pressure using polytropic EOS
     * 
     * Computes ρ = (P/k)^(1/γ)
     * 
     * @param pressure Pressure in dyne/cm²
     * @param k EOS proportionality constant
     * @param gamma Polytropic index
     * @return Density in g/cm³
     */
    double calculateDensity(double pressure, double k, double gamma) const;

    /**
     * @brief Generate a table of EOS values
     * 
     * Creates a table of density-pressure pairs according to the
     * specified parameters and writes to file.
     * 
     * @param params Parameters controlling the table generation
     * @return true if successful, false otherwise
     */
    bool generateEOSTable(const PolytropicEOSParameters& params) const;

    /**
     * @brief Get all available gas types with their parameters
     * 
     * @return Vector of all supported gas types with their EOS data
     */
    std::vector<std::pair<PolytropicGasType, PolytropicEOSData>> getAllEOSTypes() const;

    /**
     * @brief Convert gas type enum to string
     * 
     * @param type The gas type
     * @return String representation of the gas type
     */
    static std::string gasTypeToString(PolytropicGasType type);

    /**
     * @brief Convert string to gas type enum
     * 
     * @param str String representation of gas type
     * @return Gas type enum
     * @throw std::invalid_argument if string is not recognized
     */
    static PolytropicGasType stringToGasType(const std::string& str);

private:
    /**
     * @brief Create density grid for EOS table
     * 
     * @param rho_min Minimum density
     * @param rho_max Maximum density
     * @param num_points Number of points
     * @param use_log_spacing Whether to use logarithmic spacing
     * @return Vector of density values
     */
    std::vector<double> createDensityGrid(double rho_min, double rho_max, 
                                         int num_points, bool use_log_spacing) const;

    /**
     * @brief Write EOS table to file
     * 
     * @param filename Output filename
     * @param densities Vector of density values
     * @param pressures Vector of pressure values
     * @param eos_data EOS parameters used
     * @param output_log_values Whether to output logarithmic values
     * @return true if successful, false otherwise
     */
    bool writeEOSToFile(const std::string& filename,
                       const std::vector<double>& densities,
                       const std::vector<double>& pressures,
                       const PolytropicEOSData& eos_data,
                       bool output_log_values) const;
};

/**
 * @brief Standalone function to calculate EOS table
 * 
 * Convenience function that creates a PolytropicEOS instance and
 * generates an EOS table with the given parameters.
 * 
 * @param gas_type_str String representation of gas type
 * @param params EOS table generation parameters
 * @return true if successful, false otherwise
 */
bool calculatePolytropicEOS(const std::string& gas_type_str, 
                           const PolytropicEOSParameters& params);

#endif // POLYTROPIC_EOS_HPP 