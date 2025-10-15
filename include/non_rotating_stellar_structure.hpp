// non_rotating_stellar_structure.hpp
#ifndef NON_ROTATING_STELLAR_STRUCTURE_HPP
#define NON_ROTATING_STELLAR_STRUCTURE_HPP

/**
 * @file non_rotating_stellar_structure.hpp
 * @brief Header file for non-rotating stellar structure simulations.
 * @author Karan Kinariwala
 * @date 2025-05-25
 *
 * @details
 * This file defines constants, enumerations, and function declarations used for
 * simulating the structure of non-rotating stars, such as white dwarfs and neutron stars.
 * The functions solve hydrostatic equilibrium equations for both Newtonian
 * and relativistic frameworks, using different equations of state (EOS).
 *
 * @section features Features
 * @li Implements Newtonian and Tolman-Oppenheimer-Volkoff (TOV) equations for stellar structure
 * @li Supports equations of state for various degenerate gases:
 *     - Non-relativistic electron gas
 *     - Relativistic electron gas
 *     - Non-relativistic neutron gas
 *     - Relativistic neutron gas
 * @li Provides functionality for generating mass-radius relations
 *
 * @section constants Constants
 * @li G: Gravitational constant in CGS units
 * @li c: Speed of light in CGS units
 *
 * @section enums Enumerations
 * @li PolytropicGasType: Enum representing the type of polytropic gas (from polytropic_eos.hpp)
 *
 * @section functions Functions
 * @li newtonian: Computes derivatives for Newtonian hydrostatic equilibrium
 * @li tolman_oppenheimer_volkoff_derivatives: Computes derivatives for relativistic stellar
 * structure
 * @li set_eos_parameters: Sets the EOS parameters (k and gamma) based on the gas type
 * @li get_filename: Generates filenames for output data based on EOS and central density
 *
 * @section example Example Usage
 * @code{.cpp}
 * std::string eos_name;
 * double k, gamma;
 * auto eos_data = calculator.getEOSParameters(PolytropicGasType::ELECTRON_RELATIVISTIC);
 *
 * std::vector<double> state = {log_m, log_p};
 * std::vector<double> derivatives = tolman_oppenheimer_volkoff_derivatives(log_r, state, k, gamma);
 * @endcode
 *
 * @see rk4.hpp for numerical integration using Runge-Kutta methods.
 */

#include "polytropic_eos.hpp"

#include <string>
#include <vector>

// GSL includes for spline-based EOS support
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// Physical Constants in CGS

/**
 * @brief Gravitational constant in CGS units.
 */
inline constexpr double G = 6.67430e-8; // gravitational constant in CGS

/**
 * @brief Speed of light in CGS units.
 */
inline constexpr double c = 2.99792458e10; // speed of light in CGS units

/**
 * @brief Polytropic EOS types are now handled by the PolytropicEOS class.
 *
 * The stellar structure solver now uses the centralized polytropic EOS calculator
 * from polytropic_eos.hpp, which provides the same gas types with enhanced
 * functionality and consistency across the codebase.
 *
 * @see PolytropicGasType in polytropic_eos.hpp for available EOS types
 * @see PolytropicEOS class for EOS parameter calculation
 */

struct TovResult {
  int steps;
  double log10_m_surface;
  double log10_r_surface;
  bool found_surface;
};

/**
 * @brief Solves the stellar structure equations for non-rotating compact objects.
 *
 * @details
 * This function integrates the Tolman-Oppenheimer-Volkoff (TOV) equations to compute
 * the mass-radius relationship for non-rotating relativistic stars such as white dwarfs
 * and neutron stars. The integration proceeds from a small inner radius outward until
 * the stellar surface is reached (defined by pressure dropping below a threshold).
 *
 * **Physical Framework:**
 * The function solves the TOV equations in logarithmic coordinates:
 * \f{eqnarray*}{
 * \frac{d(\log m)}{d(\log r)} &=& \frac{4 \pi r^3 \rho}{m} \\
 * \frac{d(\log P)}{d(\log r)} &=& \left( - \frac{G m \rho}{P r} \right)
 *                                  \cdot \left( 1 + \frac{P}{\rho c^2} \right)
 *                                  \cdot \left( 1 + \frac{4 \pi P r^3}{m c^2} \right)
 *                                  \cdot \left( 1 - \frac{2 G m}{r c^2} \right)^{-1}
 * \f}
 *
 * **Equation of State:**
 * Uses polytropic EOS: \f$P = k \rho^\gamma\f$ where k and γ are automatically
 * determined from the gas type. Supported types include:
 * - Non-relativistic electron gas (γ = 5/3): White dwarf cores
 * - Relativistic electron gas (γ = 4/3): Massive white dwarfs
 * - Non-relativistic neutron gas (γ = 5/3): Low-density neutron star regions
 * - Relativistic neutron gas (γ = 4/3): High-density neutron star cores
 *
 * **Initial Conditions:**
 * At the starting radius r_start, the enclosed mass is computed assuming
 * uniform density: \f$m_0 = \frac{4}{3}\pi r_{start}^3 \rho_c\f$
 * The initial pressure is: \f$P_0 = k \rho_c^\gamma\f$
 *
 * **Integration Method:**
 * Uses 4th-order Runge-Kutta integration with adaptive step size control.
 * Integration terminates when:
 * - Pressure drops below 10⁻⁸ dyne/cm² (surface condition)
 * - Non-finite values are encountered (numerical instability)
 * - Maximum radius is reached
 *
 * **Output Data:**
 * Writes CSV file containing log₁₀(r), log₁₀(m), log₁₀(P) at each integration step.
 * Filename format: "data/[eos_name]_rhoc_[density].csv"
 *
 * @param[in] eos_type Type of polytropic gas defining the equation of state
 * @param[in] rho_c Central density in g/cm³ (typically 10⁶ - 10¹⁸)
 * @param[in] r_start Starting radius for integration in cm (typically ~10 cm)
 * @param[in] r_end Maximum radius for integration in cm (typically ~10⁶ cm)
 * @param[in] dlogr Step size in log₁₀(r) for integration (typically 0.01 - 0.1)
 * @param[in] mu_e Mean molecular weight per electron for electron gases (default: 2.0)
 *                 Only affects electron gas types; modifies k parameter as:
 *                 - Non-rel: k → k/μₑ^(5/3)
 *                 - Relativistic: k → k/μₑ^(4/3)
 *
 * @return std::tuple<int, double, double> containing:
 *         - Integration steps completed before reaching surface
 *         - Final log₁₀(mass) in grams at stellar surface
 *         - Final log₁₀(radius) in cm at stellar surface
 *
 * @pre rho_c > 0 (positive central density)
 * @pre r_start > 0 and r_end > r_start (valid radius range)
 * @pre dlogr > 0 (positive step size)
 * @pre mu_e > 0 (positive mean molecular weight)
 *
 * @post Creates output file in "data/" directory
 * @post Returns physically meaningful mass for successful integration
 *
 * @note For white dwarfs: typical rho_c ~ 10⁶ g/cm³, final mass ~ 0.6 M☉
 * @note For neutron stars: typical rho_c ~ 10¹⁵ g/cm³, final mass ~ 1.4 M☉
 * @note Custom μₑ allows modeling different compositions (e.g., He: μₑ=4, C: μₑ=12)
 *
 * @warning Requires "data/" directory to exist for output file creation
 * @warning May fail for extreme densities causing numerical instabilities
 * @warning Integration may not converge for inappropriate step sizes
 *
 * **Example Usage:**
 * @code{.cpp}
 * // Calculate white dwarf structure
 * auto [steps, log_mass] = non_rotating_stellar_structure(
 *     PolytropicGasType::ELECTRON_NON_RELATIVISTIC,
 *     1e6,    // Central density: 10⁶ g/cm³
 *     10.0,   // Start radius: 10 cm
 *     1e6,    // End radius: 10⁶ cm
 *     0.05    // Step size in log(r)
 * );
 *
 * double mass_grams = std::pow(10.0, log_mass);
 * double mass_solar = mass_grams / 1.989e33;  // Convert to solar masses
 * std::cout << "White dwarf mass: " << mass_solar << " M☉" << std::endl;
 *
 * // Calculate neutron star with custom composition
 * auto [ns_steps, ns_log_mass] = non_rotating_stellar_structure(
 *     PolytropicGasType::NEUTRON_RELATIVISTIC,
 *     5e14,   // Central density: 5×10¹⁴ g/cm³
 *     10.0,   // Start radius: 10 cm
 *     1e6,    // End radius: 10⁶ cm
 *     0.02,   // Smaller step for accuracy
 *     2.0     // Standard μₑ for neutron matter
 * );
 * @endcode
 *
 * @see tolman_oppenheimer_volkoff_derivatives() for the TOV equation implementation
 * @see PolytropicEOS::getEOSParameters() for EOS parameter calculation
 * @see get_filename() for output file naming convention
 */
[[nodiscard]] std::tuple<int, double, double>
non_rotating_stellar_structure(PolytropicGasType eos_type, double rho_c, double r_start,
                               double r_end, double dlogr, double mu_e = 2.0);

/**
 * @brief Solves stellar structure equations for non-rotating compact objects using spline-based
 * EOS.
 *
 * @details
 * This function integrates the Tolman-Oppenheimer-Volkoff (TOV) equations for relativistic
 * stellar structure using a tabulated equation of state provided via GSL splines. It extends
 * the existing TOV solver to support realistic, complex equations of state that cannot be
 * represented by simple polytropic relationships.
 *
 * **Physical Framework:**
 * The function solves the same TOV equations as the polytropic version:
 * \f{eqnarray*}{
 * \frac{d(\log m)}{d(\log r)} &=& \frac{4 \pi r^3 \rho}{m} \\
 * \frac{d(\log P)}{d(\log r)} &=& \left( - \frac{G m \rho}{P r} \right)
 *                                  \cdot \left( 1 + \frac{P}{\rho c^2} \right)
 *                                  \cdot \left( 1 + \frac{4 \pi P r^3}{m c^2} \right)
 *                                  \cdot \left( 1 - \frac{2 G m}{r c^2} \right)^{-1}
 * \f}
 *
 * **Equation of State:**
 * Uses tabulated EOS data via GSL spline interpolation. The density is obtained from
 * pressure using the inverse spline: \f$\rho = f^{-1}(P)\f$ where \f$f\f$ is the
 * tabulated pressure-density relationship.
 *
 * **Initial Conditions:**
 * At the starting radius r_start, the enclosed mass is computed assuming
 * uniform density: \f$m_0 = \frac{4}{3}\pi r_{start}^3 \rho_c\f$
 * The initial pressure is obtained from the forward spline: \f$P_0 = f(\rho_c)\f$
 *
 * **Integration Method:**
 * Uses 4th-order Runge-Kutta integration with optional adaptive step size control.
 * Adaptive stepping adjusts the step size based on pressure gradient magnitude to
 * maintain numerical stability in regions of rapid change.
 *
 * **Surface Detection:**
 * Integration terminates when:
 * - Pressure drops below configurable threshold (default: 10⁻⁸ dyne/cm²)
 * - Non-finite values are encountered (numerical instability)
 * - Maximum radius is reached
 * - Pressure falls outside EOS validity range
 *
 * **Output Data:**
 * Optionally writes CSV file containing log₁₀(r), log₁₀(m), log₁₀(P) at each integration step.
 * Output filename can be specified or auto-generated based on central density.
 *
 * @param[in] spline_inv GSL spline object for inverse EOS: \f$\log\rho = f^{-1}(\log P)\f$
 * @param[in] acc_inv GSL interpolation accelerator for inverse spline (for efficiency)
 * @param[in] min_logP Minimum log₁₀(pressure) for EOS validity range (dyne/cm²)
 * @param[in] max_logP Maximum log₁₀(pressure) for EOS validity range (dyne/cm²)
 * @param[in] rho_c Central density in g/cm³ (typically 10¹⁵ - 10¹⁸ for neutron stars)
 * @param[in] r_start Starting radius for integration in cm (typically ~1-10 cm)
 * @param[in] r_end Maximum radius for integration in cm (typically ~10⁶ cm)
 * @param[in] base_dlogr Base step size in log₁₀(r) for integration (typically 0.0001)
 * @param[in] use_adaptive_stepping Enable adaptive step size control (default: true)
 * @param[in] pressure_threshold Surface pressure threshold in dyne/cm² (default: 10⁻⁸)
 * @param[in] output_filename Optional output CSV filename (empty = no file output)
 *
 * @return std::tuple<int, double, double> containing:
 *         - Integration steps completed before reaching surface
 *         - Final log₁₀(mass) in grams at stellar surface
 *         - Final log₁₀(radius) in cm at stellar surface
 *
 * @pre spline_inv != nullptr (valid GSL spline object)
 * @pre acc_inv != nullptr (valid GSL accelerator object)
 * @pre min_logP < max_logP (valid pressure range)
 * @pre rho_c > 0 (positive central density)
 * @pre r_start > 0 and r_end > r_start (valid radius range)
 * @pre base_dlogr > 0 (positive step size)
 * @pre pressure_threshold > 0 (positive surface threshold)
 *
 * @post Returns physically meaningful mass and radius for successful integration
 * @post Creates output file if filename provided
 * @post Spline objects remain unchanged (read-only access)
 *
 * @note Designed for neutron star EOS: magnetic BPS, BBP, unified models, etc.
 * @note Adaptive stepping crucial for complex EOS with phase transitions
 * @note Central density typically scanned over range to generate mass-radius relation
 *
 * @warning Requires valid, initialized GSL spline objects
 * @warning May fail for extreme densities outside EOS validity range
 * @warning Integration may not converge if step size too large for EOS complexity
 * @warning GSL spline evaluation can be expensive; consider caching for repeated calls
 *
 * **Example Usage:**
 * @code{.cpp}
 * // Assuming EOS data loaded and splines initialized
 * gsl_spline* spline_inv = ...; // ρ(P) spline
 * gsl_interp_accel* acc_inv = ...;
 * double min_logP = 15.0, max_logP = 35.0; // EOS validity range
 *
 * // Calculate neutron star structure with realistic EOS
 * auto [steps, log_mass, log_radius] = non_rotating_stellar_structure_spline(
 *     spline_inv, acc_inv, min_logP, max_logP,
 *     1e15,      // Central density: 10¹⁵ g/cm³
 *     1.0,       // Start radius: 1 cm
 *     1e6,       // End radius: 10⁶ cm
 *     0.0001,    // Base step size
 *     true,      // Enable adaptive stepping
 *     1e-8,      // Surface pressure threshold
 *     "ns_structure.csv"  // Output file
 * );
 *
 * double mass_grams = std::pow(10.0, log_mass);
 * double radius_cm = std::pow(10.0, log_radius);
 * double mass_solar = mass_grams / 1.989e33;
 * double radius_km = radius_cm / 1e5;
 *
 * std::cout << "Neutron star: " << mass_solar << " M☉, "
 *           << radius_km << " km" << std::endl;
 *
 * // Scan over central densities for mass-radius relation
 * std::vector<double> densities = {1e15, 5e15, 1e16, 5e16, 1e17};
 * for (double rho_c : densities) {
 *     auto result = non_rotating_stellar_structure_spline(
 *         spline_inv, acc_inv, min_logP, max_logP, rho_c,
 *         1.0, 1e6, 0.0001, true, 1e-8, ""  // No file output
 *     );
 *     // Process results...
 * }
 * @endcode
 *
 * @see non_rotating_stellar_structure() for polytropic EOS version
 * @see tolman_oppenheimer_volkoff_derivatives() for TOV equation details
 * @see GSL documentation for spline initialization and management
 */
[[nodiscard]] TovResult non_rotating_stellar_structure_spline(
    const gsl_spline *spline_inv, gsl_interp_accel *acc_inv, double min_logP, double max_logP,
    double rho_c, double r_start, double r_end, double base_dlogr, bool use_adaptive_stepping,
    double pressure_threshold, const std::string &output_filename,
    const gsl_spline *spline_eps = nullptr, gsl_interp_accel *acc_eps = nullptr);

/**
 * @brief Computes the mass and pressure derivatives using the Newtonian equations of stellar
 * structure.
 *
 * @details
 * This function calculates the derivatives of logarithmic mass (\f$\log m\f$) and
 * logarithmic pressure (\f$\log P\f$) with respect to the logarithmic radius (\f$\log r\f$)
 * using the Newtonian framework. It assumes a polytropic equation of state (EOS)
 * defined by the constants \p k and \p gamma.
 *
 * The derivatives are computed as:
 * \f{eqnarray*}{
 * \frac{d(\log m)}{d(\log r)} &=& \frac{4\pi r^3 \rho}{m} \\
 * \frac{d(\log P)}{d(\log r)} &=& -\frac{G m \rho}{P r}
 * \f}
 *
 * where:
 * \f{description}
 * \item[$m$] Mass enclosed within radius $r$
 * \item[$\rho$] Density, calculated from pressure and EOS constants
 * \item[$P$] Pressure at radius $r$
 * \item[$G$] Gravitational constant
 * \f}
 *
 * @param[in] log_r Logarithmic radius (\f$\log_{10}(r)\f$) in centimeters
 * @param[in] state Vector containing:
 *                  \parblock
 *                  - \p state[0]: Logarithmic mass (\f$\log_{10}(m)\f$) in grams
 *                  - \p state[1]: Logarithmic pressure (\f$\log_{10}(P)\f$) in dyne/cm²
 *                  \endparblock
 * @param[in] k EOS proportionality constant (\f$k\f$), relating pressure and density
 * @param[in] gamma Polytropic index (\f$\gamma\f$), defining the EOS relationship:
 *                  \f$P = k \rho^\gamma\f$
 *
 * @return Vector of size 2 containing:
 *         \parblock
 *         - \p [0]: Derivative of logarithmic mass (\f$\frac{d(\log m)}{d(\log r)}\f$)
 *         - \p [1]: Derivative of logarithmic pressure (\f$\frac{d(\log P)}{d(\log r)}\f$)
 *         \endparblock
 *
 * @pre \p m must be greater than zero (\f$m > 0\f$) to avoid division by zero
 * @pre The EOS must be valid for the given constants \p k and \p gamma
 *
 * @note This function does not include relativistic corrections
 * @note Suitable for non-relativistic stellar structures
 *
 * @warning May encounter floating-point errors if inputs are not physically meaningful
 *
 * Example usage:
 * @code{.cpp}
 * double log_r = 0.0; // Logarithmic radius
 * std::vector<double> state = {1.0, 15.0}; // log10(m), log10(P)
 * double k = 5.3802e9; // EOS constant for neutron gas
 * double gamma = 5.0 / 3.0; // Polytropic index
 *
 * auto derivatives = newtonian(log_r, state, k, gamma);
 * std::cout << "dlogm_dlogr: " << derivatives[0] << std::endl;
 * std::cout << "dlogP_dlogr: " << derivatives[1] << std::endl;
 * @endcode
 *
 * @see tolman_oppenheimer_volkoff_derivatives() for the relativistic version
 */
[[nodiscard]] std::vector<double> newtonian(double log_r, const std::vector<double> &state,
                                            double k, double gamma);

/**
 * @brief Computes the mass and pressure derivatives using the TOV (Tolman-Oppenheimer-Volkoff)
 * equations.
 *
 * @details
 * This function calculates the derivatives of logarithmic mass (\f$\log m\f$) and logarithmic
 * pressure
 * (\f$\log P\f$) with respect to the logarithmic radius (\f$\log r\f$) in the relativistic
 * framework. The TOV equations account for general relativistic corrections, providing a more
 * accurate description of compact stellar objects like neutron stars.
 *
 * The TOV equations are:
 * \f{eqnarray*}{
 * \frac{d(\log m)}{d(\log r)} &=& \frac{4 \pi r^3 \rho}{m} \\
 * \frac{d(\log P)}{d(\log r)} &=& \left( - \frac{G m \rho}{P r} \right)
 *                                  \cdot \left( 1 + \frac{P}{\rho c^2} \right)
 *                                  \cdot \left( 1 + \frac{4 \pi P r^3}{m c^2} \right)
 *                                  \cdot \left( 1 - \frac{2 G m}{r c^2} \right)^{-1}
 * \f}
 *
 * The relativistic corrections appear as four factors:
 * \f{description}
 * \item[Factor 1] Newtonian term: $-\frac{G m \rho}{P r}$
 * \item[Factor 2] Internal energy: $1 + \frac{P}{\rho c^2}$
 * \item[Factor 3] Pressure effect: $1 + \frac{4 \pi P r^3}{m c^2}$
 * \item[Factor 4] Metric factor: $\left(1 - \frac{2 G m}{r c^2}\right)^{-1}$
 * \f}
 *
 * @param[in] log_r Logarithmic radius (\f$\log_{10}(r)\f$) in centimeters
 * @param[in] state Vector containing:
 *                  \parblock
 *                  - \p state[0]: Logarithmic mass (\f$\log_{10}(m)\f$) in grams
 *                  - \p state[1]: Logarithmic pressure (\f$\log_{10}(P)\f$) in dyne/cm²
 *                  \endparblock
 * @param[in] k EOS proportionality constant (\f$k\f$), relating pressure and density
 * @param[in] gamma Polytropic index (\f$\gamma\f$), defining the EOS relationship:
 *                  \f$P = k \rho^\gamma\f$
 *
 * @return Vector of size 2 containing:
 *         \parblock
 *         - \p [0]: Derivative of logarithmic mass (\f$\frac{d(\log m)}{d(\log r)}\f$)
 *         - \p [1]: Derivative of logarithmic pressure (\f$\frac{d(\log P)}{d(\log r)}\f$)
 *         \endparblock
 *
 * @pre \p m must be greater than zero (\f$m > 0\f$)
 * @pre Metric factor must be positive: \f$1 - \frac{2 G m}{r c^2} > 0\f$
 * @pre The EOS must be valid for the given constants \p k and \p gamma
 *
 * @note When \p m approaches zero, it is set to \f$10^{-30}\f$ to avoid division by zero
 * @note The TOV equations include relativistic corrections significant for neutron stars
 *
 * @warning May encounter floating-point errors for invalid state values
 * @warning Fails if Schwarzschild radius condition is violated (\f$2 G m / r c^2 \geq 1\f$)
 *
 * Example usage:
 * @code{.cpp}
 * double log_r = 0.0; // Logarithmic radius
 * std::vector<double> state = {1.0, 15.0}; // log10(m), log10(P)
 * double k = 1.2293e15; // EOS constant for relativistic neutron gas
 * double gamma = 4.0 / 3.0; // Polytropic index
 *
 * auto derivatives = tolman_oppenheimer_volkoff_derivatives(log_r, state, k, gamma);
 * std::cout << "dlogm_dlogr: " << derivatives[0] << std::endl;
 * std::cout << "dlogP_dlogr: " << derivatives[1] << std::endl;
 * @endcode
 *
 * @see newtonian() for the non-relativistic version
 */
[[nodiscard]] std::vector<double>
tolman_oppenheimer_volkoff_derivatives(double log_r, const std::vector<double> &state, double k,
                                       double gamma);

/**
 * @brief Computes the mass and pressure derivatives using TOV equations with spline-based EOS.
 *
 * @details
 * This function calculates the derivatives of logarithmic mass (\f$\log m\f$) and logarithmic
 * pressure
 * (\f$\log P\f$) with respect to the logarithmic radius (\f$\log r\f$) in the relativistic
 * framework, using a tabulated equation of state via GSL spline interpolation. This extends the TOV
 * solver to support realistic, complex equations of state beyond simple polytropic relationships.
 *
 * The TOV equations are identical to the polytropic version:
 * \f{eqnarray*}{
 * \frac{d(\log m)}{d(\log r)} &=& \frac{4 \pi r^3 \rho}{m} \\
 * \frac{d(\log P)}{d(\log r)} &=& \left( - \frac{G m \rho}{P r} \right)
 *                                  \cdot \left( 1 + \frac{P}{\rho c^2} \right)
 *                                  \cdot \left( 1 + \frac{4 \pi P r^3}{m c^2} \right)
 *                                  \cdot \left( 1 - \frac{2 G m}{r c^2} \right)^{-1}
 * \f}
 *
 * **Key Difference from Polytropic Version:**
 * Instead of using \f$\rho = (P/k)^{1/\gamma}\f$, this function obtains density from
 * the tabulated EOS via inverse spline interpolation: \f$\log\rho = f^{-1}(\log P)\f$
 *
 * **EOS Handling:**
 * - Pressure clamping: Ensures \f$P\f$ stays within EOS validity range [min_logP, max_logP]
 * - Spline evaluation: Uses GSL spline to get \f$\rho(P)\f$ from tabulated data
 * - Error handling: Graceful handling of extrapolation beyond EOS range
 *
 * The relativistic corrections appear as four factors:
 * \f{description}
 * \item[Factor 1] Newtonian term: $-\frac{G m \rho}{P r}$
 * \item[Factor 2] Internal energy: $1 + \frac{P}{\rho c^2}$
 * \item[Factor 3] Pressure effect: $1 + \frac{4 \pi P r^3}{m c^2}$
 * \item[Factor 4] Metric factor: $\left(1 - \frac{2 G m}{r c^2}\right)^{-1}$
 * \f}
 *
 * @param[in] log_r Logarithmic radius (\f$\log_{10}(r)\f$) in centimeters
 * @param[in] state Vector containing:
 *                  \parblock
 *                  - \p state[0]: Logarithmic mass (\f$\log_{10}(m)\f$) in grams
 *                  - \p state[1]: Logarithmic pressure (\f$\log_{10}(P)\f$) in dyne/cm²
 *                  \endparblock
 * @param[in] spline_inv GSL spline object for inverse EOS: \f$\log\rho = f^{-1}(\log P)\f$
 * @param[in] acc_inv GSL interpolation accelerator for inverse spline (for efficiency)
 * @param[in] min_logP Minimum log₁₀(pressure) for EOS validity range (dyne/cm²)
 * @param[in] max_logP Maximum log₁₀(pressure) for EOS validity range (dyne/cm²)
 *
 * @return Vector of size 2 containing:
 *         \parblock
 *         - \p [0]: Derivative of logarithmic mass (\f$\frac{d(\log m)}{d(\log r)}\f$)
 *         - \p [1]: Derivative of logarithmic pressure (\f$\frac{d(\log P)}{d(\log r)}\f$)
 *         \endparblock
 *
 * @pre \p m must be greater than zero (\f$m > 0\f$)
 * @pre Metric factor must be positive: \f$1 - \frac{2 G m}{r c^2} > 0\f$
 * @pre spline_inv != nullptr (valid GSL spline object)
 * @pre acc_inv != nullptr (valid GSL accelerator object)
 * @pre min_logP < max_logP (valid pressure range)
 *
 * @note When \p m approaches zero, it is set to \f$10^{-30}\f$ to avoid division by zero
 * @note Pressure values are clamped to EOS validity range before spline evaluation
 * @note The TOV equations include relativistic corrections significant for neutron stars
 *
 * @warning May encounter floating-point errors for invalid state values
 * @warning Fails if Schwarzschild radius condition is violated (\f$2 G m / r c^2 \geq 1\f$)
 * @warning GSL spline evaluation can fail for values outside interpolation range
 *
 * **Example Usage:**
 * @code{.cpp}
 * // Assuming splines are initialized
 * double log_r = 5.0; // log10(r) = 5 -> r = 100,000 cm
 * std::vector<double> state = {30.0, 25.0}; // log10(m), log10(P)
 *
 * auto derivatives = tolman_oppenheimer_volkoff_derivatives_spline(
 *     log_r, state, spline_inv, acc_inv, 15.0, 35.0
 * );
 *
 * std::cout << "dlogm_dlogr: " << derivatives[0] << std::endl;
 * std::cout << "dlogP_dlogr: " << derivatives[1] << std::endl;
 * @endcode
 *
 * @see tolman_oppenheimer_volkoff_derivatives() for polytropic EOS version
 * @see non_rotating_stellar_structure_spline() for full integration example
 * @see GSL documentation for spline usage and error handling
 */
[[nodiscard]] std::vector<double> tolman_oppenheimer_volkoff_derivatives_spline(
    double log_r, const std::vector<double> &state, const gsl_spline *spline_inv,
    gsl_interp_accel *acc_inv, double min_logP, double max_logP);

[[nodiscard]] std::vector<double> tolman_oppenheimer_volkoff_derivatives_spline_eps(
    double log_r, const std::vector<double> &state, const gsl_spline *rho_of_logP,
    gsl_interp_accel *acc_rho, const gsl_spline *eps_of_logP, gsl_interp_accel *acc_eps,
    double min_logP, double max_logP);

/**
 * @brief EOS parameters are now handled by the PolytropicEOS class.
 *
 * **Migration Notice:**
 * The `set_eos_parameters()` function has been removed and replaced with the
 * modular PolytropicEOS system. EOS parameters (k, γ, name) are now obtained
 * automatically within `non_rotating_stellar_structure()` using:
 *
 * ```cpp
 * PolytropicEOS eos_calculator;
 * auto eos_data = eos_calculator.getEOSParameters(eos_type);
 * double k = eos_data.k;
 * double gamma = eos_data.gamma;
 * std::string name = eos_data.name;
 * ```
 *
 * **Benefits of New Architecture:**
 * - Centralized EOS parameter management
 * - Consistent parameter values across the codebase
 * - Enhanced validation and error checking
 * - Support for custom mean molecular weights (μₑ)
 * - Extensible framework for additional EOS types
 *
 * **Parameter Values:**
 * The PolytropicEOS class provides the same parameter values that were
 * previously hardcoded:
 * - Electron non-relativistic: k = 1.0036×10¹³, γ = 5/3
 * - Electron relativistic: k = 1.2435×10¹⁵, γ = 4/3
 * - Neutron non-relativistic: k = 5.3802×10⁹, γ = 5/3
 * - Neutron relativistic: k = 1.2293×10¹⁵, γ = 4/3
 *
 * @see PolytropicEOS::getEOSParameters() for the new EOS parameter interface
 * @see non_rotating_stellar_structure() for updated function signature
 * @see PolytropicGasType for available EOS types
 */

/**
 * @brief Generates a descriptive filename for output data based on the EOS name and central
 * density.
 *
 * @details
 * This function constructs a filename for storing stellar structure data by combining
 * the EOS name and central density (\f$\rho_c\f$). The filename follows the format:
 *
 * \f{verbatim}
 * data/[eos_name]_rhoc_[formatted_density].csv
 * \f}
 *
 * **Density Formatting:**
 * The central density is formatted using scientific notation with specific transformations:
 * - Uses 2 decimal places precision (e.g., 1.23e+09)
 * - Replaces '+' with 'p' (e.g., 1.23e+09 → 1.23ep09)
 * - Replaces 'e' with 'p' (e.g., 1.23ep09 → 1.23pp09)
 * - Final result: 1.23pp09
 *
 * This encoding ensures filesystem-safe filenames while preserving density information.
 *
 * @param[in] name EOS identifier string (e.g., "electron_relativistic", "neutron_non_relativistic")
 * @param[in] rho_c Central density in g/cm³
 *
 * @return Formatted filename string in "data/" directory with ".csv" extension
 *
 * @note The function assumes the existence of a "data/" directory
 * @note No validation is performed on the input parameters
 * @note The encoding is reversible for data analysis purposes
 *
 * **Example Transformations:**
 * @code{.cpp}
 * // Example 1: Standard scientific notation
 * std::string name1 = "neutron_relativistic";
 * double rho_c1 = 1e14;
 * std::string file1 = get_filename(name1, rho_c1);
 * // Input: 1e14 → Scientific: "1.00e+14" → Replace '+': "1.00ep14" → Replace 'e': "1.00pp14"
 * // Result: "data/neutron_relativistic_rhoc_1.00pp14.csv"
 *
 * // Example 2: Fractional exponent
 * std::string name2 = "electron_non_relativistic";
 * double rho_c2 = 1.23e9;
 * std::string file2 = get_filename(name2, rho_c2);
 * // Input: 1.23e9 → Scientific: "1.23e+09" → Replace '+': "1.23ep09" → Replace 'e': "1.23pp09"
 * // Result: "data/electron_non_relativistic_rhoc_1.23pp09.csv"
 *
 * // Example 3: Negative exponent
 * std::string name3 = "electron_relativistic";
 * double rho_c3 = 1.5e-3;
 * std::string file3 = get_filename(name3, rho_c3);
 * // Input: 1.5e-3 → Scientific: "1.50e-03" → No '+' to replace → Replace 'e': "1.50p-03"
 * // Result: "data/electron_relativistic_rhoc_1.50p-03.csv"
 * @endcode
 *
 * **Filename Components:**
 * \f{description}
 * \item[Prefix] "data/" directory path
 * \item[Name] EOS identifier from PolytropicEOS system
 * \item[Separator] "\_rhoc\_" between name and density
 * \item[Density] Encoded scientific notation of central density
 * \item[Extension] ".csv" file extension
 * \f}
 *
 * @warning Directory "data/" must exist before using the generated filename
 * @warning Large density values may create very long filenames
 *
 * @see PolytropicEOS::getEOSParameters() for valid EOS names
 * @see non_rotating_stellar_structure() for usage context
 */
[[nodiscard]] std::string get_filename(const std::string &name, double rho_c);

#endif // NON_ROTATING_STELLAR_STRUCTURE_HPP
