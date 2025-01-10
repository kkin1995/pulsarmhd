// non_rotating_stellar_structure.hpp
#ifndef NON_ROTATING_STELLAR_STRUCTURE_HPP
#define NON_ROTATING_STELLAR_STRUCTURE_HPP

/**
 * @file non_rotating_stellar_structure.hpp
 * @brief Header file for non-rotating stellar structure simulations.
 * @author Karan Kinariwala
 * @date 2025-01-09
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
 * @li DegenerateGasType: Enum representing the type of degenerate gas
 *
 * @section functions Functions
 * @li newtonian: Computes derivatives for Newtonian hydrostatic equilibrium
 * @li tolman_oppenheimer_volkoff_derivatives: Computes derivatives for relativistic stellar structure
 * @li set_eos_parameters: Sets the EOS parameters (k and gamma) based on the gas type
 * @li get_filename: Generates filenames for output data based on EOS and central density
 *
 * @section example Example Usage
 * @code{.cpp}
 * std::string eos_name;
 * double k, gamma;
 * set_eos_parameters(ELECTRON_RELATIVISTIC, eos_name, k, gamma);
 * 
 * std::vector<double> state = {log_m, log_p};
 * std::vector<double> derivatives = tolman_oppenheimer_volkoff_derivatives(log_r, state, k, gamma);
 * @endcode
 *
 * @see rk4.hpp for numerical integration using Runge-Kutta methods.
 */


#include <string>
#include <vector>

// Physical Constants in CGS

/**
 * @brief Gravitational constant in CGS units.
 */
const double G = 6.67430e-8; // gravitational constant in CGS

/**
 * @brief Speed of light in CGS units.
 */
const double c = 2.99792458e10;  // speed of light in CGS units

/**
 * @brief Types of degenerate gases used in stellar structure simulations.
 * 
 * This enumeration defines different types of degenerate gases that can be used
 * in stellar structure calculations. The distinction between relativistic and
 * non-relativistic cases affects the equation of state and subsequent stellar
 * structure calculations.
 *
 * @var DegenerateGasType::ELECTRON_NON_RELATIVISTIC
 * Non-relativistic electron gas. Used for white dwarfs where electron velocities
 * are much smaller than the speed of light (v << c).
 *
 * @var DegenerateGasType::ELECTRON_RELATIVISTIC
 * Relativistic electron gas. Applicable for more massive white dwarfs where
 * electron velocities approach the speed of light.
 *
 * @var DegenerateGasType::NEUTRON_NON_RELATIVISTIC
 * Non-relativistic neutron gas. Used in low-density regions of neutron stars
 * where neutron velocities are much smaller than the speed of light.
 *
 * @var DegenerateGasType::NEUTRON_RELATIVISTIC
 * Relativistic neutron gas. Applicable for high-density regions in neutron stars
 * where neutron velocities approach the speed of light.
 */
enum DegenerateGasType {
    ELECTRON_NON_RELATIVISTIC,  ///< Non-relativistic electron degenerate gas
    ELECTRON_RELATIVISTIC,      ///< Relativistic electron degenerate gas
    NEUTRON_NON_RELATIVISTIC,   ///< Non-relativistic neutron degenerate gas
    NEUTRON_RELATIVISTIC        ///< Relativistic neutron degenerate gas
};
// Function declarations

/**
 * @brief Computes the mass and pressure derivatives using the Newtonian equations of stellar structure.
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
std::vector<double> newtonian(double log_r, const std::vector<double>& state, double k, double gamma);

/**
 * @brief Computes the mass and pressure derivatives using the TOV (Tolman-Oppenheimer-Volkoff) equations.
 *
 * @details
 * This function calculates the derivatives of logarithmic mass (\f$\log m\f$) and logarithmic pressure 
 * (\f$\log P\f$) with respect to the logarithmic radius (\f$\log r\f$) in the relativistic framework.
 * The TOV equations account for general relativistic corrections, providing a more accurate 
 * description of compact stellar objects like neutron stars.
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
std::vector<double> tolman_oppenheimer_volkoff_derivatives(double log_r, const std::vector<double>& state, double k, double gamma);

/**
 * @brief Sets the parameters for the equation of state (EOS) based on the type of degenerate gas.
 *
 * This function assigns values to the EOS parameters \p k and \p gamma, as well as a descriptive name 
 * for the EOS, based on the selected type of degenerate gas. The parameters describe the relationship 
 * between pressure and density for the gas: \f$ P = k \rho^\gamma \f$.
 *
 * ### Supported Gas Types:
 * - \p ELECTRON_NON_RELATIVISTIC: Non-relativistic electron gas.
 * - \p ELECTRON_RELATIVISTIC: Relativistic electron gas.
 * - \p NEUTRON_NON_RELATIVISTIC: Non-relativistic neutron gas.
 * - \p NEUTRON_RELATIVISTIC: Relativistic neutron gas.
 *
 * ### Parameters:
 * @param[in] type Enum representing the type of degenerate gas.
 *                 - \p ELECTRON_NON_RELATIVISTIC
 *                 - \p ELECTRON_RELATIVISTIC
 *                 - \p NEUTRON_NON_RELATIVISTIC
 *                 - \p NEUTRON_RELATIVISTIC
 * @param[out] name A reference to a string where the EOS name will be stored. 
 *                  Example values include "electron_non_relativistic" and "neutron_relativistic".
 * @param[out] k A reference to the EOS proportionality constant (\f$k\f$) describing the relationship 
 *               between pressure and density.
 * @param[out] gamma A reference to the polytropic index (\f$\gamma\f$) defining the EOS power-law relationship.
 *
 * ### Constants:
 * - \p mu_e = 2.0: Mean molecular weight per electron, assumed constant.
 *
 * ### Example Usage:
 * @code
 * std::string eos_name;
 * double k, gamma;
 * DegenerateGasType type = NEUTRON_RELATIVISTIC;
 * 
 * set_eos_parameters(type, eos_name, k, gamma);
 * 
 * std::cout << "EOS Name: " << eos_name << std::endl;
 * std::cout << "k: " << k << std::endl;
 * std::cout << "gamma: " << gamma << std::endl;
 * @endcode
 *
 * ### Notes:
 * - The EOS constants for each gas type are pre-defined based on standard physical models.
 * - The parameter \p mu_e is used only for electron gases.
 *
 * @see DegenerateGasType for available gas types.
 */
void set_eos_parameters(DegenerateGasType type, std::string &name, double &k, double &gamma);

/**
 * @brief Generates a descriptive filename for output data based on the EOS name and central density.
 *
 * @details
 * This function constructs a filename for storing stellar structure data by combining
 * the EOS name and central density (\f$\rho_c\f$). The filename follows the format:
 *
 * \f{verbatim}
 * data/[eos_name]_rhoc_[formatted_density].csv
 * \f}
 *
 * The central density is formatted with scientific notation, where:
 * - 'e' is replaced with 'p' (e.g., 1e9 → 1p9)
 * - '+' is replaced with 'p' (e.g., 1e+09 → 1p09)
 * - Two decimal places are maintained (e.g., 1.00p9)
 *
 * @param[in] name EOS identifier string (e.g., "electron_relativistic")
 * @param[in] rho_c Central density in g/cm³
 *
 * @return Formatted filename string following the pattern above
 *
 * @note The function assumes the existence of a "data/" directory
 * @note No validation is performed on the input parameters
 *
 * Example with different central densities:
 * @code{.cpp}
 * // Example 1: Standard notation
 * std::string name1 = "neutron_non_relativistic";
 * double rho_c1 = 1e14;
 * std::string file1 = get_filename(name1, rho_c1);
 * // Result: "data/neutron_non_relativistic_rhoc_1p00p14.csv"
 *
 * // Example 2: Different precision
 * std::string name2 = "electron_relativistic";
 * double rho_c2 = 1.23e9;
 * std::string file2 = get_filename(name2, rho_c2);
 * // Result: "data/electron_relativistic_rhoc_1p23p9.csv"
 * @endcode
 *
 * @par Filename Components
 * \f{description}
 * \item[Prefix] "data/" directory path
 * \item[Name] EOS identifier from input parameter
 * \item[Separator] "\_rhoc\_" between name and density
 * \item[Density] Formatted central density value
 * \item[Extension] ".csv" file extension
 * \f}
 *
 * @warning Directory "data/" must exist before using the generated filename
 * @see set_eos_parameters() for generating valid EOS names
 */
std::string get_filename(const std::string& name, double rho_c);

#endif // NON_ROTATING_STELLAR_STRUCTURE_HPP