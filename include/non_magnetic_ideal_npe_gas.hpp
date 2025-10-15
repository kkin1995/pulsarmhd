// non_magnetic_ideal_npe_gas.hpp
#ifndef NON_MAGNETIC_IDEAL_NPE_GAS_HPP
#define NON_MAGNETIC_IDEAL_NPE_GAS_HPP

/**
 * @file non_magnetic_ideal_npe_gas.hpp
 * @brief Implementation of equation of state calculations for non-magnetic neutron-proton-electron
 * gas
 * @author Karan Kinariwala
 * @date 2025-05-25
 *
 * @details This file contains the implementation of a non-magnetic neutron-proton-electron (NPE)
 * gas equation of state calculator, particularly relevant for neutron star physics. The
 * implementation considers three particle species (neutrons, protons, and electrons) and their
 * quantum statistical behavior in different density regimes.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>
#include <iomanip>
#include <iostream>
#include <memory> // for smart pointers
#include <string>
#include <vector>

/**
 * @class NonMagneticNPEGas
 * @brief A class for calculating the equation of state of non-magnetic neutron-proton-electron gas
 *
 * @details This class implements calculations for a three-component system consisting of neutrons,
 * protons, and electrons under beta equilibrium conditions. It handles three distinct density
 * regimes:
 * - Low density (electron-proton dominated)
 * - Transition region (neutron appearance threshold)
 * - High density (neutron-rich matter)
 *
 * The calculations include:
 * - Quantum statistical effects through Fermi-Dirac statistics
 * - Beta equilibrium conditions
 * - Relativistic corrections for high-density regimes
 *
 * Physical constants used in the calculations:
 * - Electron mass: 9.1094e-28 g
 * - Proton mass: 1.672622e-24 g
 * - Neutron mass: 1.674927e-24 g
 * - Reduced Planck constant: 1.0546e-27 erg·s
 * - Speed of light: 2.9979e10 cm/s
 *
 * @note All calculations are performed in CGS units
 *
 * Usage example:
 * @code{.cpp}
 * // Create an instance with debug mode enabled
 * NonMagneticNPEGas npe_gas(true);
 *
 * // Calculate EOS for a density range
 * npe_gas.calculateEOS("output.csv", 1.0e27, 1.0e35, 300);
 * @endcode
 *
 * @warning This implementation assumes ideal gas behavior and does not include
 * strong interactions between particles or magnetic field effects.
 *
 * @see BetaEquilibriumSolver for details on the beta equilibrium calculations
 */
class NonMagneticNPEGas {
private:
  /** @internal @brief Electon mass in grams */
  static constexpr double m_electron = 9.1094e-28;

  /** @internal @brief Proton mass in grams */
  static constexpr double m_proton = 1.672622e-24;

  /** @internal @brief Neutron mass in grams */
  static constexpr double m_neutron = 1.674927e-24;

  /** @internal @brief Reduced planck constant is erg*seconds */
  static constexpr double hbar = 1.0546e-27;

  /** @internal @brief Speed of light in cm/s */
  static constexpr double c = 2.9979e10;

  /**
   * @internal
   * @brief Neutron appearence threshold density in cm^{-3}
   * @details This represents the baryon number density at which neutrons
   * begin to appear.
   * */
  static constexpr double n_B_threshold = 7.37e29; // neutron appearence threshold

  /** @internal @brief Flag for enabling debug output */
  bool debug_mode;

  /** @internal @brief Electron Compton wavelength divided by 2\pi */
  double lambda_e;

  /** @internal @brief Proton Compton wavelength divided by 2\pi */
  double lambda_p;

  /** @internal @brief Neutron Compton wavelength divided by 2\pi */
  double lambda_n;

  /**
   * @internal
   * @brief Calculates the relativistic Fermi-Dirac integral for pressure calculations
   *
   * @details This function evaluates the dimensionless integral:
   * \f[
   * \phi(x) = \frac{1}{8\pi^2} \left[x\sqrt{1+x^2}(\frac{2}{3}x^2-1) + \ln(x+\sqrt{1+x^2})\right]
   * \f]
   * where x is the dimensionless Fermi momentum (p_F/mc).
   *
   * The function uses different approximations in three regimes:
   * 1. Non-relativistic (x < 0.1):
   *    \f[
   *    \phi(x) = \frac{1}{15\pi^2}(x^5 - \frac{5}{14}x^7 + \frac{5}{24}x^9)
   *    \f]
   * 2. Ultra-relativistic (x > 10):
   *    \f[
   *    \phi(x) = \frac{1}{12\pi^2}(x^4 - x^2 + \frac{3}{2}\ln(2x))
   *    \f]
   * 3. Intermediate regime (0.1 <= x <= 10):
   *    Uses the full expression
   *
   * @param x Dimensionless Fermi momentum (p_F/mc)
   * @return double The evaluated Fermi-Dirac integral for pressure
   * @note Used in calculating the degeneracy pressure of the fermionic components
   * @see chi() for the energy density counterpart
   */
  double phi(double x) const;

  /**
   * @internal
   * @brief Calculates the relativistic Fermi-Dirac integral for energy density calculations
   *
   * @details This function evaluates the dimensionless integral:
   * \f[
   * \chi(x) = \frac{1}{8\pi^2} \left[x\sqrt{1+x^2}(1 + 2x^2) - \ln(x+\sqrt{1+x^2})\right]
   * \f]
   * where x is the dimensionless Fermi momentum (p_F/mc).
   *
   * The function uses different approximations in three regimes:
   * 1. Non-relativistic (x < 0.1):
   *    \f[
   *    \chi(x) = \frac{1}{3\pi^2}(x^3 + \frac{3}{10}x^5 - \frac{3}{56}x^7)
   *    \f]
   * 2. Ultra-relativistic (x > 10):
   *    \f[
   *    \chi(x) = \frac{1}{4\pi^2}(x^4 + x^2 - \frac{1}{2}\ln(2x))
   *    \f]
   * 3. Intermediate regime (0.1 <= x <= 10):
   *    Uses the full expression
   *
   * @param x Dimensionless Fermi momentum (p_F/mc)
   * @return double The evaluated Fermi-Dirac integral for energy density
   * @note Used in calculating the energy density of the fermionic components
   * @see phi() for the pressure counterpart
   */
  double chi(double x) const;

  /**
   * @internal
   * @brief Performs cubic spline interpolation between data points
   *
   * @details Implements a cubic spline interpolation using the GSL (GNU Scientific Library)
   * for smooth transitions between different density regimes. This is particularly important
   * around the neutron appearance threshold where physical quantities need to transition
   * smoothly.
   *
   * The interpolation uses natural cubic splines with the following properties:
   * - Continuous up to second derivatives
   * - Natural boundary conditions (second derivatives are zero at endpoints)
   * - Monotonicity preservation in monotonic regions
   *
   * @param t The point at which to evaluate the interpolation
   * @param x Vector of x-coordinates of the data points (must be strictly increasing)
   * @param y Vector of y-coordinates of the data points
   * @return double The interpolated value at point t
   *
   * @throws NPEGasException if x and y vectors have different sizes
   * @throws NPEGasException if x vector is not strictly increasing
   * @throws NPEGasException if t is outside the range of x values
   *
   * @pre x and y must have the same size
   * @pre x must be strictly increasing
   * @pre t must be within the range [x[0], x[back]]
   *
   * @note This function uses GSL's cubic spline interpolation routines
   * @warning The function assumes responsibility for memory management of GSL resources
   */
  double interpolate_spline(double t, const std::vector<double> &x,
                            const std::vector<double> &y) const;

  /**
   * @class BetaEquilibriumSolver
   * @brief Solves for beta equilibrium conditions in neutron star matter
   * @private
   *
   * @details This private nested class implements a numerical solver for beta equilibrium
   * conditions in dense nuclear matter. Beta equilibrium represents the balance between the weak
   * interaction processes:
   * \f[
   * p + e^- \leftrightarrow n + \bar{\nu}_e
   * \f]
   *
   * In equilibrium, the chemical potentials satisfy:
   * \f[
   * \mu_p + \mu_e = \mu_n
   * \f]
   *
   * The solver uses GSL's root-finding algorithms to find the neutron fraction
   * that satisfies this equilibrium condition.
   */
  class BetaEquilibriumSolver {
  private:
    /**
     * @internal
     * @struct Parameters
     * @brief Container for physical parameters needed in equilibrium calculations
     *
     * @details Holds masses and dimensionless momenta required for calculating
     * chemical potentials of different particle species
     */
    struct Parameters {
      double m_e;
      double m_p;
      double m_n;
      double x_p;
      double x_e;
    };

    /**
     * @internal
     * @brief Equilibrium condition function to be solved
     *
     * @details Evaluates the function:
     * \f[
     * f(x_n) = m_n\sqrt{1 + x_n^2} - m_p\sqrt{1 + x_p^2} - m_e\sqrt{1 + x_e^2}
     * \f]
     * where x_n is the dimensionless neutron Fermi momentum
     *
     * @param x_n Dimensionless neutron Fermi momentum
     * @param p Pointer to Parameters struct
     * @return double Value of equilibrium function
     */
    static double f(double x_n, void *p);

    /**
     * @internal
     * @brief Derivative of the equilibrium condition function
     *
     * @details Evaluates:
     * \f[
     * f'(x_n) = -\frac{m_n x_n}{\sqrt{1 + x_n^2}}
     * \f]
     *
     * @param x_n Dimensionless neutron Fermi momentum
     * @param p Pointer to Parameters struct
     * @return double Derivative value
     */
    static double f_prime(double x_n, void *p);

    /**
     * @internal
     * @brief Combined evaluation of function and its derivative
     *
     * @details Evaluates both f(x_n) and f'(x_n) for efficiency in
     * Newton-type methods
     *
     * @param x_n Dimensionless neutron Fermi momentum
     * @param p Pointer to Parameters struct
     * @param[out] f_value Function value
     * @param[out] df_value Derivative value
     */
    static void f_and_df(double x_n, void *p, double *f_value, double *df_value);

  public:
    /**
     * @brief Solves the beta equilibrium condition
     *
     * @details Uses GSL's Brent algorithm to find the neutron fraction that
     * satisfies beta equilibrium. The solver handles three regimes:
     * 1. Initial guess refinement
     * 2. Root bracketing
     * 3. Iterative solution
     *
     * @param x_n Initial guess for dimensionless neutron Fermi momentum
     * @param n_B Baryon number density in cm^{-3}
     * @param m_e Electron mass in g
     * @param m_p Proton mass in g
     * @param m_n Neutron mass in g
     * @param x_p Dimensionless proton Fermi momentum
     * @param x_e Dimensionless electron Fermi momentum
     * @param debug Enable debug output
     *
     * @return double Solution for x_n that satisfies equilibrium
     *
     * @throws NPEGasException if root finding fails
     * @throws NPEGasException if maximum iterations exceeded
     *
     * @note Uses GSL's Brent algorithm for robustness
     * @warning May not converge for extremely high densities
     */
    static double solve(double x_n, double n_B, double m_e, double m_p, double m_n, double x_p,
                        double x_e, bool debug);
  };

public:
  /**
   * @brief Constructor for NonMagneticNPEGas
   *
   * @details Initializes a new instance of the NPE gas calculator with:
   * - Computation of Compton wavelengths for all particle species
   * - Optional debug mode for detailed output
   *
   * @param debug Enable debug output (default: false)
   */
  explicit NonMagneticNPEGas(bool debug = false);

  /**
   * @brief Calculates the equation of state for the NPE gas
   *
   * @details Computes the following quantities across a range of densities:
   * - Baryon number density
   * - Pressure
   * - Mass density
   *
   * The calculation handles three density regimes:
   * 1. Low density (n_B < n_B_threshold): e-p gas
   * 2. Transition region: appearance of neutrons
   * 3. High density (n_B > n_B_threshold): n-rich matter
   *
   * @param filename Output file path for EOS data
   * @param nB_min Minimum baryon number density in cm^{-3}
   * @param nB_max Maximum baryon number density in cm^{-3}
   * @param num_points Number of points in density array
   *
   * @throws NPEGasException if file operations fail
   * @throws NPEGasException if beta equilibrium solver fails
   *
   * @note Output file format: log_n,log_P,log_rho (CSV)
   * @warning Convergence not guaranteed for n_B > 1e36 cm^{-3}
   */
  void calculateEOS(const std::string &filename, double nB_min, double nB_max, int num_points);
};

#endif // NON_MAGNETIC_IDEAL_NPE_GAS_HPP
