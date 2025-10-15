#ifndef MAGNETIC_BPS_HPP
#define MAGNETIC_BPS_HPP

/**
 * @file magnetic_bps.hpp
 * @brief Implementation of the magnetic BPS (Baym-Pethick-Sutherland) equation of state
 *        for strongly magnetized neutron star crusts.
 * @author Karan Amit Kinariwala
 * @date 2025-05-25
 *
 * This class implements the BPS equation of state modified for strong magnetic fields,
 * considering the effects of Landau quantization on electrons. The implementation
 * calculates the ground state composition of matter in neutron star crusts under
 * strong magnetic fields by minimizing the Gibbs free energy.
 *
 * The implementation follows the formalism developed in:
 * Lai, D., & Shapiro, S. L. (1991). Cold equation of state in a strong magnetic field -
 * Effects of inverse beta-decay. The Astrophysical Journal, 383, 745.
 * DOI: 10.1086/170831
 *
 * Physical regime:
 * - Density range: typically 10^6 - 10^11 g/cm³
 * - Magnetic field: B > B_crit (where B_crit ≈ 4.414e13 G)
 * - Temperature: T = 0 (cold catalyzed matter)
 *
 * Units: CGS unless otherwise specified
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/**
 * @struct AtomicMass
 * @brief Structure to store nuclear mass and properties data
 *
 * Stores comprehensive nuclear data for each isotope, including mass excess,
 * binding energy, and beta-decay energy.
 */
struct AtomicMass {
  int Z;               ///< Atomic number
  int A;               ///< Mass number
  std::string element; ///< Element symbol
  double dM;           ///< Mass excess (keV)
  double BE;           ///< Binding energy (keV)
  double BDE;          ///< Beta-decay energy (keV)
  double M;            ///< Atomic mass (u)
};

/**
 * @struct EquilibriumComposition
 * @brief Structure containing the equilibrium state properties
 *
 * Stores the complete thermodynamic state and composition of matter
 * at a given density, including electron properties in quantizing magnetic fields
 * and lattice contributions.
 */
struct EquilibriumComposition {
  // Baryon Density and optimal nuclear properties
  double baryon_density;       ///< cm^-3
  int optimal_A;               ///< Optimal mass number
  int optimal_Z;               ///< Optimal atomic number
  std::string optimal_element; ///< Optimal element symbol

  // Thermodynamic quantities
  double gibbs_free_energy;    ///< Minimum Gibbs Free Energy (erg)
  double total_energy_density; ///< Total energy density (erg cm^-3)
  double total_pressure;       ///< Total pressure (dyn cm^-2)
  double total_mass_density;   ///< Total mass density (g cm^-3)

  // Electron properties
  double electron_density;        ///< Electron density (cm^-3)
  double electron_energy_density; ///< Electron energy density (erg cm^-3)
  double electron_pressure;       ///< Electron pressure (dyn cm^-2)
  double gamma_e;                 ///< Dimensionless electron chemical potential
  double max_landau_level;        ///< Maximum Occupied Landau level for electrons

  // Lattice properties
  double lattice_energy_density; ///< Lattice energy density (erg cm^-3)
  double lattice_pressure;       ///< Lattice pressure (dyn cm^-2)

  // Convergence metrics
  bool converged;        ///< Whether bisection converged
  double relative_error; ///< Final relative error
  int iterations;        ///< Number of iterations needed
};

/**
 * @class MagneticBPSEOS
 * @brief Magnetic BPS equation of state solver
 *
 * Computes the equation of state and ground state composition of neutron star
 * crustal matter in strong magnetic fields, following Lai & Shapiro (1991).
 *
 * Key equations implemented:
 *
 * 1. Electron energy density in quantizing magnetic field:
 * \f[
 * \varepsilon_e = \frac{2B}{(2\pi)^2\lambda_e^3} m_e c^2
 * \sum_{\nu=0}^{\nu_{max}} g_\nu (1 + 2\nu B/B_c)
 * \psi\left(\frac{x_e(\nu)}{\sqrt{1 + 2\nu B/B_c}}\right)
 * \f]
 * where \f$\nu_{max}\f$ is the maximum Landau level, \f$B_c\f$ is the critical field,
 * and \f$\psi(x)\f$ is defined below.
 *
 * 2. Electron pressure:
 * \f[
 * P_e = \frac{2B}{(2\pi)^2\lambda_e^3} m_e c^2
 * \sum_{\nu=0}^{\nu_{max}} g_\nu (1 + 2\nu B/B_c)
 * \eta\left(\frac{x_e(\nu)}{\sqrt{1 + 2\nu B/B_c}}\right)
 * \f]
 *
 * 3. Lattice energy:
 * \f[
 * \varepsilon_L = -1.444 Z^{2/3} e^2 n_e^{4/3}
 * \f]
 *
 * 4. Gibbs free energy per nucleon:
 * \f[
 * g = \frac{M(A,Z)}{A} + \frac{Z}{A}(\gamma_e m_e c^2 - m_e c^2)
 * + \frac{4}{3}\frac{Z}{A}\frac{\varepsilon_L}{n_e}
 * \f]
 *
 * The auxiliary functions \f$\psi(x)\f$ and \f$\eta(x)\f$ are defined as:
 * \f[
 * \psi(x) = \frac{1}{2}x\sqrt{1+x^2} + \frac{1}{2}\ln(x + \sqrt{1+x^2})
 * \f]
 * \f[
 * \eta(x) = \frac{1}{2}x\sqrt{1+x^2} - \frac{1}{2}\ln(x + \sqrt{1+x^2})
 * \f]
 *
 * Physical regime:
 * - Density range: typically \f$10^6 - 10^{11}\f$ g/cm³
 * - Magnetic field: \f$B > B_c\f$ where \f$B_c \approx 4.414 \times 10^{13}\f$ G
 * - Temperature: \f$T = 0\f$ (cold catalyzed matter)
 *
 * The equilibrium composition at each density is determined by minimizing the
 * Gibbs free energy with respect to A and Z, considering the effects of
 * Landau quantization on the electron phase space.
 */
class MagneticBPSEOS {
public:
  /**
   * @brief Constructs a magnetic BPS EOS solver
   * @param atomic_mass_file Path to the nuclear mass data file
   * @param B_ratio_electron Magnetic field strength relative to critical field (B/B_c)
   * @param rel_tol Relative tolerance for convergence (default: 1e-6)
   * @param abs_tol Absolute tolerance for convergence (default: 1e-8)
   *
   * @note The critical field \f$B_c = m_e^2c^3/(e\hbar)\f$ \f$\approx 4.414 \times 10^{13}\f$ G
   * @throws std::runtime_error if atomic mass file cannot be loaded
   */
  MagneticBPSEOS(const std::string &atomic_mass_file, double B_ratio_electron,
                 double rel_tol = 1e-6, double abs_tol = 1e-8);

  /**
   * @brief Computes the equilibrium composition at a given baryon density
   * @param nB Baryon number density in cm\f$^{-3}\f$
   * @return EquilibriumComposition containing the ground state properties
   *
   * @details For a given baryon density, this method:
   * 1. Loops through all possible nuclei in the atomic mass table
   * 2. For each nucleus, solves the electron chemical potential using bisection
   * 3. Computes the Gibbs free energy per nucleon:
   * \f[
   * g = \frac{M(A,Z)}{A} + \frac{Z}{A}(\gamma_e m_e c^2 - m_e c^2)
   * + \frac{4}{3}\frac{Z}{A}\frac{\varepsilon_L}{n_e}
   * \f]
   * 4. Selects the nucleus that minimizes the Gibbs free energy
   *
   * @note The method ensures convergence of the electron chemical potential
   *       within the specified tolerances (rel_tolerance_ and abs_tolerance_)
   */
  EquilibriumComposition computeEquilibriumComposition(double nB);

  /**
   * @brief Writes equation of state results to a file in CSV format
   * @param output_file Path to the output file
   * @param results Vector of EquilibriumComposition containing results at different densities
   *
   * @details Writes the following quantities in log10 scale:
   * - Baryon number density (cm\f$^{-3}\f$)
   * - Mass density (g cm\f$^{-3}\f$)
   * - Pressure (dyn cm\f$^{-2}\f$)
   *
   * Output format:
   * log_n,log_rho,log_P
   *
   * @throws Outputs error message to cerr if file cannot be opened
   */
  void writeEOSResults(const std::string &output_file,
                       const std::vector<EquilibriumComposition> &results) const;

  /**
   * @brief Prints a formatted table of equilibrium compositions
   * @param results Vector of EquilibriumComposition containing results at different densities
   *
   * @details Displays a table with the following columns:
   * - Baryon density (cm\f$^{-3}\f$)
   * - Mass number (A)
   * - Atomic number (Z)
   * - Element symbol
   * - Gibbs free energy (erg)
   * - Total pressure (dyn cm\f$^{-2}\f$)
   * - Mass density (g cm\f$^{-3}\f$)
   * - Maximum Landau level
   *
   * @note Results that did not achieve convergence are marked with an asterisk (*)
   */
  void printEquilibriumTable(const std::vector<EquilibriumComposition> &results) const;

  /**
   * @brief Runs the complete EOS simulation over a range of densities
   * @param nB_min Minimum baryon density in cm\f$^{-3}\f$
   * @param nB_max Maximum baryon density in cm\f$^{-3}\f$
   * @param num_points Number of density points to compute
   * @return Vector of EquilibriumComposition containing results at each density
   *
   * @details This method:
   * 1. Generates logarithmically spaced density points
   * 2. Computes equilibrium composition at each density
   * 3. Returns the complete equation of state
   *
   * @note This is the main method to generate the complete equation of state.
   *       The density range should typically be within the BPS regime
   *       (\f$10^6 - 10^{11}\f$ g cm\f$^{-3}\f$)
   */
  std::vector<EquilibriumComposition> runSimulation(double nB_min, double nB_max, int num_points);

private:
  // Configuration parameters
  double B_ratio_electron_; ///< Ratio of magnetic field to critical field (B/B_c)
  double rel_tolerance_;    ///< Relative tolerance for convergence
  double abs_tolerance_;    ///< Absolute tolerance for convergence

  /// @brief Vector storing nuclear mass and property data
  /// @details Contains information about different nuclei including mass numbers,
  ///          atomic numbers, binding energies, and atomic masses loaded from the
  ///          input file during initialization
  std::vector<AtomicMass> atomic_masses_;

  // Physical constants in CGS units
  static constexpr double hbar = 1.0546e-27;                  ///< Reduced Planck constant (erg s)
  static constexpr double c = 2.9979e10;                      ///< Speed of light (cm/s)
  static constexpr double m_electron = 9.1094e-28;            ///< Electron mass (g)
  static constexpr double e_charge = 4.80320471e-10;          ///< Elementary charge (statcoulomb)
  static constexpr double lambda_e = hbar / (m_electron * c); ///< Electron Compton wavelength (cm)

  /**
   * @brief Computes the electron energy function psi(x)
   * @param x Dimensionless momentum
   * @return \f$\psi(x) = \frac{1}{2}x\sqrt{1+x^2} + \frac{1}{2}\ln(x + \sqrt{1+x^2})\f$
   */
  static double psi(double x);

  /**
   * @brief Computes the electron pressure function eta(x)
   * @param x Dimensionless momentum
   * @return \f$\eta(x) = \frac{1}{2}x\sqrt{1+x^2} - \frac{1}{2}\ln(x + \sqrt{1+x^2})\f$
   */
  static double eta(double x);

  /**
   * @brief Reads atomic mass data from a CSV file
   * @param filename Path to the CSV file containing nuclear mass data
   * @return Vector of AtomicMass structures containing nuclear data
   *
   * @details The CSV file should have the following format:
   * Z,A,Element,dM,BE,BDE,M
   * where:
   * - Z: Atomic number
   * - A: Mass number
   * - Element: Element symbol
   * - dM: Mass excess in keV
   * - BE: Binding energy in keV
   * - BDE: Beta-decay energy in keV
   * - M: Atomic mass in atomic mass units (u)
   *
   * @note The first line of the file is assumed to be a header and is skipped
   *
   * @throws Outputs error message to cerr if file cannot be opened
   * @return Empty vector if file reading fails, populated vector otherwise
   */
  std::vector<AtomicMass> readAtomicMasses(const std::string &filename);

  /**
   * @brief Retrieves atomic mass for a given nucleus
   * @param A Mass number of the nucleus
   * @param Z Atomic number of the nucleus
   * @return Atomic mass in atomic mass units (u)
   *
   * @details Searches the loaded atomic mass data to find the matching nucleus.
   * If no matching nucleus is found, outputs an error message and returns -1.0
   *
   * @note This method assumes that atomic mass data has been successfully loaded
   *       during initialization
   *
   * @throws Outputs error message to cerr if nucleus is not found in the database
   * @return -1.0 if nucleus not found, actual mass otherwise
   */
  double getAtomicMass(int A, int Z) const;

  /**
   * @brief Generates a logarithmically spaced array of baryon densities
   * @param nB_min Minimum baryon density in cm\f$^{-3}\f$
   * @param nB_max Maximum baryon density in cm\f$^{-3}\f$
   * @param num_points Number of density points to generate
   * @return Vector of baryon densities
   *
   * @details Generates density points using the formula:
   * \f[
   * n_B(i) = n_{B,min} \times 10^{i \cdot \frac{\log_{10}(n_{B,max}/n_{B,min})}{num\_points-1}}
   * \f]
   * where i ranges from 0 to num_points-1
   *
   * @note The spacing is logarithmic to better sample the equation of state
   *       across multiple orders of magnitude in density
   */
  std::vector<double> generateBaryonDensities(double nB_min, double nB_max, int num_points);
};

#endif // MAGNETIC_BPS_HPP
