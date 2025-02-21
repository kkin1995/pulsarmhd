// magnetic_bps.hpp
#ifndef MAGNETIC_BPS_HPP
#define MAGNETIC_BPS_HPP

#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

// Structure to store atomic mass data
struct AtomicMass {
    int Z;  // Atomic number
    int A;  // Mass number
    std::string element; // Element symbol
    double dM;  // Mass excess (keV)
    double BE;  // Binding energy (keV)
    double BDE; // Beta-decay energy (keV)
    double M;   // Atomic mass (u)
};

struct EquilibriumComposition {
    // Baryon Density and optimal nuclear properties
    double baryon_density;  // cm^-3
    int optimal_A;          // Optimal mass number
    int optimal_Z;          // Optimal atomic number
    std::string optimal_element; // Optimal element symbol

    // Thermodynamic quantities
    double gibbs_free_energy;   // Minimum Gibbs Free Energy (erg)
    double total_energy_density;    // Total energy density (erg cm^-3)
    double total_pressure;          // Total pressure (dyn cm^-2)
    double total_mass_density;      // Total mass density (g cm^-3)

    // Electron properties
    double electron_density;    // Electron density (cm^-3)
    double electron_energy_density; // Electron energy density (erg cm^-3)
    double electron_pressure;   // Electron pressure (dyn cm^-2)
    double gamma_e; // Dimensionless electron chemical potential
    double max_landau_level;  // Maximum Occupied Landau level for electrons

    // Lattice properties
    double lattice_energy_density;  // Lattice energy density (erg cm^-3)
    double lattice_pressure;    // Lattice pressure (dyn cm^-2)

    // Convergence metrics
    bool converged; // Whether bisection converged
    double relative_error; // Final relative error
    int iterations; // Number of iterations needed
};

class MagneticBPSEOS {
    public:
        // Constructor that sets up the EOS with parameters and loads atomic masses.
        MagneticBPSEOS(const std::string& atomic_mass_file, double B_ratio_electron, double rel_tol = 1e-6, double abs_tol = 1e-8);
        
        EquilibriumComposition computeEquilibriumComposition(double nB);

        void writeEOSResults(const std::string& output_file, const std::vector<EquilibriumComposition>& results) const;
        
        void printEquilibriumTable(const std::vector<EquilibriumComposition>& results) const;

        std::vector<EquilibriumComposition> runSimulation(double nB_min, double nB_max, int num_points);

    private:
        // Configuration parameters
        double B_ratio_electron_;
        double rel_tolerance_;
        double abs_tolerance_;

        // Atomic mass data
        std::vector<AtomicMass> atomic_masses_;

        // Physical constants (using CGS units)
        static constexpr double hbar = 1.0546e-27;
        static constexpr double c = 2.9979e10;
        static constexpr double m_electron = 9.1094e-28;
        static constexpr double e_charge = 4.80320471e-10; // CGS units
        static constexpr double lambda_e = hbar / (m_electron * c); // Electron Compton wavelength

        static double psi(double x);
        static double eta(double x);

        std::vector<AtomicMass> readAtomicMasses(const std::string& filename);
        double getAtomicMass(int A, int Z) const;
        std::vector<double> generateBaryonDensities(double nB_min, double nB_max, int num_points);
};

#endif // MAGNETIC_BPS_HPP