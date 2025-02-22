/**
* Implementation of the magnetic BPS equation of state.
* Based on: Lai & Shapiro (1991), ApJ, 383, 745
*/

#include "magnetic_bps.hpp"
#include <vector>
#include <cmath>
#include <iostream>

// Constructor: Initializes configuration parameters and loads atomic mass data.
MagneticBPSEOS::MagneticBPSEOS(const std::string& atomic_mass_file, double B_ratio_electron, double rel_tol, double abs_tol) 
    : B_ratio_electron_(B_ratio_electron), rel_tolerance_(rel_tol),
    abs_tolerance_(abs_tol),
    atomic_masses_(readAtomicMasses(atomic_mass_file))
    {
        if (atomic_masses_.empty()) {
        std::cerr << "Error: Atomic mass data could not be loaded from " 
            << atomic_mass_file << std::endl;
        }
    }

std::vector<double> MagneticBPSEOS::generateBaryonDensities(double nB_min, double nB_max, int num_points) {
    std::vector<double> nB_list;    
    // Logarithmic spacing for smoother transitions
    for (int i = 0; i < num_points; i++) {
        double nB_i = nB_min * pow(10, i * (log10(nB_max / nB_min) / (num_points - 1)));
        nB_list.push_back(nB_i);
    }
    return nB_list;
}

std::vector<AtomicMass> MagneticBPSEOS::readAtomicMasses(const std::string& filename) {
    std::vector<AtomicMass> atomicMasses;
    std::ifstream file(filename);
    std::string line;
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return atomicMasses;
    }
    
    // Read the header line (skip it)
    std::getline(file, line);
    
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        AtomicMass entry;
        std::string temp;
        
        std::getline(ss, temp, ','); entry.Z = std::stoi(temp);
        std::getline(ss, temp, ','); entry.A = std::stoi(temp);
        std::getline(ss, entry.element, ',');
        std::getline(ss, temp, ','); entry.dM = std::stod(temp);
        std::getline(ss, temp, ','); entry.BE = std::stod(temp);
        std::getline(ss, temp, ','); entry.BDE = std::stod(temp);
        std::getline(ss, temp, ','); entry.M = std::stod(temp);
        
        atomicMasses.push_back(entry);
    }
    
    file.close();
    return atomicMasses;
}

// Function to find M(A, Z) for given A and Z
double MagneticBPSEOS::getAtomicMass(int A, int Z) const {
    for (const auto& entry : atomic_masses_) {
        if (entry.A == A && entry.Z == Z) {
            return entry.M;
        }
    }
    std::cerr << "Error: No data found for A=" << A << " Z=" << Z << std::endl;
    return -1.0; // Return an error value if not found
}

double MagneticBPSEOS::psi(double x) {
    if (fabs(x) < 1e-12) return 0.0;

    double first_term = (1.0 / 2.0) * x * sqrt(1.0 + (x * x));
    double second_term = (1.0 / 2.0) * std::log(x + sqrt(1.0 + (x * x)));

    return first_term + second_term;
}

double MagneticBPSEOS::eta(double x) {
    if (fabs(x) < 1e-12) return 0.0;

    double first_term = (1.0 / 2.0) * x * sqrt(1.0 + (x * x));
    double second_term = (1.0 / 2.0) * std::log(x + sqrt(1.0 + (x * x)));

    return first_term - second_term;
}

void MagneticBPSEOS::printEquilibriumTable(const std::vector<EquilibriumComposition>& results) const {
    std::cout << "\n=== Equilibrium Composition Results ===\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << std::setw(12) << "Density" 
              << std::setw(8) << "A" 
              << std::setw(8) << "Z" 
              << std::setw(8) << "Element"
              << std::setw(15) << "Gibbs Energy"
              << std::setw(15) << "Pressure" 
              << std::setw(15) << "Mass Density" 
              << std::setw(12) << "ν_max" << std::endl;
              
    for (const auto& comp : results) {
        std::cout << std::setw(12) << comp.baryon_density
                  << std::setw(8) << comp.optimal_A
                  << std::setw(8) << comp.optimal_Z
                  << std::setw(8) << comp.optimal_element
                  << std::setw(15) << comp.gibbs_free_energy
                  << std::setw(15) << comp.total_pressure
                  << std::setw(15) << comp.total_mass_density;
        // Print Landau level separately to ensure proper spacing
        std::cout << std::setw(12) << comp.max_landau_level
                  << (comp.converged ? "" : " *") << std::endl;
    }
}

void MagneticBPSEOS::writeEOSResults(const std::string& output_file, const std::vector<EquilibriumComposition>& results) const {
    std::ofstream outfile(output_file);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << output_file << std::endl;
        return;
    }
    outfile << "log_n,log_rho,log_P\n";

    for (const auto& comp : results) {
        outfile << std::log10(comp.baryon_density) << "," << std::log10(comp.total_mass_density) << "," << std::log10(comp.total_pressure) << "\n";
    }
    outfile.close();
}

EquilibriumComposition MagneticBPSEOS::computeEquilibriumComposition(double nB) {
    EquilibriumComposition best_composition;
    best_composition.baryon_density = nB;
    best_composition.gibbs_free_energy = std::numeric_limits<double>::max();
    best_composition.optimal_A = 0;
    best_composition.optimal_Z = 0;
    best_composition.optimal_element = "None";
    best_composition.converged = false;

    double total_pressure = 0.0;
    double total_mass_density = 0.0;

    for (const auto& nucleus : atomic_masses_) {
        int A = nucleus.A;
        int Z = nucleus.Z;
        double mass = nucleus.M;

        if (mass <= 0) {
            continue;
        }

        double mass_energy = mass * 931.494 * 1.602e-6; // Convert u to MeV, then to erg

        double calculated_electron_density;
        double calculated_electron_energy_density;
        double calculated_electron_pressure;

        double nucleon_density = nB / A;
        double electron_density = (Z * nB) / A;

        // Solve Coupled Equations
        double x_e = pow(3.0 * pow(M_PI, 2.0) * pow(lambda_e, 3.0) * electron_density, 1.0 / 3.0);

        double gamma_e_lower = 1.0 + 1e-10;
        double gamma_e_upper = std::max(gamma_e_lower, sqrt(1.0 + x_e * x_e));

        double nu_m = 0.0;
        double relative_error = 0.0;
        bool converged = false;
        const int MAX_ITERATIONS = 100;
        int iterations = 0;

        double gamma_e = (gamma_e_lower + gamma_e_upper) / 2.0;

        double rel_tolerance = 1.0e-6;
        double abs_tolerance = 1.0e-8;

        while (fabs(gamma_e_upper - gamma_e_lower) > abs_tolerance * gamma_e) {
            if (iterations > MAX_ITERATIONS) {
                converged = false;
                break;
            }

            iterations++;

            gamma_e = (gamma_e_lower + gamma_e_upper) / 2.0;

            nu_m = floor((pow(gamma_e, 2.0) - 1.0) / (2.0 * B_ratio_electron_));
            nu_m = std::max(nu_m, 0.0);

            double summation = 0.0;
            double summation_energy_density = 0.0;
            double summation_pressure = 0.0;
            for (int nu = 0; nu <= nu_m; nu++) {
                double term_inside_square_root = pow(gamma_e, 2.0) - 1.0 - 2.0 * nu * B_ratio_electron_; // x_e(\nu)^2
                if (term_inside_square_root < 0.0 && nu > 0) break;
                summation += (nu == 0) ? sqrt(term_inside_square_root) : 2.0 * sqrt(term_inside_square_root);

                summation_energy_density += (nu == 0) ? (1.0 + 2.0 * nu * B_ratio_electron_) * psi(sqrt(term_inside_square_root) / sqrt(1.0 + 2.0 * nu * B_ratio_electron_)) : 2.0 * (1.0 + 2.0 * nu * B_ratio_electron_) * psi(sqrt(term_inside_square_root) / sqrt(1.0 + 2.0 * nu * B_ratio_electron_));
                summation_pressure += (nu == 0) ? (1.0 + 2.0 * nu * B_ratio_electron_) * eta(sqrt(term_inside_square_root) / sqrt(1.0 + 2.0 * nu * B_ratio_electron_)) : 2.0 * (1.0 + 2.0 * nu * B_ratio_electron_) * eta(sqrt(term_inside_square_root) / sqrt(1.0 + 2.0 * nu * B_ratio_electron_));
            }

            calculated_electron_density = (2.0 * B_ratio_electron_) / (pow(2.0 * M_PI, 2.0) * pow(lambda_e, 3.0)) * summation;
            calculated_electron_energy_density = (2.0 * B_ratio_electron_) / (pow(2.0 * M_PI, 2.0) * pow(lambda_e, 3.0)) * (m_electron * (c * c)) * summation_energy_density;
            calculated_electron_pressure = (2.0 * B_ratio_electron_) / (pow(2.0 * M_PI, 2.0) * pow(lambda_e, 3.0)) * (m_electron * (c * c)) * summation_pressure;
            relative_error = (calculated_electron_density - electron_density) / electron_density;

            if (fabs(relative_error) < rel_tolerance) {
                converged = true;
                break;
            } else {
                if (relative_error > 0) {
                    gamma_e_upper = gamma_e;
                } else {
                    gamma_e_lower = gamma_e;
                }
            }
        }

        double lattice_energy_density = - 1.444 * pow(Z, 2.0 / 3.0) * pow(e_charge, 2.0) * pow(calculated_electron_density, 4.0 / 3.0);
        double lattice_pressure = lattice_energy_density / 3.0;

        double total_energy_density = nucleon_density * mass_energy + calculated_electron_energy_density + lattice_energy_density;
        total_pressure = calculated_electron_pressure + lattice_pressure;
        
        total_mass_density = total_energy_density / (c * c);

        double gibbs_free_energy = (mass_energy / A) + (Z / A) * ((gamma_e * m_electron * c * c) - (m_electron * c * c)) + ((4.0 / 3.0) * ((Z * lattice_energy_density) / (A * electron_density)));

        if (gibbs_free_energy < best_composition.gibbs_free_energy && total_pressure > 0 && converged) {
            // Update best composition
            best_composition.optimal_A = A;
            best_composition.optimal_Z = Z;
            best_composition.optimal_element = nucleus.element;
            best_composition.gibbs_free_energy = gibbs_free_energy;
            best_composition.total_energy_density = total_energy_density;
            best_composition.total_pressure = total_pressure;
            best_composition.total_mass_density = total_mass_density;
            best_composition.electron_density = calculated_electron_density;
            best_composition.electron_energy_density = calculated_electron_energy_density;
            best_composition.electron_pressure = calculated_electron_pressure;
            best_composition.gamma_e = gamma_e;
            best_composition.max_landau_level = nu_m;
            best_composition.lattice_energy_density = lattice_energy_density;
            best_composition.lattice_pressure = lattice_pressure;
            best_composition.converged = true;
            best_composition.relative_error = fabs(relative_error);
            best_composition.iterations = iterations;
        } 
    } // End loop over atomic_masses_
    
    return best_composition;
}

std::vector<EquilibriumComposition> MagneticBPSEOS::runSimulation(double nB_min, double nB_max, int num_points) {
    std::vector<EquilibriumComposition> results;
    std::vector<double> nB_list = generateBaryonDensities(nB_min, nB_max, num_points);
    for (const double nB : nB_list) {
        EquilibriumComposition comp = computeEquilibriumComposition(nB);
        results.push_back(comp);
    }
    return results;
}