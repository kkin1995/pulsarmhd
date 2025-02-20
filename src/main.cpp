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

void print_equilibrium_table(const std::vector<EquilibriumComposition>& results) {
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

// Function to read the atomic masses from a CSV file
std::vector<AtomicMass> read_atomic_masses(const std::string& filename) {
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
double get_atomic_mass(int A, int Z, const std::vector<AtomicMass>& atomicMasses) {
    for (const auto& entry : atomicMasses) {
        if (entry.A == A && entry.Z == Z) {
            return entry.M;
        }
    }
    std::cerr << "Error: No data found for A=" << A << " Z=" << Z << std::endl;
    return -1.0; // Return an error value if not found
}

double psi(double x) {
    if (fabs(x) < 1e-12) return 0.0;

    double first_term = (1.0 / 2.0) * x * sqrt(1.0 + (x * x));
    double second_term = (1.0 / 2.0) * std::log(x + sqrt(1.0 + (x * x)));

    return first_term + second_term;
}

double eta(double x) {
    if (fabs(x) < 1e-12) return 0.0;

    double first_term = (1.0 / 2.0) * x * sqrt(1.0 + (x * x));
    double second_term = (1.0 / 2.0) * std::log(x + sqrt(1.0 + (x * x)));

    return first_term - second_term;
}

int main() {
    constexpr double hbar = 1.0546e-27;
    constexpr double c = 2.9979e10;
    constexpr double m_electron = 9.1094e-28;
    constexpr double e_charge = 4.80320471e-10; // CGS units
    constexpr double lambda_e = hbar / (m_electron * c);

    std::string filename = "/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/atomic_masses.csv";
    auto atomicMasses = read_atomic_masses(filename);

    if (atomicMasses.empty()) {
        std::cerr << "Error: No atomic mass data loaded." << std::endl;
        return 1;
    }

    std::vector<double> n_B_list;
    double n_B_min = 1.0e22; // Before neutron appearance
    double n_B_max = 1.0e35; // Well into neutron dominance
    int num_points = 100;    // Higher resolution

    // Logarithmic spacing for smoother transitions
    for (int i = 0; i < num_points; i++) {
        double n_B_i = n_B_min * pow(10, i * (log10(n_B_max / n_B_min) / (num_points - 1)));
        n_B_list.push_back(n_B_i);
    }

    double B_ratio_electron = 1.0;

    std::vector<EquilibriumComposition> equilibrium_results;

    std::string output_filename = "magnetic_bps_equation_of_state.csv";
    std::ofstream outfile(output_filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << output_filename << std::endl;
        return 1;
    }
    outfile << "log_n,log_rho,log_P\n";

    for (const double nB : n_B_list) {
        EquilibriumComposition best_composition;
        best_composition.baryon_density = nB;
        best_composition.gibbs_free_energy = std::numeric_limits<double>::max();
        best_composition.optimal_A = 0;
        best_composition.optimal_Z = 0;
        best_composition.optimal_element = "None";
        best_composition.converged = false;

        double total_pressure = 0.0;
        double total_mass_density = 0.0;

        for (const auto& nucleus : atomicMasses) {
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
            double gamma_e_lower = 1.0 + 1e-5;
            double gamma_e_upper = std::max(gamma_e_lower, sqrt(1.0 + x_e * x_e));

            double nu_m = 0.0;
            double relative_error = 0.0;
            bool converged = false;
            const int MAX_ITERATIONS = 100;
            int iterations = 0;

            double gamma_e = (gamma_e_lower + gamma_e_upper) / 2.0;
            while (fabs(gamma_e_upper - gamma_e_lower) > 1.0e-6) {
                if (iterations > MAX_ITERATIONS) {
                    converged = false;
                    break;
                }

                iterations++;

                nu_m = floor((pow(gamma_e, 2.0) - 1.0) / (2.0 * B_ratio_electron));
                nu_m = std::max(nu_m, 0.0);

                double summation = 0.0;
                double summation_energy_density = 0.0;
                double summation_pressure = 0.0;
                for (int nu = 0; nu <= nu_m; nu++) {
                    double term_inside_square_root = pow(gamma_e, 2.0) - 1.0 - 2.0 * nu * B_ratio_electron; // x_e(\nu)^2
                    if (term_inside_square_root < 0.0 && nu > 0) break;
                    summation += (nu == 0) ? sqrt(term_inside_square_root) : 2.0 * sqrt(term_inside_square_root);

                    summation_energy_density += (nu == 0) ? (1.0 + 2.0 * nu * B_ratio_electron) * psi(sqrt(term_inside_square_root) / sqrt(1.0 + 2.0 * nu * B_ratio_electron)) : 2.0 * (1.0 + 2.0 * nu * B_ratio_electron) * psi(sqrt(term_inside_square_root) / sqrt(1.0 + 2.0 * nu * B_ratio_electron));
                    summation_pressure += (nu == 0) ? (1.0 + 2.0 * nu * B_ratio_electron) * eta(sqrt(term_inside_square_root) / sqrt(1.0 + 2.0 * nu * B_ratio_electron)) : 2.0 * (1.0 + 2.0 * nu * B_ratio_electron) * eta(sqrt(term_inside_square_root) / sqrt(1.0 + 2.0 * nu * B_ratio_electron));
                }

                calculated_electron_density = (2.0 * B_ratio_electron) / (pow(2.0 * M_PI, 2.0) * pow(lambda_e, 3.0)) * summation;
                calculated_electron_energy_density = (2.0 * B_ratio_electron) / (pow(2.0 * M_PI, 2.0) * pow(lambda_e, 3.0)) * (m_electron * (c * c)) * summation_energy_density;
                calculated_electron_pressure = (2.0 * B_ratio_electron) / (pow(2.0 * M_PI, 2.0) * pow(lambda_e, 3.0)) * (m_electron * (c * c)) * summation_pressure;
                relative_error = (calculated_electron_density - electron_density) / electron_density;

                if (fabs(relative_error) < 1.0e-6) {
                    converged = true;
                    break;
                } else {
                    if (relative_error > 0) {
                        gamma_e_upper = gamma_e;
                    } else {
                        gamma_e_lower = gamma_e;
                    }
                    gamma_e = (gamma_e_lower + gamma_e_upper) / 2.0;
                }
            }

            double lattice_energy_density = - 1.444 * pow(Z, 2.0 / 3.0) * pow(e_charge, 2.0) * pow(calculated_electron_density, 4.0 / 3.0);
            double lattice_pressure = lattice_energy_density / 3.0;

            double total_energy_density = nucleon_density * mass + calculated_electron_energy_density + lattice_energy_density;
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
        }
        // Check if any valid nucleus was found
        if (best_composition.optimal_A == 0) {
            std::cout << "  Warning: No valid nucleus found for density " << nB << "\n";
        } else {
            std::cout << "  Optimal nucleus at density " << nB << ": " 
                    << best_composition.optimal_element << "-" 
                    << best_composition.optimal_A << "\n";
        }
        // Store results for this density
        equilibrium_results.push_back(best_composition);

        outfile << std::log10(nB) << "," << std::log10(total_mass_density) << "," << std::log10(total_pressure) << "\n";
    }
    print_equilibrium_table(equilibrium_results);
    outfile.close();
}