#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

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

    int A = 56;
    int Z = 26;
    double mass = get_atomic_mass(A, Z, atomicMasses);

    if (mass != -1.0) {
        std::cout << "Atomic mass for A=" << A << " Z=" << Z << " is " << mass << " u." << std::endl;
    } else {
        std::cerr << "Error: Atomic mass not found.\n";
        return 1;
    }

    std::vector<double> n_B_list;
    double n_B_min = 1.0e22; // Before neutron appearance
    double n_B_max = 1.0e35; // Well into neutron dominance
    int num_points = 5;    // Higher resolution

    // Logarithmic spacing for smoother transitions
    for (int i = 0; i < num_points; i++) {
        double n_B_i = n_B_min * pow(10, i * (log10(n_B_max / n_B_min) / (num_points - 1)));
        n_B_list.push_back(n_B_i);
    }

    double B_ratio_electron = 0.1;
    double B_ratio_proton = 0.1;

    for (const double nB : n_B_list) {
        double calculated_electron_density;
        double calculated_electron_energy_density;
        double calculated_electron_pressure;

        double nucleon_density = nB / A;
        double electron_density = (Z * nB) / A;
        std::cout << "Baryon Density: " << nB << " cm^-3\n";
        std::cout << "Nucleon Density: " << nucleon_density << " cm^-3\n";
        std::cout << "Electron Density: " << electron_density << " cm^-3\n";

        // Solve Coupled Equations
        double x_e = pow(3.0 * pow(M_PI, 2.0) * pow(lambda_e, 3.0) * electron_density, 1.0 / 3.0);
        double gamma_e_lower = 1.0 + 1e-5;
        double gamma_e_upper = std::max(gamma_e_lower, sqrt(1.0 + x_e * x_e));
        
        std::cout << "x_e: " << x_e << "\n";
        std::cout << "Before While: gamma_e_lower: " << gamma_e_lower << "\n";
        std::cout << "Before While: gamma_e_upper: " << gamma_e_upper << "\n";

        double gamma_e = (gamma_e_lower + gamma_e_upper) / 2.0;
        while (fabs(gamma_e_upper - gamma_e_lower) > 1.0e-6) {
            std::cout << "gamma_e_lower: " << gamma_e_lower << "\n";
            std::cout << "gamma_e_upper: " << gamma_e_upper << "\n";
            std::cout << "gamma_e: " << gamma_e << "\n";

            double nu_m = floor((pow(gamma_e, 2.0) - 1.0) / (2.0 * B_ratio_electron));
            nu_m = std::max(nu_m, 0.0);
            std::cout << "nu_m: " << nu_m << "\n";

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
            std::cout << "Calculated Electron Density: " << calculated_electron_density << "\n";

            calculated_electron_energy_density = (2.0 * B_ratio_electron) / (pow(2.0 * M_PI, 2.0) * pow(lambda_e, 3.0)) * (m_electron * (c * c)) * summation_energy_density;
            std::cout << "Calculated Electron Energy Density: " << calculated_electron_energy_density << "\n";

            calculated_electron_pressure = (2.0 * B_ratio_electron) / (pow(2.0 * M_PI, 2.0) * pow(lambda_e, 3.0)) * (m_electron * (c * c)) * summation_pressure;
            std::cout << "Calculated Electron Pressure: " << calculated_electron_pressure << "\n";

            double relative_error = (calculated_electron_density - electron_density) / electron_density;
            std::cout << "Relative Error: " << relative_error << "\n";

            if (fabs(relative_error) < 1.0e-6) {
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
        std::cout << "Lattice Energy Density: " << lattice_energy_density << "\n";
        std::cout << "Lattice Pressure: " << lattice_pressure << "\n";

        double total_energy_density = nucleon_density * mass + calculated_electron_energy_density + lattice_energy_density;
        double total_pressure = calculated_electron_pressure + lattice_pressure;
        std::cout << "Total Energy Density: " << total_energy_density << "\n";
        std::cout << "Total Pressure: " << total_pressure << "\n";
    }
}