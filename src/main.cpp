#include "magnetic_bps.hpp"
#include <iostream>
#include <vector>
#include <string>

int main() {
    // Define the atomic mass data file path
    std::string atomic_mass_file = "/mnt/c/Users/karan/Dropbox/KARAN/2 Areas/Education/PhD/3 Research/pulsarmhd/data/atomic_masses.csv";

    // Magnetic field strength ratio for electrons
    constexpr double B_ratio_electron = 100.0;

    // Create an instance of MagneticBPSEOS
    MagneticBPSEOS eos(atomic_mass_file, B_ratio_electron);

    // Define baryon density range
    double nB_min = 1.0e22;  // Before neutron appearance
    double nB_max = 1.0e35;  // Well into neutron dominance
    int num_points = 250;    // Higher resolution

    // Run the equilibrium composition simulation
    std::vector<EquilibriumComposition> equilibrium_results = eos.runSimulation(nB_min, nB_max, num_points);

    // Output results to a CSV file
    std::string output_filename = "magnetic_bps_eos_B_100.csv";
    eos.writeEOSResults(output_filename, equilibrium_results);

    // Print the results in tabular format
    eos.printEquilibriumTable(equilibrium_results);

    return 0;
}
