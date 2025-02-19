// main.cpp
#include "non_magnetic_ideal_npe_gas.hpp"

int main() {
    // Create an instance of our EOS calculator
    // true enables debug output
    NonMagneticNPEGas eos(true);

    // Define the density range and resolution
    double nB_min = 1.0e27;  // Before neutron appearance
    double nB_max = 1.0e35;  // Well into neutron dominance
    int num_points = 300;    // Resolution

    // Calculate EOS and write to file
    std::string filename = "ideal_non_magnetic_npe_gas_eos.csv";
    eos.calculateEOS(filename, nB_min, nB_max, num_points);

    return 0;
}