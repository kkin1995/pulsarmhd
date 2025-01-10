#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>
#include <sstream>      // for std::ostringstream
#include <iomanip>      // for std::setprecision, std::scientific
#include <algorithm>    // for std::replace
#include <string>
#include "rk4.hpp"
#include "non_rotating_stellar_structure.hpp"

int main() {
    double pi = M_PI;
    double fraction = 4.0 / 3.0;
    double k;
    double gamma;
    std::string name;
    DegenerateGasType type = NEUTRON_RELATIVISTIC;
    set_eos_parameters(type, name, k, gamma);
    
    double r_start = 10.0; // cm
    double r_end = 1e10;// cm
    double log_r_start = log10(r_start);
    double log_r_end = log10(r_end);
    double dlogr = 0.001;

    std::vector<double> nonrelativistic_electrons_central_densities = {1e0, 1e1, 1e2, 1e3, 1e4};
    std::vector<double> relativistic_electrons_central_densities = {1e5, 1e6, 1e7, 1e8, 1e9};
    std::vector<double> nonrelativistic_neutrons_central_densities = {1e10, 1e11, 1e12, 1e13, 1e14};
    std::vector<double> relativistic_neutrons_central_densities = {1e15, 1e16, 1e17, 5e17, 1e18, 5e18, 1e19};

    std::vector<double> central_densities;
    if (type == ELECTRON_NON_RELATIVISTIC) {
        central_densities = nonrelativistic_electrons_central_densities;
    } else if (type == ELECTRON_RELATIVISTIC) {
        central_densities = relativistic_electrons_central_densities;
    } else if (type == NEUTRON_NON_RELATIVISTIC) {
        central_densities = nonrelativistic_neutrons_central_densities;
    } else if (type == NEUTRON_RELATIVISTIC) {
        central_densities = relativistic_neutrons_central_densities;
    }

    for (size_t i = 0; i < central_densities.size(); i++) {
        std::cout << "Processing " << (i+1) << " of " << central_densities.size() 
          << " densities..." << std::endl;

        double rho_c = central_densities[i];

        std::cout << "Central Density: " << rho_c << " g/cm^3" << std::endl;

        std::string filename = get_filename(name, rho_c);
        std::ofstream outfile(filename);
        if (!outfile.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            continue;
        }
        std::cout << "Writing to file: " << filename << std::endl;

        double log_m0 = log10(fraction) + log10(pi) + 3.0*log10(r_start) + log10(rho_c);
        double log_p0 = log10(k) + (gamma * log10(rho_c)); 
        std::vector<double> state = {log_m0, log_p0};

        std::cout << "Initial conditions:" << std::endl;
        std::cout << "  log_r_start = " << log_r_start << std::endl;
        std::cout << "  log_m0 = " << log_m0 << std::endl;
        std::cout << "  log_p0 = " << log_p0 << std::endl;
        std::cout << "  k = " << k << ", gamma = " << gamma << std::endl;

        outfile << "log_r[cm],log_m[g],log_P[dyne/cm^2]\n";

        int idx = 0;
        double log_r = log_r_start; 
        while (log_r < log_r_end) {
            state = rk4_step(log_r, dlogr, state, 
            [k, gamma](double r, const std::vector<double>& s) {
                return tolman_oppenheimer_volkoff_derivatives(r, s, k, gamma);
            });
            log_r += dlogr;
            outfile << log_r << "," << state[0] << "," << state[1] << "\n";

            // Stop if pressure drops below a threshold
            if (state[1] < log10(1e-8)) {
                std::cout << "Surface reached at log_r: " << log_r << std::endl;
                break;
            }

            if (!std::isfinite(state[0]) || !std::isfinite(state[1])) {
                std::cout << "Non-finite values at log_r: " << log_r << std::endl;
                break;
            }

            if (idx % 100 == 0) { 
                std::cout << "log_r: " << log_r << ", log_m: " << state[0] 
                        << ", log_P: " << state[1] << std::endl;
            }

            idx += 1;
        }
        outfile.close();
        // At end of integration
        std::cout << "Total steps: " << idx << std::endl;
        std::cout << "Final mass: " << pow(10.0, state[0]) / 1.989e33 << " M_Sun" << std::endl;
    }
}