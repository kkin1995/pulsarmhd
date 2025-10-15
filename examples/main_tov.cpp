#include "non_rotating_stellar_structure.hpp"
#include "polytropic_eos.hpp"
#include "rk4.hpp"

#include <algorithm> // for std::replace
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip> // for std::setprecision, std::scientific
#include <iostream>
#include <sstream> // for std::ostringstream
#include <string>
#include <vector>

int main() {
  double pi = M_PI;
  double fraction = 4.0 / 3.0;

  // NEW: Use polytropic EOS type instead of hardcoded parameters
  PolytropicGasType eos_type = PolytropicGasType::NEUTRON_NON_RELATIVISTIC;
  double mu_e = 2.0; // Default mean molecular weight per electron

  double r_start = 10.0; // cm
  double r_end = 1e10;   // cm
  double log_r_start = log10(r_start);
  double log_r_end = log10(r_end);
  double dlogr = 0.001;

  std::vector<double> nonrelativistic_electrons_central_densities = {1e0, 1e1, 1e2, 1e3, 1e4};
  std::vector<double> relativistic_electrons_central_densities = {1e5, 1e6, 1e7, 1e8, 1e9};
  std::vector<double> nonrelativistic_neutrons_central_densities = {1e10, 1e11, 1e12, 1e13, 1e14};
  std::vector<double> relativistic_neutrons_central_densities = {1e15, 1e16, 1e17, 5e17,
                                                                 1e18, 5e18, 1e19};

  std::vector<double> central_densities;
  int num_points = 25; // Change this to get a finer grid
  double start = 1e15;
  double end = 1e18;
  double log_start = std::log10(start);
  double log_end = std::log10(end);
  double step = (log_end - log_start) / (num_points - 1);
  for (int i = 0; i < num_points; ++i) {
    double rho = std::pow(10.0, log_start + i * step);
    central_densities.push_back(rho);
  }

  // You can uncomment and modify these to use different EOS types and density ranges:
  if (eos_type == PolytropicGasType::ELECTRON_NON_RELATIVISTIC) {
    central_densities = nonrelativistic_electrons_central_densities;
  } else if (eos_type == PolytropicGasType::ELECTRON_RELATIVISTIC) {
    central_densities = relativistic_electrons_central_densities;
  } else if (eos_type == PolytropicGasType::NEUTRON_NON_RELATIVISTIC) {
    central_densities = nonrelativistic_neutrons_central_densities;
  } else if (eos_type == PolytropicGasType::NEUTRON_RELATIVISTIC) {
    central_densities = relativistic_neutrons_central_densities;
  }

  for (size_t i = 0; i < central_densities.size(); i++) {
    std::cout << "Processing " << (i + 1) << " of " << central_densities.size() << " densities..."
              << std::endl;

    double rho_c = central_densities[i];

    // NEW: Updated function call with polytropic EOS type and mu_e
    std::tuple<int, double> result =
        non_rotating_stellar_structure(eos_type, rho_c, r_start, r_end, dlogr, mu_e);
    int idx = std::get<0>(result);
    double log_mass = std::get<1>(result);

    std::cout << "Total steps: " << idx << std::endl;
    std::cout << "Final mass: " << pow(10.0, log_mass) / 1.989e33 << " M_Sun" << std::endl;
  }
}
