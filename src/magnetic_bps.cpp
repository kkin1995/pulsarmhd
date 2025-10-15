/**
 * Implementation of the magnetic BPS equation of state.
 * Based on: Lai & Shapiro (1991), ApJ, 383, 745
 */

#include "magnetic_bps.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream> // ifstream/ofstream
#include <iomanip> // setw, setprecision
#include <iostream>
#include <sstream> // stringstream
#include <vector>
#ifdef _OPENMP
#include <cstdlib> // for getenv if you want to print env vars
#include <omp.h>
#endif

// Constructor: Initializes configuration parameters and loads atomic mass data.
MagneticBPSEOS::MagneticBPSEOS(const std::string &atomic_mass_file, double B_ratio_electron,
                               double rel_tol, double abs_tol)
    : B_ratio_electron_(B_ratio_electron), rel_tolerance_(rel_tol), abs_tolerance_(abs_tol),
      atomic_masses_(readAtomicMasses(atomic_mass_file)) {
  if (atomic_masses_.empty()) {
    std::cerr << "Error: Atomic mass data could not be loaded from " << atomic_mass_file
              << std::endl;
  }
}

std::vector<double> MagneticBPSEOS::generateBaryonDensities(double nB_min, double nB_max,
                                                            int num_points) {
  std::vector<double> nB_list;
  // Logarithmic spacing for smoother transitions
  for (int i = 0; i < num_points; i++) {
    double nB_i = nB_min * pow(10, i * (log10(nB_max / nB_min) / (num_points - 1)));
    nB_list.push_back(nB_i);
  }
  return nB_list;
}

std::vector<AtomicMass> MagneticBPSEOS::readAtomicMasses(const std::string &filename) {
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

    std::getline(ss, temp, ',');
    entry.Z = std::stoi(temp);
    std::getline(ss, temp, ',');
    entry.A = std::stoi(temp);
    std::getline(ss, entry.element, ',');
    std::getline(ss, temp, ',');
    entry.dM = std::stod(temp);
    std::getline(ss, temp, ',');
    entry.BE = std::stod(temp);
    std::getline(ss, temp, ',');
    entry.BDE = std::stod(temp);
    std::getline(ss, temp, ',');
    entry.M = std::stod(temp);

    atomicMasses.push_back(entry);
  }

  file.close();
  return atomicMasses;
}

// Function to find M(A, Z) for given A and Z
double MagneticBPSEOS::getAtomicMass(int A, int Z) const {
  for (const auto &entry : atomic_masses_) {
    if (entry.A == A && entry.Z == Z) {
      return entry.M;
    }
  }
  std::cerr << "Error: No data found for A=" << A << " Z=" << Z << std::endl;
  return -1.0; // Return an error value if not found
}

double MagneticBPSEOS::psi(double x) {
  if (fabs(x) < 1e-12)
    return 0.0;

  double first_term = (1.0 / 2.0) * x * sqrt(1.0 + (x * x));
  double second_term = (1.0 / 2.0) * std::log(x + sqrt(1.0 + (x * x)));

  return first_term + second_term;
}

double MagneticBPSEOS::eta(double x) {
  if (fabs(x) < 1e-12)
    return 0.0;

  double first_term = (1.0 / 2.0) * x * sqrt(1.0 + (x * x));
  double second_term = (1.0 / 2.0) * std::log(x + sqrt(1.0 + (x * x)));

  return first_term - second_term;
}

void MagneticBPSEOS::printEquilibriumTable(
    const std::vector<EquilibriumComposition> &results) const {
  std::cout << "\n=== Equilibrium Composition Results ===\n";
  std::cout << std::scientific << std::setprecision(4);
  std::cout << std::setw(12) << "Density" << std::setw(8) << "A" << std::setw(8) << "Z"
            << std::setw(8) << "Element" << std::setw(15) << "Gibbs Energy" << std::setw(15)
            << "Pressure" << std::setw(15) << "Mass Density" << std::setw(12) << "ν_max"
            << std::endl;

  for (const auto &comp : results) {
    std::cout << std::setw(12) << comp.baryon_density << std::setw(8) << comp.optimal_A
              << std::setw(8) << comp.optimal_Z << std::setw(8) << comp.optimal_element
              << std::setw(15) << comp.gibbs_free_energy << std::setw(15) << comp.total_pressure
              << std::setw(15) << comp.total_mass_density;
    // Print Landau level separately to ensure proper spacing
    std::cout << std::setw(12) << comp.max_landau_level << (comp.converged ? "" : " *")
              << std::endl;
  }
}

void MagneticBPSEOS::writeEOSResults(const std::string &output_file,
                                     const std::vector<EquilibriumComposition> &results) const {
  std::ofstream outfile(output_file);
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not open file " << output_file << std::endl;
    return;
  }
  outfile << "log_n,log_rho,log_P\n";

  for (const auto &comp : results) {
    outfile << std::log10(comp.baryon_density) << "," << std::log10(comp.total_mass_density) << ","
            << std::log10(comp.total_pressure) << "\n";
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

  const double rel_tolerance = 1.0e-6;
  const double abs_tolerance = 1.0e-8;
  const int MAX_ITERATIONS = 100;

  for (int i = 0; i < (int)atomic_masses_.size(); ++i) {
    const auto &nucleus = atomic_masses_[i];
    const int A = nucleus.A;
    const int Z = nucleus.Z;
    const double mass = nucleus.M;
    if (mass <= 0)
      continue;

    const double mass_energy = mass * 931.494 * 1.602e-6; // Convert u to MeV, then to erg

    const double nucleon_density = nB / A;
    const double electron_density = (Z * nB) / A;

    double calculated_electron_density = 0.0;
    double calculated_electron_energy_density = 0.0;
    double calculated_electron_pressure = 0.0;

    // Solve Coupled Equations
    const double x_e = pow(3.0 * pow(M_PI, 2.0) * pow(lambda_e, 3.0) * electron_density, 1.0 / 3.0);

    double gamma_e_lower = 1.0 + 1e-10;
    double gamma_e_upper = std::max(gamma_e_lower, sqrt(1.0 + x_e * x_e));
    double gamma_e = 0.5 * (gamma_e_lower + gamma_e_upper);

    int nu_m = 0;
    double relative_error = 0.0;
    bool converged = false;
    int iterations = 0;

    if (B_ratio_electron_ <= 0.0) {
      // Charge neutrality gives n_e directly
      calculated_electron_density = electron_density;
      const double x =
          pow(3.0 * M_PI * M_PI * pow(lambda_e, 3.0) * calculated_electron_density, 1.0 / 3.0);
      const auto asinhx = [](double z) { return log(z + sqrt(1.0 + z * z)); };

      double pref = (m_electron * c * c) / (8.0 * M_PI * M_PI * pow(lambda_e, 3.0));
      calculated_electron_pressure =
          pref * (x * ((2.0 * x * x) / 3.0 - 1.0) * sqrt(1.0 + x * x) + asinhx(x));
      calculated_electron_energy_density =
          pref * (x * (2.0 * x * x + 1.0) * sqrt(1.0 + x * x) - asinhx(x));

      gamma_e = sqrt(1.0 + x * x);

      nu_m = 0; // continuum limit, but store 0 for bookkeeping
      relative_error = 0.0;
      converged = true; // closed-form, no iteration

      iterations = 1; // for reporting
    } else {
      // Degenerate collapsed-bounds guard: treat as ν_max = 0 case
      if (std::fabs(gamma_e_upper - gamma_e_lower) <= abs_tolerance * gamma_e) {
        gamma_e = gamma_e_upper;
        nu_m = 0;
        // ν=0 sums (identical to loop, but explicit)
        const double x2 = gamma_e * gamma_e - 1.0; // x_e(ν=0)^2
        const double x = (x2 > 0.0) ? std::sqrt(x2) : 0.0;

        // Summations for ν=0
        double summation = x;                                    // (ν=0) no factor 2
        double summation_energy_state = psi(x / std::sqrt(1.0)); // psi(x)
        double summation_pressure_st = eta(x / std::sqrt(1.0));  // eta(x)

        const double pref =
            (2.0 * B_ratio_electron_) / (4.0 * M_PI * M_PI * std::pow(lambda_e, 3.0));
        calculated_electron_density = pref * summation;
        calculated_electron_energy_density = pref * (m_electron * (c * c)) * summation_energy_state;
        calculated_electron_pressure = pref * (m_electron * (c * c)) * summation_pressure_st;

        relative_error = (calculated_electron_density - electron_density) / electron_density;
        converged = true;
        iterations = 1;
      } else {
        while (fabs(gamma_e_upper - gamma_e_lower) > abs_tolerance * gamma_e) {
          if (++iterations > MAX_ITERATIONS) {
            converged = false;
            break;
          }

          gamma_e = 0.5 * (gamma_e_lower + gamma_e_upper);

          nu_m = std::max(
              0, static_cast<int>(floor((gamma_e * gamma_e - 1.0) / (2.0 * B_ratio_electron_))));

          double summation = 0.0;
          double summation_energy_density = 0.0;
          double summation_pressure = 0.0;
          for (int nu = 0; nu <= nu_m; nu++) {
            double term_inside_square_root =
                pow(gamma_e, 2.0) - 1.0 - 2.0 * nu * B_ratio_electron_; // x_e(\nu)^2
            if (term_inside_square_root < 0.0 && nu > 0)
              break;
            summation +=
                (nu == 0) ? sqrt(term_inside_square_root) : 2.0 * sqrt(term_inside_square_root);

            summation_energy_density += (nu == 0)
                                            ? (1.0 + 2.0 * nu * B_ratio_electron_) *
                                                  psi(sqrt(term_inside_square_root) /
                                                      sqrt(1.0 + 2.0 * nu * B_ratio_electron_))
                                            : 2.0 * (1.0 + 2.0 * nu * B_ratio_electron_) *
                                                  psi(sqrt(term_inside_square_root) /
                                                      sqrt(1.0 + 2.0 * nu * B_ratio_electron_));
            summation_pressure += (nu == 0) ? (1.0 + 2.0 * nu * B_ratio_electron_) *
                                                  eta(sqrt(term_inside_square_root) /
                                                      sqrt(1.0 + 2.0 * nu * B_ratio_electron_))
                                            : 2.0 * (1.0 + 2.0 * nu * B_ratio_electron_) *
                                                  eta(sqrt(term_inside_square_root) /
                                                      sqrt(1.0 + 2.0 * nu * B_ratio_electron_));
          }

          calculated_electron_density =
              (2.0 * B_ratio_electron_) / (4.0 * M_PI * M_PI * pow(lambda_e, 3.0)) * summation;
          calculated_electron_energy_density = (2.0 * B_ratio_electron_) /
                                               (4.0 * M_PI * M_PI * pow(lambda_e, 3.0)) *
                                               (m_electron * (c * c)) * summation_energy_density;
          calculated_electron_pressure = (2.0 * B_ratio_electron_) /
                                         (4.0 * M_PI * M_PI * pow(lambda_e, 3.0)) *
                                         (m_electron * (c * c)) * summation_pressure;
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
        if (!converged) {
          continue; // try next nucleus
        }
      }
    }
    double lattice_energy_density = -1.444 * pow(Z, 2.0 / 3.0) * e_charge * e_charge *
                                    pow(calculated_electron_density, 4.0 / 3.0);
    double lattice_pressure = lattice_energy_density / 3.0;

    double total_energy_density =
        nucleon_density * mass_energy + calculated_electron_energy_density + lattice_energy_density;
    const double total_pressure = calculated_electron_pressure + lattice_pressure;

    const double total_mass_density = total_energy_density / (c * c);

    double gibbs_free_energy =
        (mass_energy / A) + (Z / A) * ((gamma_e * m_electron * c * c) - (m_electron * c * c)) +
        ((4.0 / 3.0) * ((Z * lattice_energy_density) / (A * calculated_electron_density)));

    if (gibbs_free_energy < best_composition.gibbs_free_energy && total_pressure > 0 && converged) {
      // Update best composition
      best_composition.optimal_A = A;
      best_composition.optimal_Z = Z;
      best_composition.optimal_element = nucleus.element;
      best_composition.baryon_density = nB;
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
  return best_composition;
}

std::vector<EquilibriumComposition> MagneticBPSEOS::runSimulation(double nB_min, double nB_max,
                                                                  int num_points) {
  std::vector<EquilibriumComposition> results(num_points);
  std::vector<double> nB_list = generateBaryonDensities(nB_min, nB_max, num_points);

  // progress + stats
  std::atomic<int> completed(0);
  const int report_step = std::max(1, num_points / 20); // ~5% steps

  long long sum_iters = 0;
  long long sum_nu = 0;
  int max_iters = 0;
  int max_nu = 0;

#ifdef _OPENMP
  double t0 = omp_get_wtime();
#else
  auto t0 = std::chrono::high_resolution_clock::now();
#endif

#pragma omp parallel for schedule(dynamic, 2) reduction(+ : sum_iters, sum_nu)                     \
    reduction(max : max_iters, max_nu)
  for (int i = 0; i < num_points; ++i) {
    auto comp = computeEquilibriumComposition(nB_list[i]);
    results[i] = comp;

    // collect stats
    sum_iters += comp.iterations;
    sum_nu += static_cast<int>(comp.max_landau_level);
    if (comp.iterations > max_iters)
      max_iters = comp.iterations;
    // Landau level is an integer count; round to nearest and clamp at >= 0
    const int level = static_cast<int>(std::lround(comp.max_landau_level));
    if (level > max_nu)
      max_nu = level;

    int done = completed.fetch_add(1, std::memory_order_relaxed) + 1;

    if (done % report_step == 0 || done == num_points) {
#ifdef _OPENMP
      double elapsed = omp_get_wtime() - t0;
#else
      double elapsed =
          std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
#endif

#ifdef _OPENMP
#pragma omp critical
#endif
      {
        const double pct = 100.0 * double(done) / double(num_points);
        std::cout << "[progress] " << done << "/" << num_points << " (" << std::fixed
                  << std::setprecision(1) << pct << "%)"
                  << "  elapsed " << std::setprecision(2) << elapsed << " s\n";
        std::cout.unsetf(std::ios::floatfield); // optional: reset format
      }
    }
  }

#ifdef _OPENMP
  double t1 = omp_get_wtime();
  const double avg_iters = (num_points ? double(sum_iters) / num_points : 0.0);
  const double avg_nu = (num_points ? double(sum_nu) / num_points : 0.0);

  std::cout << "[timing] runSimulation took " << (t1 - t0) << " s with " << omp_get_max_threads()
            << " threads\n"
            << "[stats] iterations avg=" << avg_iters << " max=" << max_iters
            << " | nu_max avg=" << avg_nu << " max=" << max_nu << "\n";
#endif

  return results;
}
