#include "non_rotating_stellar_structure.hpp"
#include "rk4.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Helpers for reporting (avoid magic numbers)
namespace {
// CGS constants used elsewhere:
//   extern const double G; // from non_rotating_stellar_structure.hpp
//   extern const double c; // from non_rotating_stellar_structure.hpp
constexpr double kMsun_g = 1.98847e33; // g
constexpr double kCmPerKm = 1e5;       // cm per km

inline double twoGM_over_Rc2(double mass_solar, double radius_km) {
  const double M = mass_solar * kMsun_g;
  const double R = radius_km * kCmPerKm;
  return (2.0 * G * M) / (R * c * c);
}
} // namespace

bool readEoSData(const std::string &filename, std::vector<double> &log_rho,
                 std::vector<double> &log_P,
                 std::vector<double> &epsilon /* erg/cm^3; optional, may stay empty */) {
  log_rho.clear();
  log_P.clear();
  epsilon.clear();

  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << "\n";
    return false;
  }

  std::string header;
  if (!std::getline(file, header)) {
    std::cerr << "Error: Empty file " << filename << "\n";
    return false;
  }

  // tokenize header
  std::vector<std::string> cols;
  {
    std::stringstream sh(header);
    std::string tok;
    while (std::getline(sh, tok, ',')) {
      tok.erase(std::remove_if(tok.begin(), tok.end(),
                               [](unsigned char ch) { return std::isspace(ch) || ch == '\"'; }),
                tok.end());
      std::transform(tok.begin(), tok.end(), tok.begin(),
                     [](unsigned char ch) { return std::tolower(ch); });
      cols.push_back(tok);
    }
  }
  auto find_col = [&](const std::string &name) -> int {
    for (size_t i = 0; i < cols.size(); ++i) {
      if (cols[i] == name) {
        return (int)i;
      }
    }
    return -1;
  };

  int i_logrho = find_col("log_rho");
  int i_logp = find_col("log_p");
  int i_rho = find_col("rho");
  int i_p = find_col("p");

  // allow several spellings for epsilon (energy density, erg/cm^3)
  int i_eps = find_col("epsilon");
  if (i_eps < 0) {
    i_eps = find_col("energy_density");
  }
  if (i_eps < 0) {
    i_eps = find_col("eps");
  }

  enum class Mode { USE_LOGS, USE_LINEAR, FALLBACK_FIRST_TWO } mode;
  if (i_logrho >= 0 && i_logp >= 0) {
    mode = Mode::USE_LOGS;
    std::cout << "[EoS] Using columns log_rho, log_P";
  } else if (i_rho >= 0 && i_p >= 0) {
    mode = Mode::USE_LINEAR;
    std::cout << "[EoS] Using columns rho, P (computing logs)";
  } else {
    mode = Mode::FALLBACK_FIRST_TWO;
    i_logrho = 0;
    i_logp = 1;
    std::cout << "[EoS] Header unrecognized; assuming first two columns are log_rho,log_P";
  }
  if (i_eps >= 0) {
    std::cout << " + epsilon\n";
  } else {
    std::cout << "\n";
  }

  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) {
      continue;
    }
    std::vector<std::string> toks;
    {
      std::stringstream ss(line);
      std::string t;
      while (std::getline(ss, t, ',')) {
        t.erase(std::remove_if(t.begin(), t.end(),
                               [](unsigned char ch) { return std::isspace(ch) || ch == '\"'; }),
                t.end());
        toks.push_back(t);
      }
    }
    if (toks.size() < 2) {
      continue;
    }

    try {
      double lr = 0, lp = 0;
      if (mode == Mode::USE_LOGS) {
        lr = std::stod(toks[(size_t)i_logrho]);
        lp = std::stod(toks[(size_t)i_logp]);
      } else if (mode == Mode::USE_LINEAR) {
        double rho = std::stod(toks[(size_t)i_rho]);
        double P = std::stod(toks[(size_t)i_p]);
        if (!(rho > 0.0 && P > 0.0)) {
          continue;
        }
        lr = std::log10(rho);
        lp = std::log10(P);
      } else { // fallback: first two columns are log10
        lr = std::stod(toks[0]);
        lp = std::stod(toks[1]);
      }
      if (!std::isfinite(lr) || !std::isfinite(lp)) {
        continue;
      }

      log_rho.push_back(lr);
      log_P.push_back(lp);

      if (i_eps >= 0 && (size_t)i_eps < toks.size()) {
        double e = std::stod(toks[(size_t)i_eps]); // already erg/cm^3
        if (std::isfinite(e)) {
          epsilon.push_back(e);
        } else {
          epsilon.push_back(std::numeric_limits<double>::quiet_NaN());
        }
      }
    } catch (...) {
      continue;
    }
  }
  file.close();

  if (log_rho.empty()) {
    std::cerr << "No usable rows read from " << filename << "\n";
    return false;
  }

  // if epsilon present, enforce same length as log_P (truncate if needed)
  if (!epsilon.empty() && epsilon.size() != log_P.size()) {
    epsilon.resize(std::min(epsilon.size(), log_P.size()));
  }
  return true;
}

int main() {
  /*
   * MODULAR TOV STELLAR STRUCTURE SOLVER
   *
   * This program has been refactored to use the modular TOV solver from
   * non_rotating_stellar_structure.cpp instead of embedded equations.
   *
   * BENEFITS:
   * - Cleaner, more maintainable code
   * - Reuses tested TOV equation implementations
   * - Better separation of concerns
   * - Enhanced error handling and adaptive stepping
   * - Consistent physics across the codebase
   *
   * FUNCTIONALITY:
   * - Loads realistic EOS data from CSV files
   * - Sets up GSL splines for EOS interpolation
   * - Scans over central densities to generate mass-radius relation
   * - Uses modular spline-based TOV solver for each stellar model
   * - Outputs individual stellar structure files for analysis
   */

  gsl_set_error_handler_off();

  std::vector<double> log_rho;
  std::vector<double> log_P;
  std::vector<double> epsilon;
  std::vector<double> central_densities;

  // Integration parameters
  double r_start = 1.0;       // Starting radius (cm)
  double r_end = 1.0e8;       // Maximum radius (cm)
  double base_dlogr = 0.0001; // Base step size in log(r)

  // Load EOS data
  std::string eos_filename =
      "/home/karan-kinariwala/Dropbox/KARAN/2-Areas/Education/PhD/3-Research/pulsarmhd/data/"
      "sahu_basu_datta_bbp_magnetic_bps_b_0e_00.csv";

  std::string out_file = "data/tov_solution_sbd_bbp_magnetic_bps_";

  if (readEoSData(eos_filename, log_rho, log_P, epsilon)) {

    std::cout << "Successfully read " << log_rho.size() << " EOS data points." << '\n';

    // --- sort rows by log_P and keep epsilon aligned
    struct Row {
      double lp, lr, eps;
      bool has_eps;
    };
    std::vector<Row> rows;
    rows.reserve(log_P.size());
    for (size_t i = 0; i < log_P.size(); ++i) {
      Row r{log_P[i], log_rho[i], 0.0, false};
      if (!epsilon.empty() && i < epsilon.size()) {
        r.eps = epsilon[i];
        r.has_eps = true;
      }
      rows.push_back(r);
    }
    std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b) { return a.lp < b.lp; });

    // --- de-duplicate to ensure strict monotonic log_P
    std::vector<double> lp2, lr2, eps2;
    auto push_row = [&](const Row &r) {
      lp2.push_back(r.lp);
      lr2.push_back(r.lr);
      if (!epsilon.empty()) {
        eps2.push_back(r.has_eps ? r.eps : std::numeric_limits<double>::quiet_NaN());
      }
    };
    if (!rows.empty()) {
      push_row(rows.front());
    }
    for (size_t i = 1; i < rows.size(); ++i) {
      if (rows[i].lp > lp2.back()) {
        push_row(rows[i]);
      } // strict >
    }

    log_P.swap(lp2);
    log_rho.swap(lr2);
    if (!epsilon.empty()) {
      epsilon.swap(eps2);
    }

    if (log_P.size() < 2) {
      std::cerr << "EOS has fewer than 2 unique points after cleaning.\n";
      return 1;
    }
    std::cout << "After cleaning: " << log_P.size() << " unique data points.\n";

    // ---- Build central_densities safely from EOS bounds ----
    if (log_rho.size() < 2) {
      std::cerr << "[ScanBuild] EOS has < 2 points after cleaning (size=" << log_rho.size()
                << "). Cannot build rho_c scan.\n";
      return 1;
    }

    // Use actual min/max in case the vector isn’t strictly sorted
    auto [lr_min_it, lr_max_it] = std::minmax_element(log_rho.begin(), log_rho.end());
    const double lr_min = *lr_min_it;
    const double lr_max = *lr_max_it;

    if (!std::isfinite(lr_min) || !std::isfinite(lr_max)) {
      std::cerr << "[ScanBuild] Non-finite log_rho bounds: " << lr_min << ", " << lr_max << "\n";
      return 1;
    }

    // Convert to linear densities (g/cm^3)
    const double rho_min_eos = std::pow(10.0, lr_min);
    const double rho_max_eos = std::pow(10.0, lr_max);
    if (!(rho_max_eos > rho_min_eos)) {
      std::cerr << "[ScanBuild] Bad EOS range: rho_min=" << rho_min_eos
                << " rho_max=" << rho_max_eos << "\n";
      return 1;
    }

    // Stay inside the table and above some physical floor
    const double scan_lo = std::max(1e14, rho_min_eos * 1.02); // nudge in from the edge
    const double scan_hi = rho_max_eos * 0.98;

    if (!(scan_hi > scan_lo)) {
      std::cerr << "[ScanBuild] scan_hi <= scan_lo after clamping: " << scan_lo << " .. " << scan_hi
                << "\n";
      return 1;
    }

    int num_points = 30;
    central_densities.clear();
    central_densities.reserve(num_points);

    const double log_start = std::log10(scan_lo);
    const double log_end = std::log10(scan_hi);
    if (!std::isfinite(log_start) || !std::isfinite(log_end)) {
      std::cerr << "[ScanBuild] Non-finite log bounds: " << log_start << ", " << log_end << "\n";
      return 1;
    }

    const double dlog = (num_points > 1) ? (log_end - log_start) / (num_points - 1) : 0.0;
    for (int i = 0; i < num_points; ++i) {
      central_densities.push_back(std::pow(10.0, log_start + i * dlog));
    }

    std::cout << "[ScanBuild] rho_c scan: [" << central_densities.front() << ", "
              << central_densities.back() << "] g/cm^3 over " << num_points << " points.\n";
  } else {
    std::cerr << "Failed to load EOS data. Exiting." << '\n';
    return 1;
  }

  if (log_rho.empty() || log_P.empty() || log_rho.size() != log_P.size()) {
    std::cerr << "[EoS] Bad vectors: log_rho=" << log_rho.size() << " log_P=" << log_P.size()
              << "\n";
    return 1;
  }

  // Build inverse spline: log_rho(log_P)
  gsl_interp_accel *acc_inv = gsl_interp_accel_alloc();
  gsl_spline *spline_inv = gsl_spline_alloc(gsl_interp_steffen, log_P.size());
  if (gsl_spline_init(spline_inv, log_P.data(), log_rho.data(), log_P.size()) != GSL_SUCCESS) {
    std::cerr << "Failed to init inverse spline log_rho(log_P)\n";
    return 1;
  }

  // Determine EOS validity range
  double min_logP = log_P.front();
  double max_logP = log_P.back();

  std::cout << "EOS pressure range: [" << min_logP << ", " << max_logP << "] (log10 dyne/cm²)"
            << '\n';

  // Set up GSL splines for EOS interpolation
  gsl_interp_accel *acc_eps = nullptr;
  gsl_spline *spline_eps = nullptr;
  if (!epsilon.empty() && epsilon.size() == log_P.size()) {
    acc_eps = gsl_interp_accel_alloc();
    spline_eps = gsl_spline_alloc(gsl_interp_steffen, log_P.size());
    if (gsl_spline_init(spline_eps, log_P.data(), epsilon.data(), log_P.size()) != GSL_SUCCESS) {
      std::cerr << "Failed to init epsilon(log_P); proceeding with eps ≈ rho*c^2 fallback.\n";
      gsl_spline_free(spline_eps);
      spline_eps = nullptr;
      gsl_interp_accel_free(acc_eps);
      acc_eps = nullptr;
    }
  }

  // --- choose surface and compute P_stop using inverse spline
  const double rho_surface = 1e6;
  const double goal_log_rho =
      std::clamp(std::log10(rho_surface), log_rho.front() + 1e-9, log_rho.back() - 1e-9);
  auto rho_of_logP = [&](double lp) { return gsl_spline_eval(spline_inv, lp, acc_inv); };
  double a = min_logP, b = max_logP;
  for (int it = 0; it < 80; ++it) {
    double mid = 0.5 * (a + b);
    (rho_of_logP(mid) < goal_log_rho) ? (a = mid) : (b = mid);
  }
  const double logP_stop = 0.5 * (a + b);
  const double P_stop = std::pow(10.0, logP_stop);

  // Optional sanity print:
  std::cout << "Surface target: log10(rho)=" << goal_log_rho << " -> log10(P_stop)=" << logP_stop
            << " (rho_at_stop=" << rho_of_logP(logP_stop) << ")\n";

  // Generate mass-radius relation by scanning central densities
  std::cout << "\n=== GENERATING MASS-RADIUS RELATION ===" << '\n';
  std::cout << "Computing " << central_densities.size() << " stellar models..." << '\n';

  for (double rho_c : central_densities) {
    std::ostringstream oss;
    oss << out_file << std::scientific << std::setprecision(2) << rho_c << ".csv";
    std::string filename = oss.str();

    TovResult res = non_rotating_stellar_structure_spline(
        spline_inv, acc_inv, min_logP, max_logP, rho_c, r_start, r_end, base_dlogr, true, P_stop,
        filename, spline_eps, acc_eps // <-- pass ε-spline (may be nullptr)
    );

    if (!res.found_surface) {
      std::cout << "No surface bracketed before r_end; skipping. \n";
      continue;
    }

    // Convert results to physical units
    double mass_grams = std::pow(10.0, res.log10_m_surface);
    double radius_cm = std::pow(10.0, res.log10_r_surface);
    double mass_solar = mass_grams / 1.989e33; // Convert to solar masses
    double radius_km = radius_cm / 1e5;        // Convert to kilometers

    // Report results
    std::cout << "Integration completed in " << res.steps << " steps" << '\n';
    std::cout << "Final Mass = " << mass_grams << " g = " << mass_solar << " M☉" << '\n';
    std::cout << "Final Radius = " << radius_cm << " cm = " << radius_km << " km" << '\n';
    const double twoGM_Rc2 = twoGM_over_Rc2(mass_solar, radius_km);
    const double compactness = 0.5 * twoGM_Rc2; // GM/(Rc^2)
    std::cout << "2GM/(Rc^2) = " << twoGM_Rc2 << '\n';
    std::cout << "Compactness = " << compactness << '\n';
  }

  // Clean up
  if (spline_eps) {
    gsl_spline_free(spline_eps);
  }
  if (acc_eps) {
    gsl_interp_accel_free(acc_eps);
  }
  gsl_spline_free(spline_inv);
  gsl_interp_accel_free(acc_inv);

  std::cout << "\n=== COMPUTATION COMPLETE ===" << '\n';
  std::cout << "All stellar models computed successfully!" << '\n';
  std::cout << "Output files written to data/ directory." << '\n';

  return 0;
}
