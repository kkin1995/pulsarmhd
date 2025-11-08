#include "non_rotating_stellar_structure.hpp"
#include "rk4.hpp"

#include <algorithm>
#include <cctype>
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

bool readEoSData(const std::string &filename,
                 std::vector<double> &log_rho, // OUTPUT: log10(rho [g/cm^3])
                 std::vector<double> &log_P,   // OUTPUT: log10(P   [dyn/cm^2])
                 std::vector<double> &epsilon) // OUTPUT: epsilon   [erg/cm^3] (may be empty)
{
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

  auto norm = [](std::string s) {
    s.erase(std::remove_if(
                s.begin(), s.end(),
                [](unsigned char ch) { return std::isspace(ch) || ch == '\"' || ch == '\''; }),
            s.end());
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char ch) { return std::tolower(ch); });
    return s;
  };

  // Parse header → normalized column names
  std::vector<std::string> cols;
  {
    std::stringstream sh(header);
    std::string tok;
    while (std::getline(sh, tok, ','))
      cols.push_back(norm(tok));
  }
  auto find_any = [&](std::initializer_list<const char *> names) -> int {
    for (size_t i = 0; i < cols.size(); ++i)
      for (auto *n : names)
        if (cols[i] == n)
          return (int)i;
    return -1;
  };

  // Accept many aliases
  const int i_logrho = find_any({"log_rho", "log10rho", "log10_rho"});
  const int i_logp = find_any({"log_p", "log10p", "log10_p", "logpressure", "log_pressure"});
  const int i_rho = find_any({"rho", "density", "mass_density", "total_mass_density"});
  const int i_p = find_any({"p", "pressure", "total_pressure"});
  int i_eps = find_any({"epsilon", "energy_density", "eps", "total_energy_density"});

  enum class Mode { USE_LOGS, USE_LINEAR, LEGACY_LOGNRHOP, UNKNOWN } mode = Mode::UNKNOWN;

  if (i_logrho >= 0 && i_logp >= 0) {
    mode = Mode::USE_LOGS; // new/old files w/ explicit logs
  } else if (i_rho >= 0 && i_p >= 0) {
    mode = Mode::USE_LINEAR; // new magnetic-BPS format (linear rho,P,eps)
  } else if (find_any({"log_n", "lognb"}) >= 0 && find_any({"log_rho"}) >= 0 &&
             find_any({"log_p"}) >= 0) {
    mode = Mode::LEGACY_LOGNRHOP; // very old magnetic-BPS output
  } else {
    std::cerr << "[EoS] Unrecognized header in " << filename
              << " (need (log_rho,log_P) or (rho,P)).\n";
    return false;
  }

  // Read rows
  std::string line;
  std::vector<double> tmp_lr, tmp_lp, tmp_eps;
  tmp_lr.reserve(2048);
  tmp_lp.reserve(2048);
  tmp_eps.reserve(2048);

  auto push_row = [&](double lr, double lp, double e_opt, bool have_eps) {
    if (!std::isfinite(lr) || !std::isfinite(lp))
      return;
    tmp_lr.push_back(lr);
    tmp_lp.push_back(lp);
    if (have_eps)
      tmp_eps.push_back(e_opt);
    else
      tmp_eps.push_back(std::numeric_limits<double>::quiet_NaN());
  };

  while (std::getline(file, line)) {
    if (line.empty())
      continue;

    std::vector<std::string> toks;
    {
      std::stringstream ss(line);
      std::string t;
      while (std::getline(ss, t, ','))
        toks.push_back(norm(t));
    }
    if (toks.size() < 2)
      continue;

    try {
      if (mode == Mode::USE_LOGS) {
        const double lr = std::stod(toks[(size_t)i_logrho]);
        const double lp = std::stod(toks[(size_t)i_logp]);
        double e = std::numeric_limits<double>::quiet_NaN();
        bool have = false;
        if (i_eps >= 0 && (size_t)i_eps < toks.size()) {
          e = std::stod(toks[(size_t)i_eps]);
          have = std::isfinite(e);
        }
        push_row(lr, lp, e, have);

      } else if (mode == Mode::USE_LINEAR) {
        const double rho = std::stod(toks[(size_t)i_rho]);
        const double P = std::stod(toks[(size_t)i_p]);
        if (!(rho > 0.0 && P > 0.0) || !std::isfinite(rho) || !std::isfinite(P))
          continue;
        const double lr = std::log10(rho);
        const double lp = std::log10(P);
        double e = std::numeric_limits<double>::quiet_NaN();
        bool have = false;
        if (i_eps >= 0 && (size_t)i_eps < toks.size()) {
          e = std::stod(toks[(size_t)i_eps]);
          have = std::isfinite(e) && (e > 0.0);
        }
        push_row(lr, lp, e, have);

      } else { // LEGACY_LOGNRHOP
        const int i_lrho = find_any({"log_rho"});
        const int i_lp = find_any({"log_p"});
        // We ignore log_n; we just need (log_rho, log_P)
        const double lr = std::stod(toks[(size_t)i_lrho]);
        const double lp = std::stod(toks[(size_t)i_lp]);
        push_row(lr, lp, std::numeric_limits<double>::quiet_NaN(), false);
      }
    } catch (...) {
      continue; // skip malformed row
    }
  }
  file.close();

  if (tmp_lp.empty()) {
    std::cerr << "No usable rows read from " << filename << "\n";
    return false;
  }

  // Sort by log_P and deduplicate (strictly increasing for Steffen)
  struct Row {
    double lp, lr, e;
    bool have;
  };
  std::vector<Row> rows;
  rows.reserve(tmp_lp.size());
  for (size_t i = 0; i < tmp_lp.size(); ++i) {
    rows.push_back({tmp_lp[i], tmp_lr[i], tmp_eps[i], std::isfinite(tmp_eps[i])});
  }
  std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b) { return a.lp < b.lp; });

  std::vector<Row> cleaned;
  cleaned.reserve(rows.size());
  if (!rows.empty())
    cleaned.push_back(rows.front());
  for (size_t i = 1; i < rows.size(); ++i) {
    if (rows[i].lp > cleaned.back().lp)
      cleaned.push_back(rows[i]); // strict increase
  }
  if (cleaned.size() < 2) {
    std::cerr << "[EoS] Fewer than 2 unique (log_P) rows after cleaning.\n";
    return false;
  }

  // Move to outputs
  log_P.reserve(cleaned.size());
  log_rho.reserve(cleaned.size());
  epsilon.reserve(cleaned.size());
  bool all_eps = true;

  for (const auto &r : cleaned) {
    log_P.push_back(r.lp);
    log_rho.push_back(r.lr);
    if (r.have) {
      epsilon.push_back(r.e);
    } else {
      epsilon.push_back(std::numeric_limits<double>::quiet_NaN());
      all_eps = false;
    }
  }

  // If some eps are missing, and some present → either drop eps entirely or backfill using rho*c^2.
  if (!all_eps) {
    // Drop epsilon entirely; caller will fall back to ρc² consistently.
    epsilon.clear();
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
      "sahu_basu_datta_bbp_magnetic_bps_b_1e-03.csv";

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
    // free before exit
    gsl_spline_free(spline_inv);
    gsl_interp_accel_free(acc_inv);
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
  // Validate epsilon: must be present, sized like log_P, and all finite
  auto has_bad_eps = [&]() -> bool {
    if (epsilon.empty() || epsilon.size() != log_P.size())
      return true;
    for (double e : epsilon) {
      if (!std::isfinite(e))
        return true;
    }
    return false;
  }();

  if (!has_bad_eps) {
    acc_eps = gsl_interp_accel_alloc();
    spline_eps = gsl_spline_alloc(gsl_interp_steffen, log_P.size());
    if (gsl_spline_init(spline_eps, log_P.data(), epsilon.data(), log_P.size()) != GSL_SUCCESS) {
      std::cerr << "Failed to init epsilon(log_P); proceeding with eps ≈ rho*c^2 fallback.\n";
      gsl_spline_free(spline_eps);
      spline_eps = nullptr;
      gsl_interp_accel_free(acc_eps);
      acc_eps = nullptr;
    }
  } else {
    if (!epsilon.empty()) {
      std::cerr << "epsilon column has non-finite values or mismatched size; "
                   "using eps ≈ rho*c^2 fallback.\n";
    }
    epsilon.clear(); // signal "no epsilon" to downstream code
  }

  // --- choose surface and compute P_stop using inverse spline
  const double rho_surface = 1e6;
  // clamp using actual min/max of log_rho (safer than relying on front/back)
  const auto [lr_min_it2, lr_max_it2] = std::minmax_element(log_rho.begin(), log_rho.end());
  const double lr_min2 = *lr_min_it2, lr_max2 = *lr_max_it2;
  const double goal_log_rho = std::clamp(std::log10(rho_surface), lr_min2 + 1e-9, lr_max2 - 1e-9);
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
            << " (P_stop=" << P_stop << " dyne/cm^2"
            << ", rho_at_stop=" << rho_of_logP(logP_stop) << ")\n";

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
    double mass_solar = mass_grams / kMsun_g; // Convert to solar masses
    double radius_km = radius_cm / kCmPerKm;  // Convert to kilometers

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
