#pragma once
#include "rk4.hpp"
#include "tov_derivatives.hpp"

#include <fstream>
#include <functional>
#include <iostream>
#include <optional>

struct IntegrateOpts {
  double r_start;
  double r_end;
  double base_dlogr;
  bool use_adaptive_stepping = true;
  double pressure_threshold = 1e-8; // dyne/cm^2
  std::optional<std::string> output_filename;
  // adaptivity clamps
  double min_step_scale = 0.2;
  double max_step_scale = 5.0;
  int progress_stride = 1000;
};

TovResult integrate_structure(const EOSView &view, GravityModel grav, MassSource mass_src,
                              [[maybe_unused]] double rho_c, double log_m0, double log_p0,
                              const IntegrateOpts &opts) {
  std::ofstream outfile;
  const bool write_output = opts.output_filename.has_value();
  if (write_output) {
    outfile.open(*opts.output_filename);
    if (outfile.is_open()) {
      outfile << "log_r[cm],log_m[g],log_P[dyne/cm^2]\n";
    }
  }

  int idx = 0;
  const double log_r_start = std::log10(opts.r_start);
  const double log_r_end = std::log10(opts.r_end);
  double log_r = log_r_start;
  double log_r_prev = log_r;
  double current_dlogr = opts.base_dlogr;
  const double dlogr_min = opts.min_step_scale * opts.base_dlogr;
  const double dlogr_max = opts.max_step_scale * opts.base_dlogr;
  const double logP_stop = std::log10(opts.pressure_threshold);

  std::vector<double> state = {log_m0, log_p0};
  std::vector<double> state_prev = state;

  while (log_r < log_r_end) {
    state = rk4_step(log_r, current_dlogr, state,
                     // capture view by reference (or value), along with grav & mass_src
                     [&view, grav, mass_src](double r, const std::vector<double> &s) {
                       return tov_derivatives(r, s, view, grav, mass_src);
                     });

    log_r += current_dlogr;

    // Surface detection at P = threshold
    const double g_prev = state_prev[1] - logP_stop;
    const double g_curr = state[1] - logP_stop;
    if (g_prev > 0.0 && g_curr <= 0.0) {
      // linear interpolation to crossing
      const double t = g_prev / (g_prev - g_curr);
      const double log_r_surf = log_r_prev + t * (log_r - log_r_prev);
      const double log_m_surf = state_prev[0] + t * (state[0] - state_prev[0]);

      if (write_output && outfile.is_open()) {
        outfile << log_r_surf << "," << log_m_surf << "," << logP_stop << "\n";
        outfile.close();
      }
      return {idx + 1, log_m_surf, log_r_surf, true};
    }

    if (write_output && outfile.is_open()) {
      outfile << log_r << "," << state[0] << "," << state[1] << "\n";
    }

    // Stability
    if (!std::isfinite(state[0]) || !std::isfinite(state[1])) {
      break;
    }

    // Progress
    if (idx % opts.progress_stride == 0) {
      std::cout << "log_r: " << log_r << ", log_m: " << state[0] << ", log_P: " << state[1]
                << ", step: " << current_dlogr << '\n';
    }

    // Adaptivity (simple slope heuristic on logP)
    if (opts.use_adaptive_stepping && idx > 0) {
      const double slope = std::abs((state[1] - state_prev[1]) / (log_r - log_r_prev));
      const double scale = (slope > 1.0 ? 0.5 : (slope < 0.1 ? 2.0 : 1.0));
      current_dlogr = std::clamp(opts.base_dlogr * scale, dlogr_min, dlogr_max);
    }

    if (log_r + current_dlogr > log_r_end) {
      current_dlogr = log_r_end - log_r;
    }

    state_prev = state;
    log_r_prev = log_r;
    ++idx;
  }

  if (write_output && outfile.is_open()) {
    outfile.close();
  }
  return {idx, state[0], log_r, false};
}
