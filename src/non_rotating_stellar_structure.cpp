#include "non_rotating_stellar_structure.hpp"

#include "eos_view.hpp"
#include "integrate_tov.hpp"
#include "polytropic_eos.hpp"
#include "polytropic_view.hpp"
#include "spline_view.hpp"
#include "tov_derivatives.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

std::tuple<int, double, double> non_rotating_stellar_structure(PolytropicGasType eos_type,
                                                               double rho_c, double r_start,
                                                               double r_end, double dlogr,
                                                               double mu_e) {
  /*
   * STELLAR STRUCTURE INTEGRATION USING TOV EQUATIONS
   *
   * This function solves the Tolman-Oppenheimer-Volkoff (TOV) equations for relativistic
   * stellar structure in logarithmic coordinates. The approach integrates outward from
   * a small starting radius until the stellar surface is reached.
   *
   * PHYSICS:
   * - Uses polytropic EOS: P = k * ρ^γ
   * - Solves TOV equations with relativistic corrections
   * - Handles both white dwarf (electron gas) and neutron star (neutron gas) regimes
   * - Supports both non-relativistic (γ=5/3) and relativistic (γ=4/3) cases
   *
   * NUMERICAL METHOD:
   * - Logarithmic coordinates: log(r), log(m), log(P)
   * - 4th-order Runge-Kutta integration with fixed step size
   * - Surface detection via pressure threshold (P < 10^-8 dyne/cm²)
   * - Stability checks for finite values at each step
   *
   * INITIAL CONDITIONS:
   * - Assumes uniform density sphere: m₀ = (4π/3) * r_start³ * ρ_c
   * - Initial pressure from EOS: P₀ = k * ρ_c^γ
   * - Custom μₑ support for electron gases (composition effects)
   */

  // NEW: Get EOS parameters from polytropic calculator
  PolytropicEOS eos_calculator;
  auto eos_data = eos_calculator.getEOSParameters(eos_type);
  double k = eos_data.k;
  double gamma = eos_data.gamma;
  std::string name = eos_data.name;

  // Adjust k for custom mu_e (for electron gases only)
  if (eos_type == PolytropicGasType::ELECTRON_NON_RELATIVISTIC && std::abs(mu_e - 2.0) > 1e-10) {
    k = 1.0036e13 / std::pow(mu_e, 5.0 / 3.0);
  } else if (eos_type == PolytropicGasType::ELECTRON_RELATIVISTIC && std::abs(mu_e - 2.0) > 1e-10) {
    k = 1.2435e15 / std::pow(mu_e, 4.0 / 3.0);
  }

  double fraction = 4.0 / 3.0;
  double log_m0 = log10(fraction) + log10(M_PI) + 3.0 * log10(r_start) + log10(rho_c);
  double log_p0 = log10(k) + (gamma * log10(rho_c));
  std::vector<double> state = {log_m0, log_p0};

  // Adapter + integrate
  PolytropicEOSView view{k, gamma};

  IntegrateOpts opts;
  opts.r_start = r_start;
  opts.r_end = r_end;
  opts.base_dlogr = dlogr;
  opts.use_adaptive_stepping = false; // matches your fixed-step polytropic run
  opts.pressure_threshold = 1e15;     // you used 1e15 for NS surface in polytropic loop
  opts.output_filename = get_filename(name, rho_c);

  const auto res = integrate_structure(view, GravityModel::RelativisticTOV,
                                       MassSource::UseRhoForMass, rho_c, log_m0, log_p0, opts);

  // Preserve old return shape (steps, log10(m_surface), log10(r_surface))
  return {res.steps, res.log10_m_surface, res.log10_r_surface};
}

std::vector<double> newtonian(double log_r, const std::vector<double> &state, double k,
                              double gamma) {
  PolytropicEOSView view{k, gamma};
  return tov_derivatives(log_r, state, view, GravityModel::Newtonian, MassSource::UseRhoForMass);
}

std::vector<double> tolman_oppenheimer_volkoff_derivatives(double log_r,
                                                           const std::vector<double> &state,
                                                           double k, double gamma) {
  PolytropicEOSView view{k, gamma};
  return tov_derivatives(log_r, state, view, GravityModel::RelativisticTOV,
                         MassSource::UseRhoForMass);
}

std::string get_filename(const std::string &name, double rho_c) {
  std::ostringstream ss;
  ss << std::scientific << std::setprecision(2) << rho_c;
  std::string rho_str = ss.str();

  // Only replace in the scientific notation part
  std::replace(rho_str.begin(), rho_str.end(), '+', 'p');
  std::replace(rho_str.begin(), rho_str.end(), 'e', 'p');

  return "data/" + name + "_rhoc_" + rho_str + ".csv";
}

std::vector<double> tolman_oppenheimer_volkoff_derivatives_spline(
    double log_r, const std::vector<double> &state, const gsl_spline *spline_inv,
    gsl_interp_accel *acc_inv, double min_logP, double max_logP) {
  /*
   * TOV DERIVATIVES WITH SPLINE-BASED EOS
   *
   * This function computes the same TOV derivatives as the polytropic version,
   * but uses tabulated EOS data via GSL spline interpolation to get ρ(P).
   *
   * PHYSICS:
   * - Identical TOV equations to polytropic version
   * - Density obtained from spline: log_ρ = f^(-1)(log_P)
   * - Same relativistic corrections and factors
   *
   * EOS HANDLING:
   * - Pressure clamping to validity range [min_logP, max_logP]
   * - GSL spline evaluation for density lookup
   * - Error handling for extrapolation cases
   */

  SplineEOSView view{spline_inv, acc_inv, min_logP, max_logP};
  return tov_derivatives(log_r, state, view, GravityModel::RelativisticTOV,
                         MassSource::UseRhoForMass);
}

// ε-aware TOV derivatives using tabulated splines.
// - rho_of_logP: inverse EOS spline giving log10(rho) from log10(P)  [kept for diagnostics/surface
// logic]
// - eps_of_logP: spline giving epsilon (erg/cm^3) from log10(P); may be nullptr => fallback epsilon
// = rho*c^2
std::vector<double> tolman_oppenheimer_volkoff_derivatives_spline_eps(
    double log_r, const std::vector<double> &state, const gsl_spline *rho_of_logP,
    gsl_interp_accel *acc_rho, const gsl_spline *eps_of_logP, gsl_interp_accel *acc_eps,
    double min_logP, double max_logP) {

  SplineEOSView view{rho_of_logP, acc_rho, min_logP, max_logP, eps_of_logP, acc_eps};

  return tov_derivatives(log_r, state, view, GravityModel::RelativisticTOV,
                         MassSource::UseEpsilonForMass);
}

TovResult non_rotating_stellar_structure_spline(
    const gsl_spline *spline_inv, gsl_interp_accel *acc_inv, double min_logP, double max_logP,
    double rho_c, double r_start, double r_end, double base_dlogr, bool use_adaptive_stepping,
    double pressure_threshold, const std::string &output_filename, const gsl_spline *spline_eps,
    gsl_interp_accel *acc_eps) {
  /*
   * STELLAR STRUCTURE INTEGRATION WITH SPLINE-BASED EOS
   *
   * This function extends the existing TOV solver to work with tabulated EOS
   * data via GSL splines. It provides the same functionality as the polytropic
   * version but supports realistic, complex equations of state.
   *
   * PHYSICS:
   * - Same TOV equations as polytropic version
   * - Realistic EOS via spline interpolation
   * - Support for complex phase transitions and composition changes
   *
   * NUMERICAL ENHANCEMENTS:
   * - Optional adaptive step size control
   * - Configurable surface pressure threshold
   * - Optional file output (empty filename = no output)
   * - Enhanced error handling for EOS validity range
   */

  // Get forward spline for initial pressure calculation
  // We need to find the forward spline (P from ρ) - this would need to be passed
  // For now, we'll estimate initial pressure using a simple approach
  // In practice, you'd pass both forward and inverse splines

  // Initialize output file if filename provided
  // Initial mass (same as before)
  const double fraction = 4.0 / 3.0;
  const double log_m0 =
      std::log10(fraction) + std::log10(M_PI) + 3.0 * std::log10(r_start) + std::log10(rho_c);

  // Find log_p0 by bisection using inverse spline (same logic as before)
  auto rho_of_logP = [&](double lp) { return gsl_spline_eval(spline_inv, lp, acc_inv); };
  double a = min_logP, b = max_logP, target = std::log10(rho_c);
  for (int it = 0; it < 60; ++it) {
    double m = 0.5 * (a + b);
    if (rho_of_logP(m) < target)
      a = m;
    else
      b = m;
  }
  const double log_p0 = 0.5 * (a + b);

  // Adapter + integrate
  SplineEOSView view{spline_inv, acc_inv, min_logP, max_logP, spline_eps, acc_eps};

  IntegrateOpts opts;
  opts.r_start = r_start;
  opts.r_end = r_end;
  opts.base_dlogr = base_dlogr;
  opts.use_adaptive_stepping = use_adaptive_stepping;
  opts.pressure_threshold = pressure_threshold;
  if (!output_filename.empty())
    opts.output_filename = output_filename;

  return integrate_structure(view, GravityModel::RelativisticTOV, MassSource::UseEpsilonForMass,
                             rho_c, log_m0, log_p0, opts);
}
