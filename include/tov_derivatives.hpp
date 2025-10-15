#pragma once
#include "eos_view.hpp"
#include "non_rotating_stellar_structure.hpp" // G, c

#include <cmath>
#include <vector>

enum class GravityModel { Newtonian, RelativisticTOV };
enum class MassSource { UseRhoForMass, UseEpsilonForMass };

inline std::vector<double> tov_derivatives(double log_r, const std::vector<double> &state,
                                           const EOSView &eos, GravityModel grav,
                                           MassSource mass_src) {
  const double m = std::pow(10.0, state[0]);
  const double lp = eos.clamp_log10P(state[1]);
  const double r = std::pow(10.0, log_r);

  // Avoid divide-by-zero
  const double m_safe = (m < 1e-30 ? 1e-30 : m);

  const double log_rho = eos.log10_rho_from_log10P(lp);
  const double rho = std::pow(10.0, log_rho);
  const double P = std::pow(10.0, lp);

  // epsilon: from EOS if available, else fallback rho*c^2
  const double eps = eos.epsilon_from_log10P(lp).value_or(rho * c * c);

  double dlogm_dlogr = 0.0;
  double dlogP_dlogr = 0.0;

  if (grav == GravityModel::Newtonian) {
    // Newtonian:
    // dlogm/dlogr = (r/m) dm/dr = (r/m) * 4π r² ρ
    // dlogP/dlogr = (r/P) dP/dr = (r/P) * ( - G m ρ / r² )
    dlogm_dlogr = ((4.0 * M_PI * r * r * r * rho) / m_safe);
    dlogP_dlogr = (-(G * m_safe * rho) / (P * r));
  } else {
    // Relativistic TOV:
    // mass equation: choose ρ or ε (matches your two paths)
    const double dm_dr = (mass_src == MassSource::UseEpsilonForMass)
                             ? (4.0 * M_PI * r * r * (eps / (c * c)))
                             : (4.0 * M_PI * r * r * rho);

    dlogm_dlogr = (r / m_safe) * dm_dr;

    double first_factor = (-(G * m_safe * rho) / (P * r));
    double second_factor = (1.0 + P / (rho * c * c));
    double third_factor = 1.0 + ((4.0 * M_PI * P * r * r * r) / (m_safe * c * c));
    double fourth_factor = 1.0 / (1.0 - ((2.0 * G * m_safe) / (r * c * c)));

    double dlogP_dlogr = (first_factor * second_factor * third_factor * fourth_factor);
  }

  return {dlogm_dlogr, dlogP_dlogr};
}
