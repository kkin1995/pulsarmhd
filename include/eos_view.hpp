#pragma once
#include <cmath>
#include <optional>

struct EOSView {
  // Given log10(P), return log10(rho) and epsilon (erg/cm^3) if available.
  // If epsilon is not provided, we’ll fall back to rho*c^2.
  virtual double log10_rho_from_log10P(double log10P) const = 0;
  virtual std::optional<double> epsilon_from_log10P(double log10P) const { return std::nullopt; }

  // Clamp log10(P) to EOS validity range (no-op for polytropic).
  virtual double clamp_log10P(double log10P) const { return log10P; }

  virtual ~EOSView() = default;
};
