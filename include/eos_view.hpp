#pragma once
#include <optional>

struct EOSView {
  virtual ~EOSView() = default;

  [[nodiscard]] virtual double log10_rho_from_log10P(double log10P) const = 0;

  [[nodiscard]] virtual std::optional<double> epsilon_from_log10P(double) const {
    return std::nullopt;
  }

  [[nodiscard]] virtual double clamp_log10P(double log10P) const { return log10P; }
};
