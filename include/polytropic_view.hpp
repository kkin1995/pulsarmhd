#pragma once
#include "eos_view.hpp"
#include "non_rotating_stellar_structure.hpp" // for c

#include <cmath>

struct PolytropicEOSView : public EOSView {
  double k_;
  double gamma_;

  PolytropicEOSView(double k, double gamma) : k_(k), gamma_(gamma) {}

  double log10_rho_from_log10P(double log10P) const override {
    return (log10P - std::log10(k_)) / gamma_;
  }
};
