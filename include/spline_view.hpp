#pragma once
#include "eos_view.hpp"
#include "non_rotating_stellar_structure.hpp" // for c

#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <optional>

struct SplineEOSView : public EOSView {
  const gsl_spline *rho_of_logP_; // inverse EOS: log10(rho) = f(log10(P))
  gsl_interp_accel *acc_rho_;
  const gsl_spline *eps_of_logP_; // optional: epsilon(log10(P))
  gsl_interp_accel *acc_eps_;
  double min_logP_;
  double max_logP_;

  SplineEOSView(const gsl_spline *inv, gsl_interp_accel *acc_inv, double minLP, double maxLP,
                const gsl_spline *eps = nullptr, gsl_interp_accel *acc_eps = nullptr)
      : rho_of_logP_(inv), acc_rho_(acc_inv), eps_of_logP_(eps), acc_eps_(acc_eps),
        min_logP_(minLP), max_logP_(maxLP) {}

  [[nodiscard]] double clamp_log10P(double lp) const override {
    return std::clamp(lp, min_logP_, max_logP_);
  }

  [[nodiscard]] double log10_rho_from_log10P(double log10P) const override {
    return gsl_spline_eval(rho_of_logP_, clamp_log10P(log10P), acc_rho_);
  }

  [[nodiscard]] std::optional<double> epsilon_from_log10P(double log10P) const override {
    if (!eps_of_logP_) {
      return std::nullopt;
    }
    return gsl_spline_eval(eps_of_logP_, clamp_log10P(log10P), acc_eps_);
  }
};
