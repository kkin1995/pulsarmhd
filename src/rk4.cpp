#include "rk4.hpp"

std::vector<double> rk4_step(
    double t, double dt, const std::vector<double> &state,
    const std::function<std::vector<double>(double, const std::vector<double> &)> &derivatives) {

  if (state.empty()) {
    throw std::runtime_error("State vector cannot be empty");
  }

  if (dt <= 0.0) {
    throw std::runtime_error("Time step must be positive");
  }

  const std::size_t n = state.size();
  const double half_dt = 0.5 * dt;
  const double sixth_dt = dt / 6.0;
  std::vector<double> k1 = derivatives(t, state);
  if (k1.size() != n) {
    throw std::runtime_error("Derivative function must return vector of same size as state");
  }

  std::vector<double> state_k2(state.size());
  std::vector<double> state_k3(state.size());
  std::vector<double> state_k4(state.size());
  std::vector<double> final_state(state.size());

  for (std::size_t i = 0; i < n; ++i) {
    state_k2[i] = state[i] + half_dt * k1[i];
  }
  std::vector<double> k2 = derivatives(t + half_dt, state_k2);
  if (k2.size() != n) {
    throw std::runtime_error("k2 size mismatch");
  }

  for (size_t i = 0; i < n; ++i) {
    state_k3[i] = state[i] + half_dt * k2[i];
  }
  std::vector<double> k3 = derivatives(t + half_dt, state_k3);
  if (k3.size() != n) {
    throw std::runtime_error("k3 size mismatch");
  }

  for (std::size_t i = 0; i < n; ++i) {
    state_k4[i] = state[i] + dt * k3[i];
  }
  std::vector<double> k4 = derivatives(t + dt, state_k4);
  if (k4.size() != n) {
    throw std::runtime_error("k4 size mismatch");
  }

  for (std::size_t i = 0; i < n; ++i) {
    final_state[i] = state[i] + sixth_dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
  }

  return final_state;
}
