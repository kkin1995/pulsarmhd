#include "rk4.hpp"

std::vector<double> rk4_step(double t, double dt, const std::vector<double>& state,
    std::function<std::vector<double>(double, const std::vector<double>&)> derivatives) {
    
    if (state.empty()) {
        throw std::runtime_error("State vector cannot be empty");
    }

    if (dt <= 0.0) {
        throw std::runtime_error("Time step must be positive");
    }

    std::vector<double> k1 = derivatives(t, state);
    if (k1.size() != state.size()) {
        throw std::runtime_error("Derivative function must return vector of same size as state");
    }

    std::vector<double> state_k2(state.size());;
    std::vector<double> state_k3(state.size());;
    std::vector<double> state_k4(state.size());;
    std::vector<double> final_state(state.size());

    for(size_t i = 0; i < state.size(); i++) {
        state_k2[i] = state[i] + (dt/2) * k1[i];
    }
    std::vector<double> k2 = derivatives(t + dt/2, state_k2);
    if (k2.size() != state.size()) {
        throw std::runtime_error("k2 size mismatch");
    }

    for(size_t i = 0; i < state.size(); i++) {
        state_k3[i] = state[i] + (dt/2) * k2[i];
    }
    std::vector<double> k3 = derivatives(t + dt/2, state_k3);
    if (k3.size() != state.size()) {
        throw std::runtime_error("k3 size mismatch");
    }

    for(size_t i = 0; i < state.size(); i++) {
        state_k4[i] = state[i] + dt * k3[i];
    }
    std::vector<double> k4 = derivatives(t + dt, state_k4);
    if (k4.size() != state.size()) {
        throw std::runtime_error("k4 size mismatch");
    }

    for (size_t i = 0; i < state.size(); i++) {
        final_state[i] = state[i] + (dt/6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }

    return final_state;
}