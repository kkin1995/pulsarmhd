#include <cstdio>
#include <vector>
#include "rk4.hpp"

std::vector<double> calculate_derivatives(double t, const std::vector<double>& state);

const double g = 9.81; // ms-1

int main() {
    // Initial Conditions
    double y0 = 100.0;
    double v0 = 0.0;

    double t_start = 0.0;
    double dt = 0.1;
    double t_end = 10.0;

    // Initialise State Vector
    std::vector<double> state = {y0, v0};

    double t = t_start;
    while (t < t_end && state[0] > 0.0) {
        printf("t = %.2f, y = %.2f, v = %.2f\n", t, state[0], state[1]);
        if (state[0] <= 0.0) break;
        state = rk4_step(t, dt, state, calculate_derivatives);
        t += dt;
    }

    return 0;
}

std::vector<double> calculate_derivatives([[maybe_unused]] double t, const std::vector<double>& state) {
    // state[0] is position y
    // state[1] is velocity v

    return std::vector<double>{state[1], -g};
}
