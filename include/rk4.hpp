#ifndef RK4_HPP
#define RK4_HPP

#include <vector>
#include <stdexcept>
#include <functional>

/**
 * @file rk4.hpp
 * @brief Header file for the 4th-order Runge-Kutta (RK4) numerical integration method.
 * @author Karan Kinariwala
 * @date 2025-05-25
 *
 * @details
 * This file implements the classical 4th-order Runge-Kutta method for solving 
 * ordinary differential equations (ODEs). The RK4 method provides a balance between 
 * accuracy and computational efficiency in numerical simulations.
 *
 * The method solves initial value problems of the form:
 * \f{eqnarray*}{
 * \frac{dy}{dt} &=& f(t, y) \\
 * y(t_0) &=& y_0
 * \f}
 *
 * For a single step, RK4 computes:
 * \f{eqnarray*}{
 * k_1 &=& f(t_n, y_n) \\
 * k_2 &=& f(t_n + \frac{h}{2}, y_n + \frac{h}{2}k_1) \\
 * k_3 &=& f(t_n + \frac{h}{2}, y_n + \frac{h}{2}k_2) \\
 * k_4 &=& f(t_n + h, y_n + hk_3) \\
 * y_{n+1} &=& y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)
 * \f}
 *
 * where:
 * - \f$h\f$: Step size
 * - \f$t_n\f$: Current time
 * - \f$y_n\f$: Current state vector
 * - \f$k_i\f$: Intermediate derivatives
 *
 * @par Error Analysis
 * - Local truncation error: \f$O(h^5)\f$
 * - Global truncation error: \f$O(h^4)\f$
 *
 * Example usage with a simple harmonic oscillator:
* @code{.cpp}
* // System: \f$\frac{d^2x}{dt^2} = -\omega^2x\f$
* auto derivatives = [](double t, const std::vector<double>& state) 
*     -> std::vector<double> 
* {
*     const double omega = 2.0 * M_PI;  // Angular frequency
*     return {
*         state[1],                     // dx/dt = v
*         -omega * omega * state[0]     // dv/dt = -\omega^2x
*     };
* };
* @endcode
 *
 * @note The step size must be chosen carefully based on the system's characteristics
 * @warning Not optimized for stiff differential equations
 * @warning Accuracy depends on the smoothness of the derivative function
 *
 * @see non_rotating_stellar_structure.hpp for applications in stellar physics
 */

/**
 * @brief Performs one RK4 integration step to advance an ODE system.
 *
 * @details
 * Implements the classical 4th order Runge-Kutta method for solving 
 * ordinary differential equations (ODEs) of the form:
 * \f[ \frac{dy}{dt} = f(t, y) \f]
 *
 * The RK4 method computes four intermediate derivatives to advance the system:
 * \f{eqnarray*}{
 * k_1 &=& f(t_n, y_n) \\
 * k_2 &=& f(t_n + \frac{h}{2}, y_n + \frac{h}{2}k_1) \\
 * k_3 &=& f(t_n + \frac{h}{2}, y_n + \frac{h}{2}k_2) \\
 * k_4 &=& f(t_n + h, y_n + hk_3)
 * \f}
 *
 * The solution is then advanced using:
 * \f[
 * y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)
 * \f]
 *
 * @param t Current time (independent variable)
 * @param dt Time step size (must be positive)
 * @param state Current state vector of the system
 * @param derivatives Function that computes system derivatives
 * @return Updated state vector at time t + dt
 *
 * @throws std::runtime_error
 *         - If state vector is empty
 *         - If time step dt \f$\leq\f$ 0
 *         - If derivatives dimension \f$\neq\f$ state dimension
 *
 * Example using a damped harmonic oscillator:
* // System: \f$\frac{d^2x}{dt^2} + 2\zeta\omega_0\frac{dx}{dt} + \omega_0^2x = 0\f$
* auto derivatives = [](double t, const std::vector<double>& state) 
*     -> std::vector<double> 
* {
*     const double omega0 = 2.0 * M_PI;  // Natural frequency
*     const double zeta = 0.1;           // Damping ratio
*     
*     return {
*         state[1],                        // dx/dt = v
*         -2.0 * zeta * omega0 * state[1]  // dv/dt = -2\zeta\omega_0v 
*         - omega0 * omega0 * state[0]     //        -\omega_0^2x
*     };
* };
 *
 * @note Explicit, single-step method suitable for non-stiff systems
 * @warning Not efficient for stiff differential equations
 */
std::vector<double> rk4_step(double t, double dt, const std::vector<double>& state,
    const std::function<std::vector<double>(double, const std::vector<double>&)> derivatives);

#endif // RK4_HPP