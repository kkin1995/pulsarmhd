// non_magnetic_ideal_npe_gas.hpp
#ifndef NON_MAGNETIC_IDEAL_NPE_GAS_HPP
#define NON_MAGNETIC_IDEAL_NPE_GAS_HPP

#include <iostream>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <fstream>
#include <string>
#include <vector>
#include <memory>  // for smart pointers
#include <iomanip>
#include <algorithm>

class NonMagneticNPEGas {
    private:
        // Physical Constants
        static constexpr double m_electron = 9.1094e-28;
        static constexpr double m_proton = 1.672622e-24;
        static constexpr double m_neutron = 1.674927e-24;
        static constexpr double hbar = 1.0546e-27;
        static constexpr double c = 2.9979e10;
        static constexpr double n_B_threshold = 7.37e29; // neutron appearence threshold

        // Instance Members
        bool debug_mode;
        double lambda_e;
        double lambda_p;
        double lambda_n;

        // Helper Function from Fermi-Dirac Distribution
        double phi(double x) const;
        double chi(double x) const;
        double interpolate_spline(double t, const std::vector<double>& x, const std::vector<double>& y) const;

        // Beta Equilibrium Solver
        class BetaEquilibriumSolver {
            private:
                struct Parameters {
                    double m_e;
                    double m_p;
                    double m_n;
                    double x_p;
                    double x_e;
                };

                static double f(double x_n, void* p);
                static double f_prime(double x_n, void* p);
                static void f_and_df(double x_n, void* p, double* f_value, double* df_value);
            public:
                static double solve(double x_n, double n_B, double m_e, double m_p, double m_n, double x_p, double x_e, bool debug);
        };

    public:
        // Constructor
        explicit NonMagneticNPEGas(bool debug = false);

        // Main Function to Calculate EOS
        void calculateEOS(const std::string& filename, double nB_min, double nB_max, int num_points);
};

#endif