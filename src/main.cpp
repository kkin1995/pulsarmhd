#include <iostream>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

struct parameters {
    double m_e;
    double m_p;
    double m_n;
    double lambda_e;
    double lambda_p;
    double lambda_n;
    double n_B;
};

double f(double x, void* p) {
    struct parameters* params = (struct parameters*) p;
    double first_term = params->m_e * pow(1.0 + pow(x, 2.0) * ( pow(params->lambda_e, 2.0) / pow(params->lambda_p, 2.0) ), 1.0/2.0);
    double second_term = params->m_p * pow(1.0 + pow(x, 2.0), 1.0 / 2.0);
    double third_term = params->m_n * pow((1.0 + pow(pow(params->lambda_n, 3.0) * params->n_B - ( (pow(params->lambda_n, 3.0) * pow(x, 3.0)) / (pow(params->lambda_p, 3.0)) ), 2.0 / 3.0)), 1.0 / 2.0);

    double result = first_term + second_term - third_term;

    // Debugging output
    std::cout << "x: " << x << ", first_term: " << first_term << ", second_term: " << second_term << ", third_term: " << third_term << ", result: " << result << std::endl;

    return result;
}

int main() {
    double m_p = 1.6726e-24; // grams
    double m_n = 1.6750e-24; // grams
    double m_e = 9.1094e-28; // grams
    double hbar = 1.0546e-27; // ergs * seconds
    double c = 2.9979e10; // cm / s
    int status;
    int iter = 0, max_iter = 100;

    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    double r = 0.0;

    double lambda_p =  hbar / (m_p * c);
    double lambda_n =  hbar / (m_n * c);
    double lambda_e =  hbar / (m_e * c);

    double n_B = 7.37e30; // baryons / cm^3 - Neutron Appearence Density

    double x_p_initial = pow(3 * pow(M_PI, 2.0) * pow(lambda_p, 3.0) * n_B, 1.0/3.0);

    double x_lo = 0.0;
    double x_hi = x_p_initial * 2.0;

    parameters params = {
        .m_e = m_e,
        .m_p = m_p,
        .m_n = m_n,
        .lambda_e = lambda_e,
        .lambda_p = lambda_p,
        .lambda_n = lambda_n,
        .n_B = n_B
   };

    gsl_function F;
    F.function = &f;
    F.params = &params;

    if ((status = gsl_root_fsolver_set(s, &F, x_lo, x_hi))) {
        std::cout << "Error in solver initialization" << std::endl;
        return 1;
    }
    
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        r = gsl_root_fsolver_root(s);
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lo, x_hi, 0, 1e-6);
    } while (status == GSL_CONTINUE && iter < max_iter);

    if (status == GSL_SUCCESS) {
        std::cout << "Found solution x_p = " << r << std::endl;
        std::cout << "f(x_p) = " << f(r, &params) << std::endl;
    } else {
        std::cout << "Failed to converge after " << iter << " iterations" << std::endl;
        std::cout << "Current bounds: [" << x_lo << ", " << x_hi << "]" << std::endl;
    }
    
    gsl_root_fsolver_free(s);

    std::cout << "Results after " << iter << " iterations:" << std::endl;
    std::cout << "x_p = " << r << std::endl;
    return 0;
}