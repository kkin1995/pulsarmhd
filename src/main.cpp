// src/main.cpp
#include "magnetic_bps.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace fs = std::filesystem;

// --------- EDIT THESE FEW LINES TO CHANGE THE RUN ----------
struct Config {
    // Atomic mass table (CSV with header: Z,A,element,dM,BE,BDE,M,...)
    std::string masses_path = "data/atomic_masses.csv";

    // Single b value (dimensionless B/B_Q, with B_Q = 4.414e13 G)
    double b = 1e-5;

    // Baryon density grid (cm^-3), log-spaced internally by MagneticBPSEOS
    double nB_min = 1e28;
    double nB_max = 3e36;
    int    nB_pts = 250;

    // Electron solver tolerances
    double rel_tol = 1e-6;
    double abs_tol = 1e-8;

    // Output CSV
    std::string out_dir    = "data";
    std::string out_basename = "bps_single_b"; // final file: <out_dir>/<basename>_b_<sci>.csv

    // Optional: print a compact table to stdout (slower)
    bool print_table = false;
};
// -----------------------------------------------------------

static std::string tag_from_b(double b) {
    std::ostringstream os;
    os << std::scientific << std::setprecision(0) << b; // e.g., 1e-03
    std::string s = os.str();
    for (char &c : s) if (c == '+') c = '_';
    return "b_" + s;
}

int main() {
    Config cfg;

    // Create output directory if needed
    fs::create_directories(cfg.out_dir);

    std::cout << "Magnetic BPS EOS (single b)\n"
              << "Atomic masses : " << cfg.masses_path << "\n"
              << "b (B/B_Q)     : " << std::scientific << cfg.b << "\n"
              << std::defaultfloat
              << "nB grid (cm^-3): [" << cfg.nB_min << ", " << cfg.nB_max
              << "] with N=" << cfg.nB_pts << "\n"
              << "Tolerances    : rel=" << cfg.rel_tol << " abs=" << cfg.abs_tol << "\n";

    #ifdef _OPENMP
    // hard-set thread team size (overrides env for this process)
    omp_set_dynamic(0);              // keep team size fixed
    omp_set_num_threads(8);          // pick your number here
    std::cout << "[OpenMP] max_threads=" << omp_get_max_threads()
              << ", dynamic=" << omp_get_dynamic()
              << ", nested=" << omp_get_nested() << "\n";
    #endif

    // Instantiate EOS at this b
    MagneticBPSEOS eos(cfg.masses_path, /*B_ratio_electron=*/cfg.b,
                       /*rel_tol=*/cfg.rel_tol, /*abs_tol=*/cfg.abs_tol);

    // Run sweep
    auto results = eos.runSimulation(cfg.nB_min, cfg.nB_max, cfg.nB_pts);

    // Count converged/physical rows
    size_t n_ok = 0;
    for (const auto &r : results) {
        if (r.converged && r.total_pressure > 0.0) ++n_ok;
    }
    std::cout << "Converged points: " << n_ok << " / " << results.size() << "\n";

    // Write CSV
    const std::string tag = tag_from_b(cfg.b);
    fs::path out_file = fs::path(cfg.out_dir) / (cfg.out_basename + "_" + tag + ".csv");
    eos.writeEOSResults(out_file.string(), results);
    std::cout << "Wrote: " << out_file << "\n";

    if (cfg.print_table) eos.printEquilibriumTable(results);

    std::cout << "Done.\n";
    return 0;
}
