# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and this project adheres to [Semantic Versioning](https://semver.org/).
Before `1.0.0`, minor versions may include breaking changes if needed.

## [Unreleased]
### Added
### Changed
### Fixed

## [0.5.0] - 2025-10-15
### Added
- **C++20 toolchain**: use `<numbers>` for π with a portable fallback; remove reliance on `M_PI`.
- **Unified TOV core**:
  - `EOSView` abstraction + concrete `PolytropicEOSView` and `SplineEOSView`.
  - Single `tov_derivatives(...)` path with `GravityModel` (Newtonian/RelativisticTOV)
    and `MassSource` (ρ-based or ε-based) switches.
  - `integrate_structure(...)` driver centralizing RK4 stepping, surface detection, and I/O.
- **ε-aware path in unified core**: optional `epsilon(log10P)` spline, with fallback `ε ≈ ρ c²`.
- **Surface detection improvements**: bracket via inverse EOS ρ(P); linear crossing interpolation.
- **Stricter code hygiene**: `[[nodiscard]]` on EOS view queries; pre-commit formatting tidy-ups.

### Changed
- Mass–density plot: apply axis cosmetics once after plotting; consistent `fig`, `ax` lifecycle (no visual change for single EOS).
- **`main.cpp`** now routes all models through the unified TOV core; prints compactness
  in dimensionless form \(GM/(Rc^2)\) (same numeric factor as before).
- **File I/O**: consistent CSV headers and filename helper; fewer spurious logs.
- **Adaptive stepping**: gentler step-size scaling tied to \|d logP / d log r\| (smoother traces).

### Fixed
- Clang-tidy complaints: braces-around-statements, else-after-return, shadowed locals.
- Corner cases in EOS cleaning (strictly monotone `log_P` with ε kept aligned).
- Intermittent “no surface bracketed” when `ρ_surface` near table boundaries.

### Removed
- Local duplicates of TOV derivative functions in favor of the unified core.
- Ad-hoc `M_PI` macro defines (now using C++20 `<numbers>` or shim).

## [0.4.0] - 2025-08-21
### Added
- **Publication-ready mass–radius (M–R) plotting** in `scripts/stellar_plotter.py`:
  optional line connections, optional `M_max` marker/label, minor ticks, inward
  ticks, grid control, and PNG + PDF export.
- **Robust configuration + theme controls** (publication/dark/presentation/colorblind),
  curve width/alpha, legend/title toggles, marker sizes, and radius capping.
- **Paper recipe** `scripts/paper_mr.yaml` (Fig. 1).
### Changed
- Title/legend now **optional via style flags**; defaults tuned for paper figures.
- Cleaner logging and validation messages.
### Fixed
- Guarded against unphysical 1000 km end-radius artifacts in plotting path.

## [0.3.0] - 2025-08-20
### Added
- **Energy-density aware TOV integration**: support for ε(P) and using ε in dP/dr.
- **Robust EOS CSV loader**: optional ε column, building ρ(P) and ε(P) splines.
- **`mr_extractor` utility** for extracting M–R data.
### Changed
- M–R plotting: scatter-only option and artifact skipping at ~1000 km.
- Repo hygiene: moved parked mains to `examples/`, expanded `.gitignore`.
### Fixed
- Updated SciPy integration call; improved EOS join diagnostics.

## [0.2.0] - 2025-05-26
### Added
- **Spline-based TOV solver** for realistic tabulated EOS.
- **Comprehensive test suite** for spline/TOV functionality.
### Changed
- Modernized `main.cpp` and related build glue; compatibility fixes.

## [0.1.0] - 2025-05-24
### Added
- **Modular polytropic EOS** and unified EOS calculator framework.
- **Unified plotting engine with CLI**, themes, validation, and perf options.
- Initial documentation and build updates; early tests.

## [0.0.1] - 2025-02-22
### Added
- Doxygen documentation for magnetic BPS EOS (pre-alpha groundwork).

---

[Unreleased]: https://github.com/kkin1995/pulsarmhd/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/kkin1995/pulsarmhd/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/kkin1995/pulsarmhd/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/kkin1995/pulsarmhd/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/kkin1995/pulsarmhd/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/kkin1995/pulsarmhd/compare/v0.0.1...v0.1.0
