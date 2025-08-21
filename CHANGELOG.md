# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and this project adheres to [Semantic Versioning](https://semver.org/).
Before `1.0.0`, minor versions may include breaking changes if needed.

## [Unreleased]
### Added
### Changed
- Mass–density plot: apply axis cosmetics once after plotting; consistent `fig`, `ax` lifecycle (no visual change for single EOS).
### Fixed

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

[Unreleased]: https://github.com/kkin1995/pulsarmhd/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/kkin1995/pulsarmhd/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/kkin1995/pulsarmhd/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/kkin1995/pulsarmhd/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/kkin1995/pulsarmhd/compare/v0.0.1...v0.1.0
