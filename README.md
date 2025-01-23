# PulsarMHD

A C++ library for studying compact objects using both Newtonian and General Relativistic frameworks. This library implements state-of-the-art numerical methods for modeling stellar structure with emphasis on high-density regimes.

## Features
- Multiple hydrostatic equilibrium frameworks:
  - Newtonian
  - General Relativistic (TOV)
- Advanced equation of state implementations
- Numerical methods:
  - 4th order Runge-Kutta integration
  - Adaptive step size control
  - Sophisticated interpolation techniques

## Prerequisites
- C++17 or higher
- GNU Scientific Library (GSL)
- Make build system
- Doxygen (for documentation)

## Building
```bash
git clone https://github.com/yourusername/pulsarmhd.git
cd pulsarmhd
make all
```

## Documentation
- Generate user documentation:
```bash
make docs-user
```

- Generate complete documentation (for maintainers):
```bash
make docs-maintainer
```

## Contributors
This is a research project under active development. For collaboration inquiries, please contact the maintainers.