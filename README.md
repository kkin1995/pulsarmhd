# pulsarmhd

A C++ library for modelling compact objects such as neutron stars and incorporating the effects of rotation and magnetic fields using both Newtonian and General Relativistic frameworks.

## Current Features
- Newtonian and Tolman-Oppenheimer-Volkoff (TOV) equations of hydrostatic equilibrium.
- Support for various equations of state:
    - Non-relativistic electron/neutron gas (polytrope).
    - Ultra-relativistic electron/neutron gas (polytrope).
- Numerical Methods:
    - Runge-Kutta 4th Order (RK4) integration.
    - Adaptive step size control.
- Analysis Tools:
    - Python plotting utilities.

## Building
```bash
git clone https://github.com/yourusername/pulsarmhd.git
cd pulsarmhd
make all
make run
```

## Generating Documentation
```bash
doxygen Doxyfile
```