# PulsarMHD

A C++ library for studying compact objects using both Newtonian and General Relativistic frameworks. This library implements state-of-the-art numerical methods for modeling stellar structure with emphasis on high-density regimes and comprehensive equation of state support.

## Features

### Stellar Structure Modeling
- **Multiple hydrostatic equilibrium frameworks:**
  - Newtonian stellar structure equations
  - General Relativistic (TOV - Tolman-Oppenheimer-Volkoff) equations
- **Mass-radius relationship calculations** for white dwarfs and neutron stars
- **Surface detection algorithms** with adaptive integration

### Equation of State (EOS) Support
- **Modular Polytropic EOS System:**
  - Non-relativistic electron gas (white dwarf cores)
  - Relativistic electron gas (massive white dwarfs)
  - Non-relativistic neutron gas (low-density neutron star regions)
  - Relativistic neutron gas (high-density neutron star cores)
- **Magnetic field effects:** Magnetic BPS (Baym-Pethick-Sutherland) EOS
- **Realistic matter:** Non-magnetic neutron-proton-electron gas
- **EOS table generation** with configurable density ranges and spacing

### Numerical Methods
- **4th order Runge-Kutta integration** with adaptive step size
- **Logarithmic coordinate system** for improved numerical stability
- **Comprehensive error handling** and convergence monitoring
- **Surface detection** with pressure threshold algorithms

### Testing & Validation
- **Comprehensive test suite** (48 test cases with 100% pass rate)
- **Physics validation:** Chandrasekhar mass limits, neutron star maximum masses
- **Numerical regression testing** ensuring backward compatibility
- **Unit and integration testing** for all components

## Prerequisites
- **C++17 or higher**
- **GNU Scientific Library (GSL)** for advanced mathematical functions
- **Google Test (gtest)** for running test suite
- **Make build system**
- **Doxygen** (optional, for documentation generation)

## Installation

### Quick Start
```bash
git clone https://github.com/kkin1995/pulsarmhd.git
cd pulsarmhd
make all
```

### Building and Testing
```bash
# Build main program
make all

# Build and run comprehensive test suite
make test

# Build with debug symbols
make debug

# Clean build artifacts
make clean
```

## Usage Examples

### Basic Stellar Structure Calculation
```cpp
#include "non_rotating_stellar_structure.hpp"

// Calculate neutron star structure
auto result = non_rotating_stellar_structure(
    PolytropicGasType::NEUTRON_RELATIVISTIC,  // EOS type
    1e15,                                      // Central density (g/cm³)
    10.0,                                      // Start radius (cm)
    1e6,                                       // End radius (cm)
    0.05                                       // Step size
);

int steps = std::get<0>(result);              // Integration steps
double log_mass = std::get<1>(result);       // Final log₁₀(mass)
```

### Polytropic EOS Table Generation
```cpp
#include "polytropic_eos.hpp"

PolytropicEOS eos;
PolytropicEOSParameters params;
params.gas_type = PolytropicGasType::ELECTRON_RELATIVISTIC;
params.rho_min = 1e6;
params.rho_max = 1e12;
params.num_points = 1000;
params.output_file = "white_dwarf_eos.csv";

bool success = eos.generateEOSTable(params);
```

### Unified EOS Calculator Framework
```cpp
#include "eos_calculator.hpp"

// Using the factory pattern for different EOS types
auto calculator = EOSCalculatorFactory::createCalculator(
    EOSType::POLYTROPIC_NEUTRON_REL
);

EOSParameters params = getDefaultParams();
params.output_file = "neutron_star_eos.csv";
bool result = calculator->calculateEOS(params);
```

## Project Structure
```
pulsarmhd/
├── include/                    # Header files
│   ├── non_rotating_stellar_structure.hpp
│   ├── polytropic_eos.hpp
│   ├── eos_calculator.hpp
│   └── ...
├── src/                       # Source implementations
│   ├── non_rotating_stellar_structure.cpp
│   ├── polytropic_eos.cpp
│   ├── eos_calculator.cpp
│   └── ...
├── tests/                     # Comprehensive test suite
│   ├── non_rotating_stellar_structure_test.cpp
│   ├── polytropic_eos_test.cpp
│   └── ...
├── data/                      # Output data files
├── Makefile                   # Build system
└── README.md
```

## Documentation
- **Generate user documentation:**
```bash
make docs-user
```

- **Generate complete documentation (for maintainers):**
```bash
make docs-maintainer
```

## Testing
The project includes a comprehensive test suite validating:
- **Unit tests:** Individual function correctness
- **Integration tests:** Component interaction
- **Physics tests:** Astrophysical validation (Chandrasekhar limits, neutron star masses)
- **Regression tests:** Numerical equivalence with original implementation

Run tests with:
```bash
make test
```

## Contributing
This is a research project under active development. The codebase follows:
- **Modern C++17 standards**
- **Comprehensive documentation** with Doxygen
- **Test-driven development** with full validation
- **Clean architecture** with modular design

For collaboration inquiries, please contact the maintainers.

## Recent Updates
- ✨ **Modular Polytropic EOS System** - Comprehensive support for all degenerate gas types
- 🔧 **Refactored Stellar Structure** - Enhanced interface with preserved TOV physics
- 🧪 **Comprehensive Testing** - 48 tests ensuring reliability and correctness
- 🏗️ **Unified EOS Framework** - Factory pattern supporting multiple EOS implementations
- 📊 **Enhanced Build System** - Automatic detection and dependency management

## License
Research project - please contact maintainers for usage terms.