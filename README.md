# PulsarMHD

A comprehensive C++ and Python framework for studying compact objects using both Newtonian and General Relativistic approaches. This project implements state-of-the-art numerical methods for modeling stellar structure with emphasis on high-density regimes, comprehensive equation of state support, and advanced data analysis capabilities.

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
- **Unified EOS Calculator Framework** with factory pattern design

### Data Analysis & Visualization
- **Python Analysis Suite:**
  - Advanced stellar structure plotting with multiple themes
  - Mass-radius and mass-density relationship visualization
  - Comparative analysis between different EOS types
  - Publication-ready plots with customizable styling
- **CLI Interface:** Command-line tools for data processing and plotting
- **Configuration Management:** YAML-based configuration with validation
- **Data Processing:** Efficient CSV parsing and stellar model analysis

### Numerical Methods
- **4th order Runge-Kutta integration** with adaptive step size
- **Logarithmic coordinate system** for improved numerical stability
- **Comprehensive error handling** and convergence monitoring
- **Surface detection** with pressure threshold algorithms

### Testing & Validation
- **Comprehensive test suite** (89 test cases with 100% pass rate)
  - **C++ Tests:** 21/21 passing (physics validation, numerical methods)
  - **Python Tests:** 68/68 passing (data processing, visualization, CLI)
- **Physics validation:** Chandrasekhar mass limits, neutron star maximum masses
- **Numerical regression testing** ensuring backward compatibility
- **Unit and integration testing** for all components
- **Coverage reporting:** 32.63% Python coverage with detailed analysis
- **Unified testing workflow:** Single command runs all tests

## Prerequisites

### C++ Components
- **C++17 or higher**
- **GNU Scientific Library (GSL)** for advanced mathematical functions
- **Make build system**

### Python Components
- **Python 3.8+**
- **Poetry** for dependency management
- **Core packages:** numpy, pandas, matplotlib, pyyaml
- **Testing framework:** pytest with comprehensive plugin ecosystem

### Optional
- **Doxygen** for C++ documentation generation
- **Google Test (gtest)** for extended C++ testing

## Installation

### Quick Start
```bash
git clone https://github.com/kkin1995/pulsarmhd.git
cd pulsarmhd

# Install Python dependencies
poetry install

# Build C++ components
make all
```

### Building and Testing
```bash
# Build main C++ program
make all

# Run comprehensive test suite (C++ + Python)
make test-all

# Run specific test suites
make test              # C++ tests only
make test-python       # Python tests only

# Build with debug symbols
make debug

# Generate coverage reports
make coverage-python

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

### Python Data Analysis
```bash
# Generate stellar structure plots
python3 -m scripts.stellar_plotter profile --file data/neutron_star.csv

# Create mass-radius relationships
python3 -m scripts.stellar_plotter mass-radius --eos-types neutron_relativistic

# Comparative analysis
python3 -m scripts.stellar_plotter compare --eos-types electron_rel neutron_rel

# Custom configuration
python3 -m scripts.stellar_plotter all --config custom_config.yaml --theme publication
```

## Project Structure
```
pulsarmhd/
├── include/                    # C++ Header files
│   ├── non_rotating_stellar_structure.hpp
│   ├── polytropic_eos.hpp
│   ├── eos_calculator.hpp
│   └── ...
├── src/                       # C++ Source implementations
│   ├── non_rotating_stellar_structure.cpp
│   ├── polytropic_eos.cpp
│   ├── eos_calculator.cpp
│   └── ...
├── scripts/                   # Python analysis tools
│   ├── stellar_plotter.py    # Main plotting framework
│   ├── config_manager.py     # Configuration management
│   └── ...
├── tests/                     # Comprehensive test suite
│   ├── C++ Tests (21/21 passing)
│   │   ├── non_rotating_stellar_structure_test.cpp
│   │   ├── polytropic_eos_test.cpp
│   │   └── eos_calculator_test.cpp
│   └── python/               # Python test suite (68/68 passing)
│       ├── conftest.py       # Pytest configuration
│       └── unit/             # Unit tests
│           ├── test_cli_interface.py
│           ├── test_config_manager.py
│           ├── test_eos_data_processor.py
│           └── test_stellar_plotter.py
├── data/                      # Output data files
├── configs/                   # Configuration files
├── Makefile                   # Unified build system
├── pyproject.toml            # Python dependencies
├── TESTING_PROCESS.md        # Comprehensive testing documentation
└── README.md
```

## Documentation

### C++ Documentation
```bash
# Generate user documentation
make docs-user

# Generate complete documentation (for maintainers)
make docs-maintainer
```

### Testing Documentation
- **[TESTING_PROCESS.md](TESTING_PROCESS.md)** - Comprehensive testing strategy and procedures
- **Test Coverage:** 89/89 tests passing (21 C++ + 68 Python)
- **Coverage Reports:** Automated HTML coverage reporting

## Testing

The project includes a comprehensive test suite validating:

### C++ Tests (21/21 passing)
- **EOS Calculator Tests:** Magnetic BPS, non-magnetic NPE gas, factory patterns
- **Polytropic EOS Tests:** Four gas types, calculations, table generation
- **Stellar Structure Tests:** TOV equations, Newtonian physics, numerical stability

### Python Tests (68/68 passing)
- **CLI Interface Tests:** Argument parsing, command execution, error handling
- **Configuration Tests:** YAML loading, path validation, defaults
- **Data Processing Tests:** File discovery, CSV parsing, dataset loading
- **Visualization Tests:** Plot generation, themes, comparative analysis

### Running Tests
```bash
# Run all tests (C++ + Python)
make test-all

# Run specific test suites
make test                    # C++ tests only
make test-python            # Python tests only
make test-python-unit       # Python unit tests only

# Generate coverage reports
make coverage-python
```

## Development Workflow

### Prerequisites Setup
```bash
# Install Poetry (Python dependency manager)
curl -sSL https://install.python-poetry.org | python3 -

# Install project dependencies
poetry install

# Activate virtual environment
poetry shell
```

### Testing During Development
```bash
# Quick validation
make test-python-unit

# Full validation before commit
make test-all

# Performance monitoring
python3 -m pytest tests/python/unit/ --benchmark-only
```

## Contributing
This is a research project under active development. The codebase follows:
- **Modern C++17 standards** with comprehensive testing
- **Python best practices** with Poetry dependency management
- **Test-driven development** with 89/89 tests passing
- **Clean architecture** with modular design
- **Comprehensive documentation** with detailed testing procedures

For collaboration inquiries, please contact the maintainers.

## Recent Updates
- 🧪 **Comprehensive Testing Framework** - 89 tests (21 C++ + 68 Python) with unified workflow
- 🐍 **Python Analysis Suite** - Advanced data processing and visualization capabilities
- 🎨 **Publication-Ready Plotting** - Multiple themes, customizable styling, comparative analysis
- ⚙️ **CLI Interface** - Command-line tools for data analysis and plotting
- 📊 **Coverage Reporting** - Automated test coverage with HTML reports
- 🏗️ **Unified Build System** - Single Makefile for C++ and Python components
- ✨ **Modular Polytropic EOS System** - Comprehensive support for all degenerate gas types
- 🔧 **Enhanced EOS Framework** - Factory pattern with polytropic integration
- 📋 **Testing Documentation** - Complete testing process and procedures guide

## License
Research project - please contact maintainers for usage terms.
