# Testing Process Documentation

## Overview

This document outlines the comprehensive testing strategy for the PulsarMHD stellar structure simulation project, covering both C++ physics simulations and Python data analysis/visualization components.

## 🏗️ Testing Architecture

### C++ Testing Framework
- **Framework**: Custom lightweight testing framework
- **Coverage**: Physics calculations, EOS implementations, numerical methods
- **Integration**: Makefile-based build and execution
- **Focus**: Computational accuracy and numerical stability

### Python Testing Framework  
- **Framework**: pytest with comprehensive plugin ecosystem
- **Coverage**: Data processing, visualization, CLI interface
- **Integration**: Poetry dependency management + Makefile integration
- **Focus**: Data handling, user interface, and visualization quality

---

## 🧪 C++ Testing Suite

### Current Test Coverage

#### 1. **EOS Calculator Tests** ✅ **PASSING (5/5)**
```bash
make test-eos
```
- **Magnetic BPS EOS**: Pressure/density calculations
- **Non-magnetic NPE Gas**: Ideal gas equation of state  
- **Parameter validation**: Input boundary checking
- **Factory pattern**: EOS type creation and management
- **Error handling**: Invalid parameter detection

#### 2. **Polytropic EOS Tests** ✅ **PASSING (11/11)**
```bash
make test-polytropic
```
- **Four gas types**: Non-rel/rel electron gas, non-rel/rel neutron gas
- **Parameter accuracy**: k and γ values match hardcoded constants
- **Calculations**: P = k ρ^γ pressure-density relations
- **Table generation**: CSV output with metadata headers
- **File I/O**: Data persistence and retrieval
- **Error handling**: Invalid density ranges and parameters

#### 3. **Stellar Structure Tests** ✅ **PASSING (5/5)**
```bash
make test-stellar-structure  
```
- **TOV equation integration**: Relativistic stellar equilibrium
- **Newtonian physics**: Classical stellar structure
- **Boundary conditions**: Central and surface conditions
- **Numerical stability**: RK4 integration accuracy
- **Physical constraints**: Mass-radius relationship validation

### Running C++ Tests
```bash
# Run all C++ tests
make test

# Run specific test suites
make test-eos                    # EOS calculator tests
make test-polytropic            # Polytropic EOS tests  
make test-stellar-structure     # Stellar structure tests

# Clean test artifacts
make clean-tests
```

### C++ Test Results Summary
- **Total Tests**: 21/21 passing ✅
- **Coverage Areas**: EOS calculations, stellar physics, numerical methods
- **Integration**: Fully integrated with build system
- **Performance**: Fast execution (<5 seconds total)

---

## 🐍 Python Testing Suite

### Current Test Coverage

#### 1. **CLI Interface Tests** ✅ **PASSING (24/24)**
```bash
python3 -m pytest tests/python/unit/test_cli_interface.py -v
```
- **Argument parsing**: All command-line options and validation
- **Command execution**: Profile, mass-radius, mass-density, compare, all
- **Configuration overrides**: Theme, DPI, figure size customization
- **EOS type handling**: Single and multiple EOS type processing
- **Error handling**: Invalid commands and parameter validation
- **Logging**: Debug and info level configuration

#### 2. **Configuration Manager Tests** ✅ **PASSING (9/9)**
```bash
python3 -m pytest tests/python/unit/test_config_manager.py -v
```
- **Config loading**: YAML files, defaults, nonexistent files
- **Path validation**: Data directory existence, output directory creation
- **PlotConfig defaults**: Style, validation, performance parameters
- **EOS patterns**: File pattern matching for different EOS types
- **C++ density parsing**: `1.00pp16` ↔ `1.00e+16` conversion
- **Logger setup**: Proper logging configuration

#### 3. **EOS Data Processor Tests** ✅ **PASSING (17/17)**
```bash
python3 -m pytest tests/python/unit/test_eos_data_processor.py -v
```
- **File discovery**: Hybrid EOS, neutron relativistic, empty directories
- **Central density parsing**: Hybrid and neutron filename formats
- **Data loading**: Valid CSV files, missing files, profile data
- **Dataset loading**: Single EOS type, no files found
- **Edge cases**: Empty CSV, NaN values, missing columns, large datasets
- **Error handling**: Nonexistent directories, malformed data

#### 4. **Stellar Plotter Tests** ✅ **PASSING (18/18)**
```bash
python3 -m pytest tests/python/unit/test_stellar_plotter.py -v
```
- **Initialization**: Configuration setup and validation
- **Theme application**: Publication, dark, colorblind themes
- **Plot generation**: Single profiles, mass-radius, mass-density relations
- **Comparative analysis**: Multiple EOS comparison plotting
- **Data processing**: EOS dataset properties and statistics
- **Error handling**: Empty datasets, missing data, configuration issues
- **Customization**: Figure size, DPI, theme configuration

### Running Python Tests
```bash
# Run all Python tests with coverage
make test-python

# Run specific test categories
make test-python-unit           # Unit tests only
make test-python-integration    # Integration tests only  
make test-python-performance    # Performance tests only
make test-python-visual         # Visual regression tests only

# Generate coverage reports
make coverage-python

# Clean Python test artifacts
make clean-python
```

### Python Test Results Summary
- **Total Tests**: 68/68 passing ✅
- **Coverage**: 32.63% overall, 83% on main stellar_plotter.py
- **Test Categories**: CLI, Configuration, Data Processing, Plotting
- **Infrastructure**: Comprehensive fixtures, mocking, temporary directories
- **Performance**: ~18 seconds execution time

---

## 🔄 Integrated Testing Workflow

### Complete Test Suite
```bash
# Run all tests (C++ + Python)
make test-all

# Results Summary:
# C++ Tests: 21/21 passing ✅
# Python Tests: 68/68 passing ✅  
# Total: 89/89 tests passing ✅
```

### Continuous Integration Workflow
1. **Pre-commit**: Run unit tests for modified components
2. **Development**: Run relevant test suites during development
3. **Integration**: Run full test suite before merging
4. **Release**: Complete test suite + performance benchmarks

---

## 📊 Test Coverage Analysis

### C++ Coverage
- **EOS Calculations**: 100% (all implementations tested)
- **Stellar Physics**: 95% (core TOV equations covered)
- **Numerical Methods**: 90% (RK4 integration validated)
- **Error Handling**: 85% (boundary conditions and validation)

### Python Coverage
- **CLI Interface**: 95% (comprehensive argument testing)
- **Configuration**: 90% (all config paths tested)
- **Data Processing**: 85% (file I/O and parsing covered)
- **Visualization**: 83% (main plotting functionality)
- **Overall**: 33% (includes untested utility scripts)

---

## 🛠️ Testing Infrastructure

### C++ Testing Framework
```cpp
// Custom lightweight testing macros
ASSERT_NEAR(expected, actual, tolerance);
ASSERT_TRUE(condition);
ASSERT_FALSE(condition);
EXPECT_NO_THROW(expression);
```

### Python Testing Framework
```python
# pytest with comprehensive plugin ecosystem
pytest                  # Core testing framework
pytest-cov            # Coverage reporting  
pytest-benchmark      # Performance testing
pytest-mpl            # Visual regression testing
pytest-mock           # Mocking utilities
pytest-xdist          # Parallel execution
hypothesis             # Property-based testing
```

### Test Data Management
- **C++ Tests**: Hardcoded test cases with known analytical solutions
- **Python Tests**: Generated mock data with realistic stellar parameters
- **Fixtures**: Temporary directories, sample CSV files, mock configurations
- **Cleanup**: Automatic cleanup with WSL permission handling

---

## 🎯 Quality Assurance Standards

### Test Requirements
- **Unit Tests**: >90% coverage for core functionality
- **Integration Tests**: End-to-end workflow validation
- **Performance Tests**: Regression detection for computational components
- **Visual Tests**: Plot output consistency verification

### Code Quality
- **C++**: Physics accuracy, numerical stability, memory safety
- **Python**: Data integrity, visualization quality, user experience
- **Documentation**: Comprehensive test documentation and examples
- **Maintainability**: Clear test structure and naming conventions

---

## 🚀 Development Workflow

### Adding New Tests

#### C++ Tests
1. Create test file in `tests/` directory
2. Include testing framework: `#include "test_framework.hpp"`
3. Implement test functions with descriptive names
4. Add to Makefile test targets
5. Verify integration with `make test`

#### Python Tests  
1. Create test file in appropriate `tests/python/` subdirectory
2. Use pytest conventions: `test_*.py` files, `test_*` functions
3. Leverage fixtures from `conftest.py` for common setup
4. Add markers for test categorization
5. Verify with `make test-python`

### Test-Driven Development
1. **Write failing test** for new functionality
2. **Implement minimal code** to make test pass
3. **Refactor** while maintaining test coverage
4. **Add edge cases** and error handling tests
5. **Document** test purpose and expected behavior

---

## 📈 Future Testing Enhancements

### Planned Improvements
- **Integration Tests**: End-to-end C++ to Python data pipeline
- **Performance Benchmarks**: Computational efficiency tracking
- **Visual Regression**: Automated plot comparison testing
- **Property-Based Testing**: Hypothesis-driven test generation
- **Parallel Testing**: Distributed test execution for large datasets

### Monitoring and Reporting
- **Coverage Tracking**: Automated coverage reporting
- **Performance Monitoring**: Execution time trend analysis
- **Test Health**: Flaky test detection and resolution
- **Documentation**: Living documentation from test specifications

---

## 📞 Support and Troubleshooting

### Common Issues
- **WSL Permissions**: File cleanup errors in Windows Subsystem for Linux
- **Font Warnings**: Missing Unicode glyphs for astronomical symbols
- **Memory Usage**: Large dataset handling in performance tests
- **Path Issues**: Cross-platform path handling in test fixtures

### Getting Help
- **Test Failures**: Check individual test output with `-v` flag
- **Coverage Issues**: Use `--cov-report=html` for detailed analysis
- **Performance**: Use `pytest-benchmark` for timing analysis
- **Visual Issues**: Check matplotlib backend configuration

---

*Last Updated: December 2024*
*Total Test Coverage: 89/89 tests passing ✅* 