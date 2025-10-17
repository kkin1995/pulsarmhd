# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -Wpedantic -O3 -march=native -fopenmp
DEBUG_FLAGS := -g -DDEBUG

# Promote warnings to errors in CI only (keeps local builds flexible)
ifdef CI
CXXFLAGS += -Werror
endif


# Directory structure
SRC_DIR := src
INC_DIR := include
TEST_DIR := tests
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin
TEST_BIN_DIR := $(BUILD_DIR)/test
DEP_DIR := $(BUILD_DIR)/deps

# Python testing directories
PYTHON_TEST_DIR := tests/python
SCRIPTS_DIR := scripts

# Find source files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
TEST_SRCS := $(wildcard $(TEST_DIR)/*.cpp)
HEADERS := $(wildcard $(INC_DIR)/*.hpp)

# Object files
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
MAIN_OBJS := $(filter-out $(OBJ_DIR)/main.o $(OBJ_DIR)/main_%.o, $(OBJS))
TEST_OBJS := $(TEST_SRCS:$(TEST_DIR)/%.cpp=$(OBJ_DIR)/test_%.o)
DEPS := $(OBJS:$(OBJ_DIR)/%.o=$(DEP_DIR)/%.d)
TEST_DEPS := $(TEST_OBJS:$(OBJ_DIR)/%.o=$(DEP_DIR)/%.d)

# Main binary name
MAIN_BIN := $(BIN_DIR)/stellar_structure
TEST_BIN := $(TEST_BIN_DIR)/eos_calculator_test

# Libraries
LIBS := -lm -lgsl -lgslcblas
TEST_LIBS := $(LIBS) -lgtest -lgtest_main -pthread

# Include directories
INC_FLAGS := -I$(INC_DIR)

# Make sure the build directory and its subdirectories exist
REQUIRED_DIRS := $(sort $(BUILD_DIR) $(OBJ_DIR) $(BIN_DIR) $(TEST_BIN_DIR) $(DEP_DIR))

# Default target
all: dirs $(MAIN_BIN)

# Test target - C++ only
test: dirs $(TEST_BIN)
	@echo "Running C++ tests..."
	@$(TEST_BIN)

# Python testing targets
.PHONY: test-python test-python-unit test-python-integration test-python-performance test-python-visual
.PHONY: test-all install-python-deps coverage-python

# Install Python testing dependencies
install-python-deps:
	@echo "Installing Python testing dependencies..."
	uv lock
	uv sync --group test

# Run all Python tests
test-python: install-python-deps
	@echo "Running Python tests..."
	. .venv/bin/activate && cd $(PYTHON_TEST_DIR) && python -m pytest . -v --cov=$(SCRIPTS_DIR) --cov-report=term-missing --cov-report=html:../coverage_html

# Run Python unit tests only
test-python-unit: install-python-deps
	@echo "Running Python unit tests..."
	. .venv/bin/activate && cd $(PYTHON_TEST_DIR) && python -m pytest unit/ -v -m unit --cov=$(SCRIPTS_DIR)/stellar_plotter.py --cov-report=term-missing

# Run Python integration tests only
test-python-integration: install-python-deps
	@echo "Running Python integration tests..."
	. .venv/bin/activate && cd $(PYTHON_TEST_DIR) && python -m pytest integration/ -v -m integration --cov=$(SCRIPTS_DIR) --cov-report=term-missing

# Run Python performance tests only
test-python-performance: install-python-deps
	@echo "Running Python performance tests..."
	. .venv/bin/activate && cd $(PYTHON_TEST_DIR) && python -m pytest performance/ -v -m performance --benchmark-only

# Run Python visual regression tests only
test-python-visual: install-python-deps
	@echo "Running Python visual regression tests..."
	. .venv/bin/activate && cd $(PYTHON_TEST_DIR) && python -m pytest visual/ -v -m visual --mpl

# Run all tests (C++ and Python)
test-all: test test-python
	@echo "All tests completed successfully!"

# Run integration tests between C++ and Python
test-integration: test test-python
	@echo "Running C++ → Python integration validation..."
	cd $(PYTHON_TEST_DIR)/scripts && python validate_integration.py

# Generate Python coverage report
coverage-python: install-python-deps
	@echo "Generating Python coverage report..."
	. .venv/bin/activate && cd $(PYTHON_TEST_DIR) && python -m pytest . --cov=$(SCRIPTS_DIR) --cov-report=html:../coverage_html --cov-report=term-missing
	@echo "Coverage report generated in tests/coverage_html/"

# Create all required directories
dirs: $(REQUIRED_DIRS)

$(REQUIRED_DIRS):
	mkdir -p $@

# Main program
$(MAIN_BIN): $(OBJS) | $(BIN_DIR)
	@echo "Linking $@..."
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LIBS)

# Test program
$(TEST_BIN): $(MAIN_OBJS) $(TEST_OBJS) | $(TEST_BIN_DIR)
	@echo "Linking $@..."
	$(CXX) $(CXXFLAGS) $(MAIN_OBJS) $(TEST_OBJS) -o $@ $(TEST_LIBS)

# Compile main source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR) $(DEP_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) -MMD -MP -MF $(DEP_DIR)/$*.d -c $< -o $@

# Compile test source files
$(OBJ_DIR)/test_%.o: $(TEST_DIR)/%.cpp | $(OBJ_DIR) $(DEP_DIR)
	@echo "Compiling test $<..."
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) -MMD -MP -MF $(DEP_DIR)/test_$*.d -c $< -o $@

# Include dependency files
-include $(DEPS)
-include $(TEST_DEPS)

# Debug build
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: clean all

# Clean build files
clean:
	rm -rf $(BUILD_DIR)

# Clean Python test artifacts
clean-python:
	@echo "Cleaning Python test artifacts..."
	find $(PYTHON_TEST_DIR) -name "*.pyc" -delete
	find $(PYTHON_TEST_DIR) -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	rm -rf $(PYTHON_TEST_DIR)/.pytest_cache
	rm -rf tests/coverage_html
	rm -rf $(PYTHON_TEST_DIR)/.coverage

# Clean all artifacts
clean-all: clean clean-python

# List all source files (useful for debugging)
list-sources:
	@echo "Main source files:"
	@echo $(SRCS)
	@echo "Test source files:"
	@echo $(TEST_SRCS)
	@echo "Header files:"
	@echo $(HEADERS)
	@echo "Object files to be created:"
	@echo $(OBJS)
	@echo "Test object files to be created:"
	@echo $(TEST_OBJS)
	@echo "Binaries to be created:"
	@echo $(MAIN_BIN)
	@echo $(TEST_BIN)

# Opt-in strict warnings (keeps your defaults untouched)
STRICT_WARNINGS := -Wshadow -Wconversion -Wdouble-promotion

strict:
	@echo "Building with strict warnings..."
	$(MAKE) clean || true
	$(MAKE) CXXFLAGS="$(CXXFLAGS) $(STRICT_WARNINGS)"

asan:
	@echo "Building with AddressSanitizer..."
	$(MAKE) clean || true
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -fsanitize=address -fno-omit-frame-pointer" \
	        LIBS="$(LIBS)"

ubsan:
	@echo "Building with UndefinedBehaviorSanitizer..."
	$(MAKE) clean || true
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -fsanitize=undefined -fno-omit-frame-pointer" \
	        LIBS="$(LIBS)"


# Static analysis helpers (non-fatal). Requires compile_commands.json.
tidy:
	@if ! command -v clang-tidy >/dev/null 2>&1; then \
		echo "clang-tidy not found. Install: sudo apt install clang-tidy"; exit 1; fi
	@if [ ! -f compile_commands.json ]; then echo "Run 'make compdb' first."; exit 1; fi
	clang-tidy -p . $(SRCS) || true

tidy-all:
	@if ! command -v clang-tidy >/dev/null 2>&1; then \
		echo "clang-tidy not found. Install: sudo apt install clang-tidy"; exit 1; fi
	@if [ ! -f compile_commands.json ]; then echo "Run 'make compdb' first."; exit 1; fi
	clang-tidy -p . $(SRCS) $(TEST_SRCS) || true


cppcheck:
	@echo "Running cppcheck (non-fatal)..."
	@if [ ! -f compile_commands.json ]; then \
		echo "compile_commands.json not found. Run 'make compdb' first." ; \
		exit 1 ; \
	fi
	cppcheck --project=compile_commands.json --enable=warning,performance,portability --inline-suppr || true

# Generate compile_commands.json using bear (if installed)
compdb:
	@echo "Generating compile_commands.json with bear..."
	@command -v bear >/dev/null 2>&1 || { echo "Please install 'bear' (sudo apt install bear)"; exit 1; }
	@rm -f compile_commands.json
	bear -- $(MAKE) clean all
	@echo "compile_commands.json generated."

# Run the main program
run: $(MAIN_BIN)
	@echo "Running main program..."
	@$(MAIN_BIN)

# Help
help:
	@echo "Available targets:"
	@echo "  all                    - Build main program"
	@echo "  test                   - Build and run C++ tests"
	@echo "  test-python            - Run all Python tests with coverage"
	@echo "  test-python-unit       - Run Python unit tests only"
	@echo "  test-python-integration - Run Python integration tests only"
	@echo "  test-python-performance - Run Python performance tests only"
	@echo "  test-python-visual     - Run Python visual regression tests only"
	@echo "  test-all               - Run both C++ and Python tests"
	@echo "  test-integration       - Run full integration validation"
	@echo "  coverage-python        - Generate Python coverage report"
	@echo "  install-python-deps    - Install Python testing dependencies"
	@echo "  debug                  - Build with debug symbols"
	@echo "  clean                  - Remove C++ build files"
	@echo "  clean-python           - Remove Python test artifacts"
	@echo "  clean-all              - Remove all build and test artifacts"
	@echo "  list-sources           - List all source files and targets"
	@echo "  run                    - Run main program"
	@echo "  help                   - Show this help message"

.PHONY: all test clean debug help list-sources dirs run clean-all

# ---- FUTURE: enable recursive discovery if you add subdirs ----
# SRCS := $(shell find $(SRC_DIR) -type f -name '*.cpp' | sort)
# TEST_SRCS := $(shell find $(TEST_DIR) -type f -name '*.cpp' | sort)
#
# # Map to objects preserving subdirectory layout under $(OBJ_DIR)
# OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
# TEST_OBJS := $(patsubst $(TEST_DIR)/%.cpp,$(OBJ_DIR)/test_%.o,$(TEST_SRCS))
#
# # Put .d files next to .o to avoid name collisions
# DEPS := $(OBJS:.o=.d)
# TEST_DEPS := $(TEST_OBJS:.o=.d)
#
# # Ensure subdirs exist when compiling
# $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
# 	@echo "Compiling $<..."
# 	@mkdir -p $(dir $@)
# 	$(CXX) $(CXXFLAGS) $(INC_FLAGS) -MMD -MP -c $< -o $@
#
# # Include dependency files
# -include $(DEPS) $(TEST_DEPS)
