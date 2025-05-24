# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O3
DEBUG_FLAGS := -g -DDEBUG

# Directory structure
SRC_DIR := src
INC_DIR := include
TEST_DIR := tests
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin
TEST_BIN_DIR := $(BUILD_DIR)/test
DEP_DIR := $(BUILD_DIR)/deps

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

# Test target
test: dirs $(TEST_BIN)
	@echo "Running tests..."
	@$(TEST_BIN)

# Create all required directories
dirs: $(REQUIRED_DIRS)

$(REQUIRED_DIRS):
	mkdir -p $@

# Main program
$(MAIN_BIN): $(OBJS) | $(BIN_DIR)
	@echo "Linking $@..."
	$(CXX) $(OBJS) -o $@ $(LIBS)

# Test program
$(TEST_BIN): $(MAIN_OBJS) $(TEST_OBJS) | $(TEST_BIN_DIR)
	@echo "Linking $@..."
	$(CXX) $(MAIN_OBJS) $(TEST_OBJS) -o $@ $(TEST_LIBS)

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

# Run the main program
run: $(MAIN_BIN)
	@echo "Running main program..."
	@$(MAIN_BIN)

# Help
help:
	@echo "Available targets:"
	@echo "  all           - Build main program"
	@echo "  test          - Build and run tests"
	@echo "  debug         - Build with debug symbols"
	@echo "  clean         - Remove build files"
	@echo "  list-sources  - List all source files and targets"
	@echo "  run           - Run main program"
	@echo "  help          - Show this help message"

.PHONY: all test clean debug help list-sources dirs run