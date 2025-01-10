# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O3
DEBUG_FLAGS := -g -DDEBUG

# Directory structure
SRC_DIR := src
INC_DIR := include
EXAMPLES_DIR := examples
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin
DEP_DIR := $(BUILD_DIR)/deps
EXAMPLES_BIN_DIR := $(BUILD_DIR)/examples

# Find source files
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
HEADERS := $(wildcard $(INC_DIR)/*.hpp)
EXAMPLE_SRCS := $(wildcard $(EXAMPLES_DIR)/*.cpp)

# Object files
OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
EXAMPLE_OBJS := $(EXAMPLE_SRCS:$(EXAMPLES_DIR)/%.cpp=$(OBJ_DIR)/example_%.o)
DEPS := $(OBJS:$(OBJ_DIR)/%.o=$(DEP_DIR)/%.d)
EXAMPLE_DEPS := $(EXAMPLE_OBJS:$(OBJ_DIR)/%.o=$(DEP_DIR)/%.d)

# Main binary name
MAIN_BIN := $(BIN_DIR)/stellar_structure

# Example binaries
EXAMPLE_BINS := $(EXAMPLE_SRCS:$(EXAMPLES_DIR)/%.cpp=$(EXAMPLES_BIN_DIR)/%)

# Libraries
LIBS := -lm

# Include directories
INC_FLAGS := -I$(INC_DIR)

# Make sure the build directory and its subdirectories exist
REQUIRED_DIRS := $(sort $(BUILD_DIR) $(OBJ_DIR) $(BIN_DIR) $(DEP_DIR) $(EXAMPLES_BIN_DIR))

# Default target
all: dirs $(MAIN_BIN) $(EXAMPLE_BINS)

# Create all required directories
dirs: $(REQUIRED_DIRS)

$(REQUIRED_DIRS):
	mkdir -p $@

# Main program
$(MAIN_BIN): $(OBJS) | $(BIN_DIR)
	@echo "Linking $@..."
	$(CXX) $(OBJS) -o $@ $(LIBS)

# Example programs
$(EXAMPLES_BIN_DIR)/%: $(OBJ_DIR)/example_%.o | $(EXAMPLES_BIN_DIR)
	@echo "Linking $@..."
	$(CXX) $< -o $@ $(LIBS)

# Compile main source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR) $(DEP_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) -MMD -MP -MF $(DEP_DIR)/$*.d -c $< -o $@

# Compile example source files
$(OBJ_DIR)/example_%.o: $(EXAMPLES_DIR)/%.cpp | $(OBJ_DIR) $(DEP_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) -MMD -MP -MF $(DEP_DIR)/example_$*.d -c $< -o $@

# Include dependency files
-include $(DEPS)
-include $(EXAMPLE_DEPS)

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
	@echo "Header files:"
	@echo $(HEADERS)
	@echo "Example files:"
	@echo $(EXAMPLE_SRCS)
	@echo "Object files to be created:"
	@echo $(OBJS)
	@echo "Example object files to be created:"
	@echo $(EXAMPLE_OBJS)
	@echo "Binaries to be created:"
	@echo $(MAIN_BIN)
	@echo $(EXAMPLE_BINS)

# Run specific example (usage: make run-example EXAMPLE=free_fall)
run-example:
	@if [ -z "$(EXAMPLE)" ]; then \
		echo "Please specify an example: make run-example EXAMPLE=free_fall"; \
	else \
		$(EXAMPLES_BIN_DIR)/$(EXAMPLE); \
	fi

# Run the main program
run: $(MAIN_BIN)
	@echo "Running main program..."
	@$(MAIN_BIN)

# Run all examples
run-all-examples: $(EXAMPLE_BINS)
	@echo "Running all examples..."
	@for example in $(EXAMPLE_BINS); do \
		echo "\nRunning ${example}..."; \
		$example; \
	done

docs:
	doxygen Doxyfile

# Help
help:
	@echo "Available targets:"
	@echo "  all           - Build main program and examples"
	@echo "  debug         - Build with debug symbols"
	@echo "  clean         - Remove build files"
	@echo "  list-sources  - List all source files and targets"
	@echo "  run-example   - Run specific example (make run-example EXAMPLE=free_fall)"
	@echo "  help          - Show this help message"

.PHONY: all clean debug help list-sources dirs run-example docs