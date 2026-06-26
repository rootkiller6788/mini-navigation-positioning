# Makefile for mini-indoor-positioning
# Indoor Positioning and Navigation Module
#
# Targets:
#   make           - Build all examples
#   make test      - Build and run test suite
#   make clean     - Remove build artifacts
#   make all       - Build everything (examples + tests)

CC      = gcc
CFLAGS  = -std=c11 -Wall -Wextra -O2 -g
LDFLAGS = -lm -Wl,--stack,16777216

SRC_DIR   = src
INC_DIR   = include
TEST_DIR  = tests
EX_DIR    = examples
BUILD_DIR = build

# Source files
SRCS = $(SRC_DIR)/indoor_positioning.c \
       $(SRC_DIR)/fingerprint_positioning.c \
       $(SRC_DIR)/inertial_navigation.c \
       $(SRC_DIR)/sensor_fusion.c \
       $(SRC_DIR)/tof_tdoa_positioning.c \
       $(SRC_DIR)/positioning_error.c

# Example files
EXAMPLES = $(EX_DIR)/example_rssi_trilateration \
           $(EX_DIR)/example_imu_fusion \
           $(EX_DIR)/example_fingerprint_navigation

# Test file
TEST_BIN = $(BUILD_DIR)/test_positioning

.PHONY: all clean test examples

all: examples test

# Build all examples
examples: $(EXAMPLES)

$(EX_DIR)/example_rssi_trilateration: $(EX_DIR)/example_rssi_trilateration.c $(SRCS)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< $(SRCS) $(LDFLAGS) -o $@

$(EX_DIR)/example_imu_fusion: $(EX_DIR)/example_imu_fusion.c $(SRCS)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< $(SRCS) $(LDFLAGS) -o $@

$(EX_DIR)/example_fingerprint_navigation: $(EX_DIR)/example_fingerprint_navigation.c $(SRCS)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< $(SRCS) $(LDFLAGS) -o $@

# Build and run test suite
test: $(TEST_BIN)
	@echo "=== Running Indoor Positioning Test Suite ==="
	$(TEST_BIN)

$(TEST_BIN): $(TEST_DIR)/test_positioning.c $(SRCS)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) $< $(SRCS) $(LDFLAGS) -o $@

# Build test only (don't run)
build-test: $(TEST_BIN)

# Clean build artifacts
clean:
	rm -f $(EX_DIR)/example_rssi_trilateration
	rm -f $(EX_DIR)/example_imu_fusion
	rm -f $(EX_DIR)/example_fingerprint_navigation
	rm -rf $(BUILD_DIR)

# Help
help:
	@echo "Indoor Positioning Module — Build System"
	@echo ""
	@echo "Targets:"
	@echo "  make           Build all examples + test"
	@echo "  make test      Build and run test suite"
	@echo "  make examples  Build examples only"
	@echo "  make clean     Remove build artifacts"
	@echo "  make help      Show this help"
