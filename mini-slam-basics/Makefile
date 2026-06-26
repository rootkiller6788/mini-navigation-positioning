# Makefile for mini-slam-basics
# Simultaneous Localization and Mapping (SLAM) Basics
#
# Reference:
#   Thrun, Burgard & Fox (2005) "Probabilistic Robotics"
#   Dissanayake et al. (2001) "A Solution to the SLAM Problem"
#   Montemerlo et al. (2002) "FastSLAM"
#   Dellaert & Kaess (2006) "Square Root SAM"
#
# Targets:
#   make          - build library and tests
#   make test     - build and run tests
#   make examples - build all examples
#   make clean    - remove build artifacts

CC      := gcc
CFLAGS  := -Wall -Wextra -std=c11 -O2 -I./include -D_POSIX_C_SOURCE=200809L
LDFLAGS := -lm

# Source files
SRC := src/slam_core.c       \
       src/slam_ekf.c        \
       src/slam_fastslam.c   \
       src/slam_graph.c      \
       src/slam_sensor.c     \
       src/slam_data_assoc.c \
       src/slam_loop_closure.c

# Object files
OBJ := $(SRC:.c=.o)

# Test
TEST_SRC := tests/test_slam.c
TEST_BIN := test_slam

# Examples
EX1_SRC := examples/example_ekf_slam.c
EX2_SRC := examples/example_graph_slam.c
EX3_SRC := examples/example_fastslam.c
EX1_BIN := example_ekf_slam
EX2_BIN := example_graph_slam
EX3_BIN := example_fastslam

.PHONY: all test examples clean

all: $(OBJ) $(TEST_BIN)

# Compile C source to object
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Link test binary
$(TEST_BIN): $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(OBJ) -o $(TEST_BIN) $(LDFLAGS)

# Build examples
examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

$(EX1_BIN): $(OBJ) $(EX1_SRC)
	$(CC) $(CFLAGS) $(EX1_SRC) $(OBJ) -o $(EX1_BIN) $(LDFLAGS)

$(EX2_BIN): $(OBJ) $(EX2_SRC)
	$(CC) $(CFLAGS) $(EX2_SRC) $(OBJ) -o $(EX2_BIN) $(LDFLAGS)

$(EX3_BIN): $(OBJ) $(EX3_SRC)
	$(CC) $(CFLAGS) $(EX3_SRC) $(OBJ) -o $(EX3_BIN) $(LDFLAGS)

# Run tests
test: $(TEST_BIN)
	./$(TEST_BIN)

# Clean
clean:
	rm -f $(OBJ) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)
	rm -f src/*.o tests/*.o examples/*.o
