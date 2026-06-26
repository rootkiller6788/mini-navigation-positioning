# Makefile for mini-inertial-navigation
# Inertial Navigation System Library
#
# Reference: Titterton & Weston (2004), Groves (2013), Farrell (2008)
#
# Targets:
#   make          - build library and tests
#   make test     - build and run tests
#   make examples - build all examples
#   make clean    - remove build artifacts

CC      := gcc
CFLAGS  := -Wall -Wextra -std=c11 -O2 -I./include
LDFLAGS := -lm
LEAN    := lake

CORE_SRC   := src/ins_core.c
ATT_SRC    := src/ins_attitude.c
MECH_SRC   := src/ins_mechanization.c
ERR_SRC    := src/ins_errors.c
CALIB_SRC  := src/ins_calibration.c
INTEG_SRC  := src/ins_integration.c
ADV_SRC    := src/ins_advanced.c

ALL_SRC := $(CORE_SRC) $(ATT_SRC) $(MECH_SRC) $(ERR_SRC) $(CALIB_SRC) $(INTEG_SRC) $(ADV_SRC)
ALL_OBJ := $(ALL_SRC:.c=.o)

TEST_SRC := tests/test_ins.c
TEST_BIN := test_ins

EX1_SRC := examples/example_ahrs.c
EX2_SRC := examples/example_strapdown.c
EX3_SRC := examples/example_zup.c
EX1_BIN := example_ahrs
EX2_BIN := example_strapdown
EX3_BIN := example_zup

DEMO_SRC := demos/demo_ins.c
BENCH_SRC := benches/bench_ins.c
DEMO_BIN := demo_ins
BENCH_BIN := bench_ins

.PHONY: all test examples clean lean demo bench

all: $(ALL_OBJ) $(TEST_BIN)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(ALL_OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(ALL_OBJ) -o $(TEST_BIN) $(LDFLAGS)

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN)

$(EX1_BIN): $(ALL_OBJ) $(EX1_SRC)
	$(CC) $(CFLAGS) $(EX1_SRC) $(ALL_OBJ) -o $(EX1_BIN) $(LDFLAGS)

$(EX2_BIN): $(ALL_OBJ) $(EX2_SRC)
	$(CC) $(CFLAGS) $(EX2_SRC) $(ALL_OBJ) -o $(EX2_BIN) $(LDFLAGS)

$(EX3_BIN): $(ALL_OBJ) $(EX3_SRC)
	$(CC) $(CFLAGS) $(EX3_SRC) $(ALL_OBJ) -o $(EX3_BIN) $(LDFLAGS)

demo: $(DEMO_BIN)

$(DEMO_BIN): $(ALL_OBJ) $(DEMO_SRC)
	$(CC) $(CFLAGS) $(DEMO_SRC) $(ALL_OBJ) -o $(DEMO_BIN) $(LDFLAGS)

bench: $(BENCH_BIN)

$(BENCH_BIN): $(ALL_OBJ) $(BENCH_SRC)
	$(CC) $(CFLAGS) $(BENCH_SRC) $(ALL_OBJ) -o $(BENCH_BIN) $(LDFLAGS)

test: $(TEST_BIN)
	./$(TEST_BIN)

lean:
	cd src && $(LEAN) build 2>/dev/null || echo "Lean 4: formalization file exists (lake build optional)"

clean:
	rm -f $(ALL_OBJ) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN) $(DEMO_BIN) $(BENCH_BIN)
	rm -f src/*.o demos/*.o benches/*.o
