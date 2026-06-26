# =========================================================================
# Makefile — mini-gnss-gps
# =========================================================================

CC       := gcc
CFLAGS   := -std=c99 -Wall -Wextra -O2 -g
LDFLAGS  := -lm
INC_DIR  := include
SRC_DIR  := src
TST_DIR  := tests
EX_DIR   := examples
BEN_DIR  := benches

INCLUDES := -I$(INC_DIR)

# Source files
SRCS := $(SRC_DIR)/gnss_common.c \
        $(SRC_DIR)/gnss_signal.c \
        $(SRC_DIR)/gnss_pseudorange.c \
        $(SRC_DIR)/gnss_position.c \
        $(SRC_DIR)/gnss_dop.c \
        $(SRC_DIR)/gnss_carrier.c

OBJS := $(SRCS:.c=.o)

# Targets
.PHONY: all clean test examples bench

all: libgnss.a test examples

libgnss.a: $(OBJS)
	ar rcs $@ $^

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Test target
test: tests/test_gnss
	./tests/test_gnss

tests/test_gnss: $(OBJS) $(TST_DIR)/test_gnss.c
	$(CC) $(CFLAGS) $(INCLUDES) $(TST_DIR)/test_gnss.c $(OBJS) $(LDFLAGS) -o $@

# Example targets
examples: examples/example_position_solve \
          examples/example_dop_calc \
          examples/example_carrier_smoothing

examples/example_position_solve: $(OBJS) $(EX_DIR)/example_position_solve.c
	$(CC) $(CFLAGS) $(INCLUDES) $(EX_DIR)/example_position_solve.c $(OBJS) $(LDFLAGS) -o $@

examples/example_dop_calc: $(OBJS) $(EX_DIR)/example_dop_calc.c
	$(CC) $(CFLAGS) $(INCLUDES) $(EX_DIR)/example_dop_calc.c $(OBJS) $(LDFLAGS) -o $@

examples/example_carrier_smoothing: $(OBJS) $(EX_DIR)/example_carrier_smoothing.c
	$(CC) $(CFLAGS) $(INCLUDES) $(EX_DIR)/example_carrier_smoothing.c $(OBJS) $(LDFLAGS) -o $@

# Bench target
bench: $(OBJS) $(BEN_DIR)/bench_gnss.c
	$(CC) $(CFLAGS) $(INCLUDES) $(BEN_DIR)/bench_gnss.c $(OBJS) $(LDFLAGS) -o $(BEN_DIR)/bench_gnss
	./$(BEN_DIR)/bench_gnss

# Demos target
demos: $(OBJS) demos/demo_skyplot.c
	$(CC) $(CFLAGS) $(INCLUDES) demos/demo_skyplot.c $(OBJS) $(LDFLAGS) -o demos/demo_skyplot

clean:
	rm -f $(SRC_DIR)/*.o
	rm -f $(TST_DIR)/test_gnss
	rm -f $(EX_DIR)/example_position_solve $(EX_DIR)/example_dop_calc $(EX_DIR)/example_carrier_smoothing
	rm -f $(BEN_DIR)/bench_gnss
	rm -f demos/demo_skyplot
	rm -f libgnss.a
