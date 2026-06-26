CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -g -Iinclude
LDFLAGS = -lm
AR = ar
ARFLAGS = rcs

SRCS = src/geomag_core.c src/geomag_model.c src/geomag_math.c src/geomag_sensor.c src/geomag_navigation.c src/geomag_kalman.c
OBJS = $(SRCS:.c=.o)
LIB = libgeomag.a

TEST_SRCS = tests/test_core.c tests/test_model.c tests/test_navigation.c tests/test_kalman.c tests/test_sensor.c
TEST_BINS = $(TEST_SRCS:.c=.exe)

EXAMPLE_SRCS = examples/example_magcom.c examples/example_igrf_position.c examples/example_compass_nav.c
EXAMPLE_BINS = $(EXAMPLE_SRCS:.c=.exe)

.PHONY: all clean test examples

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(LIB) $(TEST_BINS)
	@echo "Running tests..."
	@for t in $(TEST_BINS); do ./$$t && echo "  PASS: $$t" || echo "  FAIL: $$t"; done

tests/%.exe: tests/%.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lgeomag $(LDFLAGS) -o $@

examples: $(LIB) $(EXAMPLE_BINS)

examples/%.exe: examples/%.c $(LIB)
	$(CC) $(CFLAGS) $< -L. -lgeomag $(LDFLAGS) -o $@

clean:
	rm -f src/*.o $(LIB) tests/*.exe examples/*.exe
