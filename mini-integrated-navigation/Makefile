# Makefile for mini-integrated-navigation
# Integrated Navigation: INS/GNSS Kalman Filter Integration

CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -O2 -I./include
LDFLAGS  = -lm

SRCS     = $(wildcard src/*.c)
OBJS     = $(SRCS:.c=.o)
DEPS     = $(OBJS:.o=.d)

TEST_SRCS = $(wildcard tests/*.c)
TEST_BINS = $(TEST_SRCS:.c=.exe)

EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(EXAMPLE_SRCS:.c=.exe)

.PHONY: all clean test examples

all: libnav.a

libnav.a: $(OBJS)
	ar rcs $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@

-include $(DEPS)

test: libnav.a $(TEST_BINS)
	@for t in $(TEST_BINS); do echo "Running $$t..."; ./$$t || exit 1; done
	@echo "All tests passed."

examples: libnav.a $(EXAMPLE_BINS)
	@echo "Examples built successfully."

tests/%.exe: tests/%.c libnav.a
	$(CC) $(CFLAGS) $< -L. -lnav $(LDFLAGS) -o $@

examples/%.exe: examples/%.c libnav.a
	$(CC) $(CFLAGS) $< -L. -lnav $(LDFLAGS) -o $@

clean:
	rm -f src/*.o src/*.d tests/*.o tests/*.d tests/*.exe
	rm -f examples/*.o examples/*.d examples/*.exe
	rm -f libnav.a
