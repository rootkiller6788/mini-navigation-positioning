CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
TEST_SRC = $(wildcard tests/*.c)
TEST_BIN = $(TEST_SRC:.c=.exe)
EXAMPLE_SRC = $(wildcard examples/*.c)
EXAMPLE_BIN = $(EXAMPLE_SRC:.c=.exe)

.PHONY: all clean test examples

all: test examples

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

tests/%.exe: tests/%.c $(filter-out src/uwb_cli.o,$(OBJ))
	$(CC) $(CFLAGS) $< $(filter-out src/uwb_cli.o,$(OBJ)) -o $@ $(LDFLAGS)

examples/%.exe: examples/%.c $(filter-out src/uwb_cli.o,$(OBJ))
	$(CC) $(CFLAGS) $< $(filter-out src/uwb_cli.o,$(OBJ)) -o $@ $(LDFLAGS)

test: $(TEST_BIN)
	@for t in $(TEST_BIN); do echo "Running $$t..."; ./$$t || exit 1; done
	@echo "All tests passed."

examples: $(EXAMPLE_BIN)

clean:
	rm -f src/*.o tests/*.exe examples/*.exe
