# Makefile for mini-timing-sync
# Knowledge: Clock synchronization, PTP, NTP, Allan variance, time transfer

CC       := gcc
CFLAGS   := -Wall -Wextra -std=gnu11 -O2 -g
INCLUDES := -Iinclude
LDFLAGS  := -lm

SRC_DIR  := src
INC_DIR  := include
TST_DIR  := tests
EXA_DIR  := examples
BLD_DIR  := build

SRCS := $(SRC_DIR)/timing_sync.c \
        $(SRC_DIR)/clock_model.c \
        $(SRC_DIR)/phase_detector.c \
        $(SRC_DIR)/ptp_engine.c \
        $(SRC_DIR)/ntp_client.c \
        $(SRC_DIR)/allan_variance.c \
        $(SRC_DIR)/time_transfer.c

OBJS := $(SRCS:$(SRC_DIR)/%.c=$(BLD_DIR)/%.o)
TEST_SRC := $(TST_DIR)/test_timing_sync.c
TEST_BIN := $(BLD_DIR)/test_timing_sync
EXAMPLE_SRCS := $(wildcard $(EXA_DIR)/*.c)
EXAMPLE_BINS := $(EXAMPLE_SRCS:$(EXA_DIR)/%.c=$(BLD_DIR)/%)
HEADERS := $(wildcard $(INC_DIR)/*.h)

.PHONY: all
all: test examples

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(TEST_BIN): $(OBJS) $(TEST_SRC) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $(TEST_SRC) $(OBJS) $(LDFLAGS) -o $@

$(BLD_DIR)/%: $(EXA_DIR)/%.c $(OBJS) | $(BLD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) $< $(OBJS) $(LDFLAGS) -o $@

.PHONY: test
test: $(TEST_BIN)
	@echo "Running tests..."
	$(TEST_BIN)

.PHONY: examples
examples: $(EXAMPLE_BINS)

.PHONY: run-ptp
run-ptp: $(BLD_DIR)/example_ptp_sync
	$(BLD_DIR)/example_ptp_sync

.PHONY: run-allan
run-allan: $(BLD_DIR)/example_allan_variance
	$(BLD_DIR)/example_allan_variance

.PHONY: run-ntp
run-ntp: $(BLD_DIR)/example_ntp_time_transfer
	$(BLD_DIR)/example_ntp_time_transfer

.PHONY: run-examples
run-examples: examples
	@echo "=== Example: PTP Sync ==="
	$(BLD_DIR)/example_ptp_sync
	@echo "=== Example: Allan Variance ==="
	$(BLD_DIR)/example_allan_variance
	@echo "=== Example: NTP and Time Transfer ==="
	$(BLD_DIR)/example_ntp_time_transfer

.PHONY: lines
lines:
	@echo "include/ lines:"
	@find include -name "*.h" -exec cat {} + | wc -l
	@echo "src/ lines:"
	@find src -name "*.c" -exec cat {} + | wc -l
	@echo "src/ Lean lines:"
	@find src -name "*.lean" -exec cat {} + | wc -l

.PHONY: count
count:
	@echo "include/ + src/ total:"
	@find include src -type f | xargs cat | wc -l

.PHONY: audit
audit:
	@echo "=== Filler Detection ==="
	@grep -rn "_fn[0-9]" --include="*.c" --include="*.h" . || echo "  None found"
	@grep -rn "_aux[0-9]" --include="*.c" --include="*.h" . || echo "  None found"
	@grep -rn "sorry" --include="*.lean" . || echo "  No sorry found"
	@echo "=== Stub Detection ==="
	@find . -name "*.c" -size -200c -exec echo "SUSPECT: {}" \; || echo "  No small files"
	@echo "=== Audit Complete ==="

.PHONY: clean
clean:
	rm -rf $(BLD_DIR)

.PHONY: rebuild
rebuild: clean all

.PHONY: help
help:
	@echo "mini-timing-sync Makefile targets:"
	@echo "  make all           - Build tests and examples"
	@echo "  make test          - Build and run tests"
	@echo "  make examples      - Build all examples"
	@echo "  make run-ptp       - Run PTP sync example"
	@echo "  make run-allan     - Run Allan variance example"
	@echo "  make run-ntp       - Run NTP/time transfer example"
	@echo "  make run-examples  - Run all examples"
	@echo "  make lines         - Count include/src lines"
	@echo "  make audit         - Safety audit (filler/stub/sorry)"
	@echo "  make clean         - Remove build artifacts"
	@echo "  make rebuild       - Clean + rebuild all"
