CC := gcc
CFLAGS := -std=c11 -Wall -Wextra -Wpedantic -O2 -Iinclude
LDFLAGS := -lm

SRC := src/main.c src/solver.c
OBJ := $(SRC:.c=.o)
TARGET := sextic_solver

.PHONY: all clean test test-diverse

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(TARGET)
	printf "1 -21 175 -735 1624 -1764 720\nn\n" | ./$(TARGET)

test-diverse: $(TARGET)
	./tests/run_cases.sh ./$(TARGET)

clean:
	rm -f $(OBJ) $(TARGET)
