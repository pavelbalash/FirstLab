CC := gcc
CFLAGS := -std=c11 -Wall -Wextra -Wpedantic -O2 -Iinclude
LDFLAGS := -lm

SRC := src/main.c src/solver.c
OBJ := $(SRC:.c=.o)
TARGET := sextic_solver

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(TARGET)
	./$(TARGET) --demo

clean:
	rm -f $(OBJ) $(TARGET)
