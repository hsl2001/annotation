CC ?= cc
CFLAGS ?= -O3
CFLAGS += -std=c11 -Wall -Wextra -fopenmp -Wunused-result -Wunused-function
LIBS = -lz -lm

.PHONY: all clean test

all: anno

anno: anno.c anno.h ann.c ann.h
	$(CC) $(CFLAGS) -o $@ anno.c ann.c $(LIBS)

test-anno: test.c anno.h ann.c ann.h
	$(CC) $(CFLAGS) -o $@ test.c ann.c $(LIBS)

test: anno test-anno
	./test-anno

clean:
	rm -f anno test-anno
