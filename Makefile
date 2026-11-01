CC ?= cc
CFLAGS ?= -O3
CFLAGS += -std=c11 -Wall -Wextra -pthread
LIBS = -lz -lm

.PHONY: all clean test

all: anno

anno: anno.c anno.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ anno.c $(LIBS)

test-anno: test.c anno.c anno.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ test.c $(LIBS)

test: anno test-anno
	./test-anno

clean:
	rm -f anno test-anno
