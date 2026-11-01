CC ?= cc
CFLAGS ?= -O3
CFLAGS += -std=c11 -Wall -Wextra -pthread -Wunused-result -Wunused-function -DHAVE_PTHREAD
LIBS = -lz -lm

.PHONY: all clean test

all: anno

KANN = klib/kann.c klib/kautodiff.c

anno: anno.c anno.h klib/kthread.c $(KANN)
	$(CC) $(CFLAGS) -o $@ anno.c klib/kthread.c $(KANN) $(LIBS)

test-anno: test.c anno.c anno.h klib/kthread.c $(KANN)
	$(CC) $(CFLAGS) -o $@ test.c klib/kthread.c $(KANN) $(LIBS)

test: anno test-anno
	./test-anno

clean:
	rm -f anno test-anno
