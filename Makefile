CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra -pthread
LIBS = -lz -lm

PROGRAM = anno

.PHONY: all clean test

all: $(PROGRAM)

$(PROGRAM): anno.c anno.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ anno.c $(LIBS)

test: $(PROGRAM)
	./$(PROGRAM) --self-test

clean:
	rm -f $(PROGRAM)
