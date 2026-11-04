CC = clang
CFLAGS = -O3 -std=c11 -Wall -Wextra -Xpreprocessor -fopenmp \
         -I/opt/homebrew/opt/libomp/include
LIBS = -lz -lm -L/opt/homebrew/opt/libomp/lib -lomp

.PHONY: all clean test

all: anno

anno: anno.c anno.h unet.c
	$(CC) $(CFLAGS) -o $@ anno.c unet.c $(LIBS)

test-anno: test.c anno.c anno.h unet.c
	$(CC) $(CFLAGS) -o $@ test.c $(LIBS)

test: anno test-anno
	./test-anno

clean:
	rm -f anno test-anno
