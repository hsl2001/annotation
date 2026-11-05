CC = clang
CFLAGS = -O3 -std=c11 -Wall -Wextra -Xpreprocessor -fopenmp \
         -I/opt/homebrew/opt/libomp/include
LIBS = -lz -lm -L/opt/homebrew/opt/libomp/lib -lomp

SRC_DIR = src/

LIB_SOURCES = $(SRC_DIR)util.c $(SRC_DIR)sequence.c $(SRC_DIR)labels.c $(SRC_DIR)segments.c $(SRC_DIR)features.c $(SRC_DIR)crf.c \
              $(SRC_DIR)crf_pipeline.c $(SRC_DIR)decoder.c $(SRC_DIR)gff.c $(SRC_DIR)prediction.c $(SRC_DIR)lbfgs.c
HEADERS = $(SRC_DIR)anno.h $(SRC_DIR)anno_internal.h $(SRC_DIR)crf.h $(SRC_DIR)crf_internal.h $(SRC_DIR)lbfgs.h

.PHONY: all clean test

all: anno

anno: $(SRC_DIR)anno.c $(LIB_SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $(SRC_DIR)anno.c $(LIB_SOURCES) $(LIBS)

test-bin: test/test.c $(LIB_SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ test/test.c $(LIB_SOURCES) $(LIBS)

test: anno test-bin
	./test-bin

clean:
	rm -f anno test-bin
