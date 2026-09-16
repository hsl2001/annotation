CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra
LIBS = -lz -lm

TARGET = anno
HEADERS = anno.h kseq.h ketopt.h

.PHONY: all clean test

all: $(TARGET)

$(TARGET): anno.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ anno.c $(LIBS)

test: $(TARGET)
	uv run --no-project python3 test_exon.py

clean:
	rm -f *.o $(TARGET)