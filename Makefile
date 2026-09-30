CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra
LIBS = -lz -lm

TARGET = anno
HEADERS = anno.h kseq.h ketopt.h

.PHONY: all clean test

all: $(TARGET) anno_cwt

$(TARGET): anno.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ anno.c $(LIBS)

anno_cwt: anno_cwt.c kseq.h rfft.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ anno_cwt.c $(LIBS)

test: $(TARGET)
	uv run --no-project python3 test_exon.py

clean:
	rm -f *.o $(TARGET) anno_cwt