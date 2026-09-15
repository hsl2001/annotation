CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra
LIBS = -lz -lm

TARGET = anno
HEADERS = anno.h kseq.h ketopt.h

.PHONY: all clean

all: $(TARGET)

$(TARGET): anno.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ anno.c $(LIBS)

clean:
	rm -f *.o $(TARGET)