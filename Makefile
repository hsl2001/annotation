CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra
LIBS = -lz -lm

TARGET = anno
HEADERS = anno.h kseq.h
TESTS = test_signal

.PHONY: all clean test

all: $(TARGET)

$(TARGET): anno.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ anno.c $(LIBS)

test_signal: test_signal.c anno.c anno.h
	$(CC) $(CFLAGS) $(LDFLAGS) -DANNO_NO_MAIN -o $@ test_signal.c anno.c $(LIBS)

test: $(TARGET) $(TESTS)
	./test_signal

clean:
	rm -f *.o $(TARGET) $(TESTS)