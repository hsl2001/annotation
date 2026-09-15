CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra
LIBS = -lz -lm

TARGET = anno
SRCS = anno.c anno_signal.c
OBJS = $(SRCS:.c=.o)
HEADERS = anno_signal.h kseq.h
TESTS = test_signal

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
	$(RM) $(OBJS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

test_signal: test_signal.c anno_signal.c anno_signal.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test_signal.c anno_signal.c $(LIBS)

test: $(TARGET) $(TESTS)
	./test_signal

clean:
	rm -f $(OBJS) *.o $(TARGET) $(TESTS)