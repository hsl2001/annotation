CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra
LIBS = -lz -lm

TARGET = anno
SRCS = anno.c anno_signal.c anno_cnn.c anno_data.c anno_hmm.c genann.c
OBJS = $(SRCS:.c=.o)
HEADERS = anno_signal.h anno_cnn.h anno_data.h anno_hmm.h genann.h kseq.h
TESTS = test_signal test_cnn test_data test_hmm test_overfit

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
	$(RM) $(OBJS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

test_signal: test_signal.c anno_signal.c anno_signal.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test_signal.c anno_signal.c $(LIBS)

test_cnn: test_cnn.c anno_cnn.c anno_signal.c genann.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test_cnn.c anno_cnn.c anno_signal.c genann.c $(LIBS)

test_data: test_data.c anno_data.c anno_signal.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test_data.c anno_data.c anno_signal.c $(LIBS)

test_hmm: test_hmm.c anno_hmm.c anno_signal.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test_hmm.c anno_hmm.c anno_signal.c $(LIBS)

test_overfit: test_overfit.c $(SRCS) $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ test_overfit.c $(filter-out anno.c,$(SRCS)) $(LIBS)

test: $(TARGET) $(TESTS)
	./test_signal
	./test_cnn
	./test_data
	./test_hmm
	./test_overfit --self-test
	python3 test_pipeline.py

clean:
	rm -f $(OBJS) *.o $(TARGET) $(TESTS)