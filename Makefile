CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra -pthread
LIBS = -lz -lm

PROGRAMS = anno
TEST_PROGRAM = test-runner
HEADERS = anno.h crf.h cwt.h kseq.h ketopt.h lbfgs.h rfft.h

.PHONY: all clean test

all: $(PROGRAMS) anno_cwt

$(PROGRAMS): %: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

anno_cwt: anno.c $(HEADERS)
	$(CC) $(CFLAGS) -DCWT_BRANCH_MAPPING $(LDFLAGS) -o $@ $< $(LIBS)

$(TEST_PROGRAM): test.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

test: $(TEST_PROGRAM)
	ulimit -s unlimited && ./$(TEST_PROGRAM)

clean:
	rm -f *.o $(PROGRAMS) $(TEST_PROGRAM)