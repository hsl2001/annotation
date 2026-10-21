CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra -pthread
LIBS = -lz -lm

PROGRAMS = anno
TEST_PROGRAM = test-runner
HEADERS = anno.h crf.h cwt.h kseq.h ketopt.h lbfgs.h rfft.h

.PHONY: all clean test

all: $(PROGRAMS) anno_cwt_general anno_cwt_intron_low_energy anno_cwt_branchpoint anno_cwt

$(PROGRAMS): %: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

anno_cwt_general: anno.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

anno_cwt_intron_low_energy: anno.c $(HEADERS)
	$(CC) $(CFLAGS) -DCWT_MAPPING_INTRON_LOW_ENERGY $(LDFLAGS) -o $@ $< $(LIBS)

anno_cwt_branchpoint: anno.c $(HEADERS)
	$(CC) $(CFLAGS) -DCWT_MAPPING_BRANCHPOINT $(LDFLAGS) -o $@ $< $(LIBS)

anno_cwt: anno_cwt_branchpoint
	cp $< $@

$(TEST_PROGRAM): test.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

test: $(TEST_PROGRAM)
	ulimit -s unlimited && ./$(TEST_PROGRAM)

clean:
	rm -f *.o $(PROGRAMS) anno_cwt anno_cwt_general anno_cwt_intron_low_energy anno_cwt_branchpoint $(TEST_PROGRAM)