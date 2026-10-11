CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra -pthread
LIBS = -lz -lm

PROGRAMS = anno anno_cwt
HEADERS = anno.h crf.h cwt.h kseq.h ketopt.h lbfgs.h rfft.h

.PHONY: all clean

all: $(PROGRAMS)

$(PROGRAMS): %: %.c $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

anno_cwt: anno.c

clean:
	rm -f *.o $(PROGRAMS)