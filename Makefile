CC = cc
CFLAGS = -O3 -std=c11 -Wall -Wextra -pthread
LIBS = -lz -lm

NVC := $(shell command -v nvc 2>/dev/null)
ifneq ($(NVC),)
  GPU_CC = nvc
  GPU_CFLAGS = -O3 -std=c11 -acc=gpu -Minfo=accel -pthread
else
  GPU_CC = $(shell command -v gcc-acc 2>/dev/null || echo "gcc -fopenacc -foffload=nvptx-none")
  GPU_CFLAGS = -O3 -fcx-limited-range -foffload-options=nvptx-none="-fcx-limited-range" -std=c11 -Wall -Wextra -pthread
endif

CPU_DIR = build/cpu
GPU_DIR = build/gpu

PROGRAMS = anno
HEADERS = anno.h crf.h cwt.h kseq.h ketopt.h lbfgs.h rfft.h

.PHONY: all cpu gpu test clean

all: $(PROGRAMS)

cpu: $(CPU_DIR)/anno

gpu: $(GPU_DIR)/anno

$(PROGRAMS): $(CPU_DIR)/anno
	cp $< $@

$(CPU_DIR)/anno: anno.c $(HEADERS)
	@mkdir -p $(CPU_DIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

$(GPU_DIR)/anno: anno.c $(HEADERS)
	@mkdir -p $(GPU_DIR)
	$(GPU_CC) $(GPU_CFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

$(CPU_DIR)/test_compute: tests/test_compute.c anno.c $(HEADERS)
	@mkdir -p $(CPU_DIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -I. -o $@ $< $(LIBS)

$(GPU_DIR)/test_compute: tests/test_compute.c anno.c $(HEADERS)
	@mkdir -p $(GPU_DIR)
	$(GPU_CC) $(GPU_CFLAGS) $(LDFLAGS) -I. -o $@ $< $(LIBS)

test: $(CPU_DIR)/test_compute $(GPU_DIR)/test_compute
	./$(CPU_DIR)/test_compute --write $(CPU_DIR)/ref.bin
	./$(GPU_DIR)/test_compute --verify $(CPU_DIR)/ref.bin

clean:
	rm -rf *.o $(PROGRAMS) build