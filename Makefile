CC ?= cc
CFLAGS += -O3 -Wall -Wextra -Isrc
LDLIBS += -lz -lm

OS := $(shell uname -s)

ifeq ($(OS),Darwin)
LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null)
CFLAGS += -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
LDFLAGS += -L$(LIBOMP_PREFIX)/lib
LDLIBS += -lomp
else
CFLAGS += -fopenmp
LDLIBS += -fopenmp
endif

SRC_DIR = src
INCLUDE_DIR = $(SRC_DIR)/include

APP_SOURCES = $(SRC_DIR)/surf.c
LIB_SOURCES = $(SRC_DIR)/sequence.c $(SRC_DIR)/cwt.c \
              $(SRC_DIR)/trainer.c $(SRC_DIR)/decoder.c $(SRC_DIR)/gff.c \
              $(SRC_DIR)/intervals.c $(SRC_DIR)/lbfgs.c $(SRC_DIR)/utils.c
HEADERS = $(wildcard $(INCLUDE_DIR)/*.h)

.PHONY: all clean test itest

all: surf

test: surf
	$(CC) $(CFLAGS) $(LDFLAGS) -o surf-test test/test.c $(LIB_SOURCES) $(LDLIBS)
	./surf-test

itest: surf
	python3 test/integration_test.py

surf: $(APP_SOURCES) $(LIB_SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(APP_SOURCES) $(LIB_SOURCES) $(LDLIBS)

clean:
	rm -f surf surf-test
