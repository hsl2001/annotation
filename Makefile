CC ?= clang
CFLAGS ?= -O3
CFLAGS += -std=c11 -Wall -Wextra -D_POSIX_C_SOURCE=200809L -Isrc -fopenmp
CPPFLAGS +=
LDLIBS += -lz -lm -fopenmp

SRC_DIR = src
INCLUDE_DIR = $(SRC_DIR)/include

APP_SOURCES = $(SRC_DIR)/surf.c
LIB_SOURCES = $(SRC_DIR)/sequence.c $(SRC_DIR)/cwt.c \
          $(SRC_DIR)/crf.c $(SRC_DIR)/decoder.c $(SRC_DIR)/gff.c \
          $(SRC_DIR)/lbfgs.c $(SRC_DIR)/utils.c
HEADERS = $(wildcard $(INCLUDE_DIR)/*.h)

.PHONY: all clean

all: surf

surf: $(APP_SOURCES) $(LIB_SOURCES) $(HEADERS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $(APP_SOURCES) $(LIB_SOURCES) $(LDLIBS)

clean:
	rm -f surf
