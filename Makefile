CC = clang
CFLAGS = -O3 -Wall -Wextra -Xpreprocessor -fopenmp \
         -I/opt/homebrew/opt/libomp/include -Isrc
LDFLAGS = -L/opt/homebrew/opt/libomp/lib
LDLIBS = -lomp -lz -lm

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
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(APP_SOURCES) $(LIB_SOURCES) $(LDLIBS)

clean:
	rm -f surf
