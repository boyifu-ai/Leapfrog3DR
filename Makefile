CC ?= gcc
CFLAGS ?= -O2 -std=gnu99 -Wall -Wextra
LDFLAGS ?= -lm

.PHONY: all clean

all: leapfrog accuracy

leapfrog: src/leapfrog.c src/leapfrog.h src/util.h
	$(CC) $(CFLAGS) src/leapfrog.c $(LDFLAGS) -o $@

accuracy: src/accuracy.c src/util.h
	$(CC) $(CFLAGS) src/accuracy.c $(LDFLAGS) -o $@

clean:
	rm -f leapfrog accuracy leapfrog.exe accuracy.exe
