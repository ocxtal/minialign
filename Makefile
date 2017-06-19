CC = gcc
GIT = git
VERSION = $(shell $(GIT) describe --tags || grep "define MM_VERSION" minialign.c | grep -o '".*"' | sed 's/"//g')
CFLAGS = -O3 -Wall -Wno-unused-function -march=native -std=c99
LDFLAGS = -lm -lz -lpthread
PREFIX = /usr/local

all: minialign

minialign: minialign.c gaba_wrap.c gaba_linear.o gaba_affine.o
	$(CC) -o $@ $(CFLAGS) -DMM_VERSION=\"$(VERSION)\" $^ $(LDFLAGS)

gaba_linear.o: gaba.c
	$(CC) -c -o $@ $(CFLAGS) -DMODEL=LINEAR -DSUFFIX $<

gaba_affine.o: gaba.c
	$(CC) -c -o $@ $(CFLAGS) -DMODEL=AFFINE -DSUFFIX $<

clean:
	rm -fr gmon.out *.o a.out minialign *~ *.a *.dSYM session*

install:
	mkdir -p $(PREFIX)/bin
	cp minialign $(PREFIX)/bin/minialign

uninstall:
	rm -f $(PREFIX)/bin/minialign

gaba.c: gaba.h log.h lmm.h unittest.h sassert.h
gaba_wrap.c: gaba.h log.h unittest.h sassert.h
minialign.c: kvec.h ksort.h gaba.h lmm.h unittest.h sassert.h
