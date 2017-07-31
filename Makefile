CC = gcc
GIT = git
VERSION = $(shell $(GIT) describe --tags || grep "define MM_VERSION" minialign.c | grep -o '".*"' | sed 's/"//g')
CFLAGS = -O3 -Wall -Wno-unused-function -march=native -std=c99 -pipe -DMM_VERSION=\"$(VERSION)\"
LDFLAGS = -lm -lz -lpthread
PREFIX = /usr/local
TARGET = minialign

all: native

native:
	$(MAKE) -f Makefile.core CC=$(CC) CFLAGS='$(CFLAGS)'
	$(CC) -o $(TARGET) $(CFLAGS) minialign.o gaba_linear_16.o gaba_affine_16.o gaba_linear_32.o gaba_affine_32.o $(LDFLAGS)

sse41 avx2:
	$(MAKE) -f Makefile.core CC=$(CC) CFLAGS='$(CFLAGS) -DUNITTEST=0 -DNAMESPACE=$@' NAMESPACE=$@

universal: sse41 avx2
	$(CC) -o $(TARGET) -march=native main.c minialign.*.o gaba_linear.*.o gaba_affine.*.o $(LDFLAGS)

clean:
	rm -fr gmon.out *.o a.out $(TARGET) *~ *.a *.dSYM session*

install:
	mkdir -p $(PREFIX)/bin
	cp $(TARGET) $(PREFIX)/bin/$(TARGET)

uninstall:
	rm -f $(PREFIX)/bin/$(TARGET)

gaba.c: gaba.h log.h lmm.h unittest.h sassert.h
gaba_wrap.c: gaba.h log.h unittest.h sassert.h
minialign.c: kvec.h ksort.h gaba.h lmm.h unittest.h sassert.h
