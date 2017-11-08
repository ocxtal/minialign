CC = gcc
GIT = git
VERSION = $(shell $(GIT) describe --tags || grep "define MM_VERSION" minialign.c | grep -o '".*"' | sed 's/"//g')
ARCHFLAGS = -march=native
CFLAGS = -O3 -Wall -Wno-unused-function -Wno-unused-label -std=c99 -pipe -DMM_VERSION=\"$(VERSION)\"
LDFLAGS = -lm -lz -lpthread
PREFIX = /usr/local
TARGET = minialign

all: native

native:
	$(MAKE) -f Makefile.core CC=$(CC) CFLAGS='$(CFLAGS) $(ARCHFLAGS)'
	$(CC) -o $(TARGET) $(CFLAGS) minialign.o gaba.*.o $(LDFLAGS)

sse41 avx2:
	$(MAKE) -f Makefile.core CC=$(CC) CFLAGS='$(CFLAGS) $(ARCHFLAGS) -DUNITTEST=0 -DNAMESPACE=$@' NAMESPACE=$@

universal: sse41 avx2
	$(CC) -o $(TARGET) $(CFLAGS) -march=generic universal.c minialign.*.o gaba.*.o $(LDFLAGS)

clean:
	rm -fr gmon.out *.o a.out $(TARGET) *~ *.a *.dSYM session*

install:
	mkdir -p $(PREFIX)/bin
	cp $(TARGET) $(PREFIX)/bin/$(TARGET)

uninstall:
	rm -f $(PREFIX)/bin/$(TARGET)

gaba.c: gaba.h log.h unittest.h sassert.h
gaba_wrap.h: gaba.h log.h unittest.h sassert.h
minialign.c: kvec.h ksort.h gaba_wrap.h lmm.h unittest.h sassert.h
