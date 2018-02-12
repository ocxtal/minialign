
# compilers
CC = gcc
GIT = git

# compiler flags
CFLAGS = -O3 -Wall -Wno-unused-function -std=c99 -pipe -DMM_VERSION=\"$(VERSION)\"
LDFLAGS = -lm -lz -lpthread

# default version string is parsed from git tags, otherwise extracted from the source
VERSION = $(shell $(GIT) describe --tags || grep "define MM_VERSION" minialign.c | grep -o '".*"' | sed 's/"//g')

# install directory and binary name
PREFIX = /usr/local
TARGET = minialign

all: native

native:
	$(MAKE) -f Makefile.core CC=$(CC) CFLAGS='$(CFLAGS)' all
	$(CC) -o $(TARGET) $(CFLAGS) minialign.o gaba.*.o $(LDFLAGS)

sse41 avx2:
	$(MAKE) -f Makefile.core CC=$(CC) CFLAGS='$(CFLAGS) -DUNITTEST=0' ARCH=`echo $@ | tr a-z A-Z` NAMESPACE=$@ all

universal: sse41 avx2
	$(CC) -o $(TARGET) $(CFLAGS) -mtune=generic universal.c minialign.*.o gaba.*.o $(LDFLAGS)

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
