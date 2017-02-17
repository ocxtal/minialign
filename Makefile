CC = gcc
CFLAGS = -O3 -Wall -Wno-unused-function -march=native -std=c99
LIBS = -lm -lz -lpthread
SRCS = minialign.c gaba_wrap.c
PREFIX = /usr/local

all: minialign

minialign: $(SRCS) gaba_linear.o gaba_affine.o
	$(CC) -o minialign $(CFLAGS) $^ $(LIBS)

gaba_linear.o: gaba.c gaba.h unittest.h sassert.h
	$(CC) -c -o gaba_linear.o $(CFLAGS) -DMODEL=LINEAR -DSUFFIX gaba.c

gaba_affine.o: gaba.c gaba.h unittest.h sassert.h
	$(CC) -c -o gaba_affine.o $(CFLAGS) -DMODEL=AFFINE -DSUFFIX gaba.c

clean:
	rm -fr gmon.out *.o a.out minialign *~ *.a *.dSYM session*

install:
	mkdir -p $(PREFIX)/bin
	cp minialign $(PREFIX)/bin/minialign

uninstall:
	rm -f $(PREFIX)/bin/minialign

minialign.c: kvec.h kseq.h ksort.h gaba.h lmm.h sassert.h
gaba_wrap.c: gaba.h unittest.h sassert.h

