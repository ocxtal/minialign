CC = gcc
CFLAGS = -O3 -Wall -Wno-unused-function -march=native -std=c99 -D_POSIX_C_SOURCE=200112L
LIBS = -lm -lz -lpthread
SRCS = minialign.c ptask.c queue.c queue_internal.c gaba_wrap.c
PREFIX = /usr/local

all: minialign samsplit

minialign: $(SRCS) gaba_linear.o gaba_affine.o
	$(CC) -o minialign $(CFLAGS) $^ $(LIBS)

gaba_linear.o: gaba.c gaba.h unittest.h sassert.h
	$(CC) -c -o gaba_linear.o $(CFLAGS) -DMODEL=LINEAR -DSUFFIX gaba.c

gaba_affine.o: gaba.c gaba.h unittest.h sassert.h
	$(CC) -c -o gaba_affine.o $(CFLAGS) -DMODEL=AFFINE -DSUFFIX gaba.c

samsplit: samsplit.c
	$(CC) -o samsplit $(CFLAGS) samsplit.c

clean:
	rm -fr gmon.out *.o a.out minialign samsplit *~ *.a *.dSYM session*

install:
	mkdir -p $(PREFIX)/bin
	cp minialign $(PREFIX)/bin/minialign

install.all:
	mkdir -p $(PREFIX)/bin
	cp minialign $(PREFIX)/bin/minialign
	cp samsplit $(PREFIX)/bin/samsplit

uninstall:
	rm $(PREFIX)/bin/minialign $(PREFIX)/bin/samsplit

minialign.c: kvec.h ptask.h gaba.h
ptask.c: ptask.h queue.h unittest.h sassert.h
queue.c: queue.h queue_internal.h
gaba_wrap.c: gaba.h unittest.h sassert.h
