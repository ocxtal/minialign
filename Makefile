CC = gcc
CFLAGS = -O3 -Wall -Wno-unused-function -march=native -std=c99 -D_POSIX_C_SOURCE=200112
LIBS = -lm -lz -lpthread
SRCS = minialign.c psort.c ptask.c queue.c queue_internal.c gaba_wrap.c
PREFIX = /usr/local

all: minialign

minialign: $(SRCS) gaba_linear.o gaba_affine.o
	$(CC) -o minialign $(CFLAGS) $^ $(LIBS)

gaba_linear.o:
	$(CC) $(CFLAGS) -std=c99 -DMODEL=LINEAR -DSUFFIX gaba.c -c -o gaba_linear.o

gaba_affine.o:
	$(CC) $(CFLAGS) -std=c99 -DMODEL=AFFINE -DSUFFIX gaba.c -c -o gaba_affine.o

clean:
	rm -fr gmon.out *.o a.out minialign *~ *.a *.dSYM session*

install:
	mkdir -p $(PREFIX)/bin
	cp minialign $(PREFIX)/bin/minialign

uninstall:
	rm $(PREFIX)/bin/minialign

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

minialign.c: kvec.h ptask.h psort.h gaba.h
psort.c: ptask.h psort.h unittest.h sassert.h
ptask.c: ptask.h queue.h unittest.h sassert.h
queue.c: queue.h queue_internal.h
gaba_wrap.c: gaba.h unittest.h sassert.h
