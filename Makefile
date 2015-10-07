CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wno-unused-function
CPPFLAGS=
INCLUDES=	-I.
OBJS=		kthread.o misc.o bseq.o sketch.o sdust.o index.o map.o
PROG=		minimap
PROG_EXTRA=	pa-sel sdust
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

extra:all $(PROG_EXTRA)

minimap:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

sdust:sdust.c kdq.h kvec.h kseq.h sdust.h
		$(CC) -D_SDUST_MAIN $(CFLAGS) $< -o $@ -lz

pa-sel:pa-sel.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bseq.o: bseq.h kseq.h
index.o: minimap.h bseq.h kvec.h khash.h
main.o: minimap.h bseq.h
map.o: bseq.h kvec.h minimap.h sdust.h ksort.h
misc.o: minimap.h bseq.h ksort.h
pa-sel.o: kseq.h khash.h
sdust.o: kdq.h kvec.h sdust.h
sketch.o: kvec.h minimap.h bseq.h
