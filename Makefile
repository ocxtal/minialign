CC=			gcc
CFLAGS=		-g -Wall -O2 -Wc++-compat -Wno-unused-function
CPPFLAGS=
INCLUDES=	
OBJS=		kthread.o misc.o bseq.o sketch.o index.o map.o
PROG=		minimap
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

minimap:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bseq.o: bseq.h kseq.h
index.o: minimap.h bseq.h kvec.h khash.h ksort.h
main.o: minimap.h bseq.h
map.o: bseq.h kvec.h minimap.h ksort.h
misc.o: minimap.h bseq.h ksort.h
sketch.o: kvec.h minimap.h bseq.h
