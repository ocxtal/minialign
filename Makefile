CC=			gcc
CFLAGS=		-g -Wall -O2 -Wno-unused-function #-fno-inline-functions -fno-inline-functions-called-once
CPPFLAGS=
INCLUDES=	
OBJS=		kthread.o bseq.o sketch.o index.o genidx.o
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
genidx.o: minimap.h bseq.h
index.o: minimap.h kvec.h khash.h ksort.h
sketch.o: kvec.h minimap.h
