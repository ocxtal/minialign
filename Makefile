CC=			gcc
CFLAGS=		-g -O3 -Wall -Wc++-compat -Wno-unused-function
CPPFLAGS=
INCLUDES=	-I.
OBJS=		kthread.o misc.o bseq.o sketch.o sdust.o index.o map.o
PROG=		minialign
PROG_EXTRA=	sdust minialign-lite
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

extra:all $(PROG_EXTRA)

minialign:main.o libminimap.a libgaba.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminimap -lgaba $(LIBS)

minialign-lite:example.o libminimap.a libgaba.a
		$(CC) $(CFLAGS) $< -o $@ -L. -lminimap -lgaba $(LIBS)

libminimap.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)

libgaba.a:gaba.c gaba.h gaba_wrap.c
		$(CC) -g -O3 -Wall -Wno-unused-function -march=native -std=c99 -DMODEL=LINEAR gaba.c -c -o gaba_linear.o
		$(CC) -g -O3 -Wall -Wno-unused-function -march=native -std=c99 -DMODEL=AFFINE gaba.c -c -o gaba_affine.o
		$(CC) -g -O3 -Wall -Wno-unused-function -march=native -std=c99 gaba_wrap.c -c -o gaba_wrap.o
		$(AR) -csr $@ gaba_linear.o gaba_affine.o gaba_wrap.o

sdust:sdust.c kdq.h kvec.h kseq.h sdust.h
		$(CC) -D_SDUST_MAIN $(CFLAGS) $< -o $@ -lz

clean:
		rm -fr gmon.out *.o a.out $(PROG) $(PROG_EXTRA) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bseq.o: bseq.h kseq.h
example.o: minimap.h bseq.h kseq.h
index.o: minimap.h bseq.h kvec.h khash.h
main.o: minimap.h bseq.h
map.o: bseq.h kvec.h minimap.h sdust.h ksort.h
misc.o: minimap.h bseq.h ksort.h
sdust.o: kdq.h kvec.h sdust.h
sketch.o: kvec.h minimap.h bseq.h
