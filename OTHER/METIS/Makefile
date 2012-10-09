include ../../Makefile.def


.c.o:
	@$(ECHO) Making $@ from $<
	$(CC) $(CFLAGS) -I. -c $< -o $@

LD = $(CC) -L. 


OBJS = coarsen.o fm.o initpart.o match.o ccgraph.o memory.o \
       pmetis.o pqueue.o refine.o util.o timing.o debug.o \
       bucketsort.o graph.o stat.o kmetis.o kwayrefine.o \
       kwayfm.o balance.o ometis.o srefine.o sfm.o separator.o \
       mincover.o mmd.o mesh.o meshpart.o frename.o fortran.o \
       myqsort.o compress.o parmetis.o estmem.o \
       mpmetis.o mcoarsen.o mmatch.o minitpart.o mbalance.o \
       mrefine.o mutil.o mfm.o mkmetis.o mkwayrefine.o mkwayfmh.o \
       mrefine2.o minitpart2.o mbalance2.o mfm2.o \
       kvmetis.o kwayvolrefine.o kwayvolfm.o subdomains.o


all: $(OBJS)
	@$(RM) $(RMFLAGS) $(METIS_LIBRARY)
	@$(AR) $(ARFLAGS) $(METIS_LIBRARY) $(OBJS)
	@$(RANLIB) $(RANLIBFLAGS) $(METIS_LIBRARY)

clean:
	rm -f *.o

realclean:
	rm -f *.o ; rm -f ../libmetis.a
