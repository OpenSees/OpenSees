
include ../../Makefile.def

BINDIR = $(HOME)/bin
LIBS = -lm 

INCLUDES = -I. 
OPTIMIZATION = -O

DEBUGFLAGS = 
EXTRADEFINES =  $(DEBUGFLAGS) 
CFLAGS = $(OPTIMIZATION) $(EXTRADEFINES) $(INCLUDES)
LD = $(CC) -L. $(LIBDIR)

MLEVEL_OBJS = main.o io.o util.o match.o coarsen.o memory.o fm.o  \
              GKlib.o htable.o list.o refine.o bucketlist.o stat.o \
              mlevelpart.o separator.o mlevelorder.o smbfactor.o \
              smbfct.o mincover.o initpart.o lanczos.o mmd.o balpart.o \
              kwaypart.o kwayfm.o kwaygreedy.o 

lib:	
	$(MAKE) -f Makefile.lib

all: $(MLEVEL_OBJS)
	touch *.c
	make $(BINDIR)/metis
	touch *.c
	$(MAKE) -f Makefile.lib

$(BINDIR)/metis: $(MLEVEL_OBJS)
	$(LD) -o $@ $(MLEVEL_OBJS) $(LIBS)
	rm -f $(BINDIR)/pmetis $(BINDIR)/ometis $(BINDIR)/kmetis
	ln -s $@ $(BINDIR)/pmetis
	ln -s $@ $(BINDIR)/ometis
	ln -s $@ $(BINDIR)/kmetis


clean:
	@$(RM) $(RMFLAGS) *.o

realclean:
	$(RM) $(RMFLAGS) *.o 
	$(MAKE) -f Makefile.lib realclean

wipe: realclean

checkin:
	@for file in *.[c,h]; \
	  do \
	  ci -u -m'Maintance' $$file;\
	  done 

checkin2:
	@for file in *.[c,h]; \
	  do \
	  ci $$file;\
	  rcs -U $$file;\
	  co $$file;\
	  done 






