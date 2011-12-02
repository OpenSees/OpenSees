include ../../Makefile.def

UMFD	= umd2fb.o umd2co.o umd2fa.o umd2f0.o umd2f1.o umd2f2.o \
	  umd2fg.o umd2in.o umd2of.o umd2s2.o umd2sl.o umd2so.o umd2su.o \
	  umd2er.o umd2p1.o umd2p2.o umd2lt.o umd2ut.o \
	  umd2rf.o umd2ra.o umd2r0.o umd2r2.o umd2rg.o umd21i.o

HARWELL	= mc21b.o mc13e.o


all: $(UMFD) $(HARWELL)
	@$(RM) $(RMFLAGS) $(UMFPACK_LIBRARY)
	@$(AR) $(ARFLAGS) $(UMFPACK_LIBRARY) $(UMFD) $(UMFS) $(HARWELL)
	@$(RANLIB) $(UMFPACK_LIBRARY)

clean:
	@$(RM) $(RMFLAGS) *.o

wipe: clean

