include ../../Makefile.def



OBJS       = CrdTransf.o CrdTransf2d.o CrdTransf3d.o \
	LinearCrdTransf2d.o LinearCrdTransf3d.o \
	PDeltaCrdTransf2d.o PDeltaCrdTransf3d.o \
	TclGeomTransfCommand.o

# Compilation control

all:         $(OBJS)


lib: CorotCrdTransf2d.o CorotCrdTransf3d.o 
	@$(AR) $(ARFLAGS) libCorot.a  CorotCrdTransf2d.o CorotCrdTransf3d.o 
	@$(RANLIB) libCorot.a

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
