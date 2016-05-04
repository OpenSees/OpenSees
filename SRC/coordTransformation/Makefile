include ../../Makefile.def



OBJS       = CrdTransf.o \
	LinearCrdTransf2d.o \
	LinearCrdTransf3d.o \
	PDeltaCrdTransf2d.o \
	PDeltaCrdTransf3d.o \
	ModerateDispCrdTransf2d.o \
	CorotCrdTransf2d.o \
	CorotCrdTransf3d.o \
	CorotCrdTransfWarping2d.o \
	TclGeomTransfCommand.o

# Compilation control

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
