include ../../../Makefile.def

OBJS       = DispBeamColumn2dInt.o \
	FiberSection2dInt.o \
	LinearCrdTransf2dInt.o \
	TclDispBeamColumnIntCommand.o


all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
