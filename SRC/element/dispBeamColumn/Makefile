include ../../../Makefile.def

OBJS       = DispBeamColumn2d.o DispBeamColumn3d.o \
	TclDispBeamColumnCommand.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
