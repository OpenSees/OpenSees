include ../../../../Makefile.def

OBJS       = QuadRule.o QuadRule1d01.o GaussQuadRule1d01.o \
	GaussLobattoQuadRule1d01.o GaussQuadRule1d.o QuadRule1d.o


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
