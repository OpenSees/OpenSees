include ../../../../Makefile.def

OBJS       = NullYS2D.o Attalla2D.o  ElTawil2D.o  ElTawil2DUnSym.o  Hajjar2D.o \
	Orbison2D.o  YieldSurface_BC2D.o  YieldSurface_BC.o \
	TclModelBuilderYieldSurfaceBCCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
