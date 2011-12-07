include ../../../Makefile.def

OBJS       = SimpleContact2D.o \
             SimpleContact3D.o \
 		     BeamContact2D.o \
			 BeamContact2Dp.o \
			 BeamContact3D.o \
			 BeamEndContact3D.o \
             SSPquad.o \
	         SSPquadUP.o \
			 SSPbrick.o \
			 SSPbrickUP.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
