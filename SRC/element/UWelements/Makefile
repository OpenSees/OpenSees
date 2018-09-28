include ../../../Makefile.def

OBJS       = SimpleContact2D.o \
             SimpleContact3D.o \
             BeamContact2D.o \
             BeamContact2Dp.o \
             BeamContact3D.o \
             BeamContact3Dp.o \
             BeamEndContact3D.o \
             BeamEndContact3Dp.o \
             SSPquad.o \
             SSPquadUP.o \
             SSPbrick.o \
             SSPbrickUP.o \
	     Quad4FiberOverlay.o \
	     Brick8FiberOverlay.o \
             QuadBeamEmbedContact.o \
	     PileToe3D.o \
			Tcl_generateInterfacePoints.o \
			EmbeddedBeamInterfaceL.o \
			EmbeddedBeamInterfaceP.o \
			EmbeddedEPBeamInterface.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
