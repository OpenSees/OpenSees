include ../../../Makefile.def

OBJS       = FourNodeQuad.o \
	TclFourNodeQuadCommand.o \
	ConstantPressureVolumeQuad.o \
	EnhancedQuad.o \
	FourNodeQuad3d.o \
	NineNodeMixedQuad.o \
	NineNodeQuad.o \
	EightNodeQuad.o \
	SixNodeTri.o \
	FourNodeQuadWithSensitivity.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
