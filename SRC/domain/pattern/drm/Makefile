include ../../../../Makefile.def

OBJS       = Mesh3DSubdomain.o \
	GeometricBrickDecorator.o \
	DRMBoundaryLayerDecorator.o \
	DRMLoadPattern.o \
	DRMLoadPatternWrapper.o \
	DRMInputHandler.o \
	PlaneDRMInputHandler.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ 

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
