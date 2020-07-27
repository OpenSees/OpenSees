include ../../../../Makefile.def

ifdef H5DRM_FLAG
	H5_FILE = H5DRM.o
else
	H5_FIle = 
endif

OBJS       = Mesh3DSubdomain.o \
	GeometricBrickDecorator.o \
	DRMBoundaryLayerDecorator.o \
	DRMLoadPattern.o \
	DRMLoadPatternWrapper.o \
	DRMInputHandler.o \
	PlaneDRMInputHandler.o \
	$(H5_FILE)

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ 

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
