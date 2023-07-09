include ../../../../Makefile.def

OBJS       = FluidSolidPorousMaterial.o \
	MultiYieldSurface.o \
	MultiYieldSurfaceClay.o \
	PressureDependMultiYield.o \
	TclUpdateMaterialStageCommand.o \
	PressureIndependMultiYield.o \
	T2Vector.o \
	PressureDependMultiYield02.o \
	PressureDependMultiYield03.o

all:	$(OBJS)

# Miscellaneous

tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
