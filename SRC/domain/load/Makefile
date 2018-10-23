include ../../../Makefile.def

OBJS       = Load.o \
	ElementalLoadIter.o \
	NodalLoadIter.o \
	Beam2dUniformLoad.o \
	Beam2dPointLoad.o \
	Beam2dPartialUniformLoad.o \
	Beam3dUniformLoad.o \
	Beam3dPointLoad.o \
	Beam3dPartialUniformLoad.o \
	BrickSelfWeight.o \
	Beam2dTempLoad.o \
	Beam2dThermalAction.o \
	NodalThermalAction.o \
	ThermalActionWrapper.o \
	Beam3dThermalAction.o \
	ShellThermalAction.o \
	SurfaceLoader.o \
	SelfWeight.o \
	LysmerVelocityLoader.o 

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
