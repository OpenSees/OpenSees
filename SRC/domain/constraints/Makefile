include ../../../Makefile.def

OBJS       = SP_Constraint.o MP_Constraint.o RigidDiaphragm.o \
	RigidRod.o RigidBeam.o ImposedMotionSP.o ImposedMotionSP1.o \
	Pressure_Constraint.o

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
