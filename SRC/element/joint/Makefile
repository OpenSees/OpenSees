include ../../../Makefile.def

OBJS       = TclJoint2dCommand.o \
	Joint2D.o \
	MP_Joint2D.o \
	BeamColumnJoint2d.o \
	Joint3D.o \
	LehighJoint2d.o \
	MP_Joint3D.o \
	BeamColumnJoint3d.o \
	ElasticTubularJoint.o \
	TclBeamColumnJointCommand.o \
	TclJoint3dCommand.o 

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
