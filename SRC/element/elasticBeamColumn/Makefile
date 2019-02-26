include ../../../Makefile.def

OBJS       = 	ElasticBeam2d.o \
	ModElasticBeam2d.o \
	ElasticBeam3d.o \
	ElasticTimoshenkoBeam2d.o \
	ElasticTimoshenkoBeam3d.o \
	WheelRail.o \
	TclElasticBeamCommand.o

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
