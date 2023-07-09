include ../../Makefile.def



OBJS       = Damping.o \
	SecStifDamping.o \
	UniformDamping.o \
	URDDamping.o \
	URDDampingbeta.o \
	TclDampingCommand.o

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
