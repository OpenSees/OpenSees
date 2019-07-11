include ../../../../Makefile.def

OBJS       = SectionIntegration.o \
	WideFlangeSectionIntegration.o RCSectionIntegration.o \
	RCTBeamSectionIntegration.o RCCircularSectionIntegration.o \
	TubeSectionIntegration.o RCTunnelSectionIntegration.o

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
