include ../../../Makefile.def

OBJS       = Truss.o \
	TrussSection.o \
	CorotTruss.o \
	CorotTrussSection.o \
	Truss2.o \
	N4BiaxialTruss.o \
	CorotTruss2.o \
	InertiaTruss.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
