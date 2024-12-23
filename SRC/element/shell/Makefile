include ../../../Makefile.def

OBJS       = ShellMITC4.o \
	ShellMITC9.o \
	ShellDKGQ.o \
	ShellNLDKGQ.o \
	ShellDKGT.o \
	ShellNLDKGT.o \
	ShellMITC4Thermal.o \
	ShellNLDKGQThermal.o \
    ShellNLDKGTThermal.o \
	ASDShellQ4.o \
	ASDShellT3.o \
	R3vectors.o 

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
