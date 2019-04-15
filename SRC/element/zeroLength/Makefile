include ../../../Makefile.def

OBJS       = ZeroLength.o \
	ZeroLengthSection.o \
	ZeroLengthContact2D.o \
	ZeroLengthContact3D.o \
	ZeroLengthND.o \
	ZeroLengthContactNTS2D.o \
	ZeroLengthInterface2D.o \
	ZeroLengthRocking.o \
	CoupledZeroLength.o \
	ZeroLengthImpact3D.o \
	ZeroLengthVG_HG.o \
	TclZeroLength.o


#


all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core example main.o

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
