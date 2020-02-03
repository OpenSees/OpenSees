include ../../../../Makefile.def

OBJS       = 	BetaRV.o \
		GumbelRV.o \
		ParetoRV.o \
		Type1LargestValueRV.o \
		UniformRV.o \
		ChiSquareRV.o \
		LaplaceRV.o \
		RayleighRV.o \
		Type1SmallestValueRV.o \
		WeibullRV.o \
		ExponentialRV.o \
		LognormalRV.o \
		ShiftedExponentialRV.o \
		Type2LargestValueRV.o \
		GammaRV.o \
		NormalRV.o \
		ShiftedRayleighRV.o \
		Type3SmallestValueRV.o \
		UserDefinedRV.o \
		PythonRV.o

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
