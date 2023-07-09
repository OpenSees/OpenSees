include ../../../../../Makefile.def

OBJS       = LimitCurve.o \
	AxialCurve.o \
	ShearCurve.o \
	RotationShearCurve.o \
	ThreePointCurve.o \
	WrapperLimitCurve.o


all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
