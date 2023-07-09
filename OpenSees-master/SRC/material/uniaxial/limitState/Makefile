include ../../../../Makefile.def

OBJS       = LimitStateMaterial.o \
	PinchingLimitStateMaterial.o \
	TclLimitState.o


all:         $(OBJS)
	@$(CD) $(FE)/material/uniaxial/limitState/limitCurve; $(MAKE);

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/material/uniaxial/limitState/limitCurve; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
