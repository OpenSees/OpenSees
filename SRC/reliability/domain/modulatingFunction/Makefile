include ../../../../Makefile.def

OBJS       = ConstantModulatingFunction.o \
	GammaModulatingFunction.o \
	ModulatingFunction.o ModulatingFunctionIter.o \
	TrapezoidalModulatingFunction.o \
	KooModulatingFunction.o

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
