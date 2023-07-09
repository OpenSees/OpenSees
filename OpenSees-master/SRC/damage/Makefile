include ../../Makefile.def

OBJS       = DamageModel.o \
	HystereticEnergy.o \
	Kratzig.o \
	Mehanny.o \
	ParkAng.o \
	NormalizedPeak.o \
	DamageResponse.o \
	TclModelBuilderDamageModelCommand.o

# Compilation control

all:         $(OBJS)
	@$(CD) $(FE)/recorder/response; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
