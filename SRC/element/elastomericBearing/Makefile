include ../../../Makefile.def

OBJS = ElastomericBearingPlasticity2d.o \
	ElastomericBearingPlasticity3d.o \
	ElastomericBearingBoucWen2d.o \
	ElastomericBearingBoucWen3d.o \
	ElastomericBearingBoucWenMod3d.o \
	ElastomericBearingUFRP2d.o \
	ElastomericX.o \
	HDR.o \
	LeadRubberX.o \
	TclElastomericBearingBoucWenCommand.o \
	TclElastomericBearingUFRPCommand.o \
	TclElastomericBearingPlasticityCommand.o

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
