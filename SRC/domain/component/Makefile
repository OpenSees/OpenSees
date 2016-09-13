include ../../../Makefile.def

OBJS       = DomainComponent.o \
	Parameter.o \
	ElementParameter.o \
	RVParameter.o \
	DVParameter.o \
	NodeResponseParameter.o \
	LoadFactorParameter.o \
	TclParameterCommands.o \
	MaterialStageParameter.o \
	MatParameter.o \
	InitialStateParameter.o \
	ElementStateParameter.o \
	TclUpdateMaterialCommand.o

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
