include ../../../../Makefile.def

OBJS       = FedeasHardeningMaterial.o Hard1.o \
	FedeasBond1Material.o \
	FedeasBond2Material.o \
	FedeasConcr1Material.o \
	FedeasConcr2Material.o \
	FedeasConcr3Material.o \
	FedeasHyster1Material.o \
	FedeasHyster2Material.o \
	FedeasSteel1Material.o \
	FedeasSteel2Material.o Steel2.o Concr2.o \
	TclFedeasMaterialCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
