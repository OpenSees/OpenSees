include ../../Makefile.def

OBJS       = Timer.o FileIter.o File.o SimulationInformation.o StringContainer.o PeerNGA.o

# Compilation control

all:         $(OBJS)


test: test.o $(OBJS)
	$(LINKER) $(LINKFLAGS) test.o $(OBJS) $(FE_LIBRARY) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(METIS_LIBRARY) \
	 -o test

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
