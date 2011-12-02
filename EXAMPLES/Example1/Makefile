include ../../Makefile.def

PROGRAM         = example1

all:         $(PROGRAM)

$(PROGRAM):  main.o
	$(LINKER) $(LINKFLAGS) main.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o $(PROGRAM)

rohit:  rohit.o
	$(LINKER) $(LINKFLAGS) rohit.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o rohit

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
