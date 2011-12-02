include ../../../Makefile.def

OBJS       = ArrayOfTaggedObjects.o ArrayOfTaggedObjectsIter.o \
	MapOfTaggedObjectsIter.o MapOfTaggedObjects.o

# Compilation control

all:         $(OBJS)

test: main.o
	$(LINKER) $(LINKFLAGS) main.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o test

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o test

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
