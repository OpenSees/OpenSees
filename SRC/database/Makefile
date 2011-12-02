#  Written: fmk 
#  Created: 10/97

include ../../Makefile.def

OBJS       = FE_Datastore.o FileDatastore.o TclDatabaseCommands.o MySqlDatastore.o \
	NEESData.o

# Compilation control

all:         $(OBJS)

test:  NeesDataTest.o 
	$(LINKER) $(LINKFLAGS) NeesDataTest.o NEESData.o $(FE_LIBRARY) \
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



