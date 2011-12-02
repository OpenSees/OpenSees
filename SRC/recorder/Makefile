include ../../Makefile.def

OBJS       = DatastoreRecorder.o MaxNodeDispRecorder.o \
	ElementRecorder.o TclRecorderCommands.o FilePlotter.o \
	NodeRecorder.o

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
