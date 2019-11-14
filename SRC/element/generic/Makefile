include ../../../Makefile.def

OBJS       = GenericClient.o \
	GenericCopy.o \
	TclGenericClientCommand.o \
	TclGenericCopyCommand.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
