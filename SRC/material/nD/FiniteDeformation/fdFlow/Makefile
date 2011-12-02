include ../../../../../Makefile.def

OBJS       = fdFlow.o \
	     fdFlowVM.o \
	     fdFlowDP.o \
	     fdFlowDP.o

all: $(OBJS)

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS)  *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.

