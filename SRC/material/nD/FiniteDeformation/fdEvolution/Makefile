include ../../../../../Makefile.def

OBJS       = fdEvolution_S.o \
	     fdEvolution_SLS.o \
	     fdEvolution_T.o \
	     fdEvolution_TL.o 

all: $(OBJS)

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS)  *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.

