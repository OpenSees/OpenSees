include ../../../../Makefile.def

OBJS = stressDensity.o SDM-UC.o

all: $(OBJS) 

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o *.mod

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
