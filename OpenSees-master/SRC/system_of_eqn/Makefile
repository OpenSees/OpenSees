# Makefile for fe objects

include ../../Makefile.def

OBJS       = 

all: $(OBJS)
	@$(CD) $(FE)/system_of_eqn/linearSOE; $(MAKE);
	@$(CD) $(FE)/system_of_eqn/eigenSOE; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/system_of_eqn/linearSOE; $(MAKE) wipe;
	@$(CD) $(FE)/system_of_eqn/eigenSOE; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
